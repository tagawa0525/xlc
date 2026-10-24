"""分岐のテーブルの交差の項を落とす誤差を測る（ADR-0200）。軸の直積で解いたテーブルから、公称の状態と、公称の状態から 1 軸ずつ
変えた分岐と、選んだ軸の組（交差の項）の部分格子だけを取り出して直積のすべての状態を組み立て（cut-HDMR）、直積の状態と比べる。

使い方: python3 compare_branch_cross_terms.py --json <要約の出力先> <入力（TOML）>
        python3 compare_branch_cross_terms.py --report <要約の JSON>
        python3 compare_branch_cross_terms.py --subgrids <直積のテーブル> <出力> <軸の組>...

`--subgrids` は、lattice が書いた直積のテーブル（軸の公称の値を持つ）から、公称の状態と、1 軸ずつの分岐と、軸の組（`boron,moderator`
のようにカンマで区切る）の部分格子の点だけを残したテーブル（ADR-0200 の形）を書く。燃焼度と制御棒の軸はすべての値を残す。core が
これから組み立てた状態を、ここでの組み立て（直積のテーブルから同じ点を取り出す）と比べる。

入力の TOML（パスは TOML のあるディレクトリから）:

    [[table]]
    name = "2A（40 MWd/kg）"
    file = "table-2A-cross.msgpack"
    rods = "table-2A-cross-rods.msgpack"   # 省略できる。同じ軸の、制御棒を入れた格子のテーブル
    nominal = { boron = 1300, moderator = 565, fuel = 565 }

    [[model]]
    name = "ホウ素 × 減速材"
    keep = [["boron", "moderator"]]        # 部分格子を解く軸の組。1 軸ずつの分岐はどのモデルも持つ

- 組み立て: 公称の状態を x0 とし、モデルの軸の組とその部分集合の全体を U（空集合と 1 軸ずつを含む）として、状態 x の値を
  Σ_{u∈U} Σ_{v⊆u} (−1)^{|u|−|v|} f(x_v, x0) とする（f(x_v, x0) は v の軸だけ x の値にし、ほかを公称にした状態。cut-HDMR。
  H. Rabitz, Ö. F. Aliş, J. Math. Chem. 25 (1999) 197）。U に入る軸の組の部分格子の点では、テーブルの値をそのまま返す
- 制御棒は、`rods` のテーブルを離散の軸 `rods`（公称は "out"、もう 1 点が "in"）として足す。モデルの組に `rods` があり
  テーブルに `rods` がなければ、そのモデルを飛ばす
- 燃焼度の軸があれば、燃焼度の点ごとにその点の公称の状態から組み立てる（燃焼度は組み立てる軸にしない）
- 状態の数: モデルが解く状態の数 Σ_{u∈U} Π_{a∈u} (n_a − 1)（n_a は軸 a の点の数）と、直積の状態の数 Π n_a
- 誤差: k∞ は少数群定数の無限媒質の固有値（compare_burnup_interpolation.py の `k_inf`）で、Δρ = (1/k_直積 − 1/k_組み立て) × 1e5 pcm。
  項目ごとの相対の差 |組み立て − 直積| / |直積| の、群（ADF は面と群）と状態についての最大も書く
- M² = D1 / (Σt1 − Σs1→1) + D2 / (Σt2 − Σs2→2)（公称の状態の 2 群の移動面積。出力の分布の誤差の見積もりに使う。ADR-0200）
- `inputs` は読んだファイルのパス（実行した cwd から）と md5
"""

import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_burnup_interpolation import interpolate, k_inf, unvector, vector

# 相対の差を書く項目
RELATIVE = ["fission", "nu_fission", "absorption", "diffusion", "adf"]
RODS = ["out", "in"]
# 分岐と公称の差を燃焼度について内挿する方式（compare_burnup_interpolation.py の名前。hermite は core の内挿。ADR-0192）
THINNING_METHODS = ["linear", "hermite"]


def closure(keep: list[list[str]], names: list[str]) -> list[tuple[str, ...]]:
    """モデルの軸の組と、その部分集合と、1 軸ずつ（空集合を含む）。軸は `names` の順"""
    sets = {()} | {(n,) for n in names}
    for k in keep:
        ordered = tuple(n for n in names if n in k)
        for r in range(len(ordered) + 1):
            sets.update(itertools.combinations(ordered, r))
    return sorted(sets, key=lambda s: (len(s), [names.index(n) for n in s]))


def rebuild(
    grid: dict, nominal: tuple, index: tuple, sets: list, names: list
) -> np.ndarray:
    """cut-HDMR の組み立て（モジュールの説明）。`grid` は軸の点の番号の組 → 状態の値の列"""
    total = np.zeros_like(grid[nominal])
    for u in sets:
        positions = [names.index(n) for n in u]
        for r in range(len(positions) + 1):
            for v in itertools.combinations(positions, r):
                at = tuple(
                    index[a] if a in v else nominal[a] for a in range(len(names))
                )
                total = total + (-1) ** (len(positions) - r) * grid[at]
    return total


def migration_area(state: dict) -> float:
    total = np.asarray(state["total"], dtype=float)
    scatter = np.asarray(state["scatter"], dtype=float)
    diffusion = np.asarray(state["diffusion"], dtype=float)
    if len(total) != 2:
        raise SystemExit(f"the migration area needs two groups, not {len(total)}")
    return float(sum(diffusion[g] / (total[g] - scatter[g][g]) for g in range(2)))


def relative_errors(rebuilt: dict, truth: dict) -> dict:
    out = {}
    for item in RELATIVE:
        t = np.asarray(truth[item], dtype=float)
        r = np.asarray(rebuilt[item], dtype=float)
        nonzero = t != 0
        if (r[~nonzero] != 0).any():
            raise SystemExit(f"{item}: the rebuilt state is not zero where the grid is")
        out[item] = (
            float(np.max(np.abs(r[nonzero] - t[nonzero]) / np.abs(t[nonzero])))
            if nonzero.any()
            else 0.0
        )
    return out


def table_grid(name: str, table: dict, nominal: dict) -> dict:
    """テーブルの燃焼度の点（なければ None）、連続の軸の名前と点、公称の状態の点の番号、燃焼度の点ごとの状態の数を確かめて返す"""
    if table.get("format") != "xlc-table" or table.get("split") != 1:
        raise SystemExit(
            f"{name}: the table must be an xlc-table that does not split the assembly"
        )
    axes = table["axes"]
    if any("points" not in a for a in axes):
        raise SystemExit(f"{name}: give the rods as a second table, not as a rod axis")
    burnup = None
    if axes[0]["name"] == "burnup":
        burnup = [float(b) for b in axes[0]["points"]]
        axes = axes[1:]
    names = [a["name"] for a in axes]
    points = [[float(p) for p in a["points"]] for a in axes]
    if sorted(nominal) != sorted(names):
        raise SystemExit(
            f"{name}: give the nominal value of every axis {names}, not {sorted(nominal)}"
        )
    for n, p in zip(names, points, strict=True):
        if float(nominal[n]) not in p:
            raise SystemExit(
                f"{name}: the nominal {n} {nominal[n]} is not a point of the axis {p}"
            )
    per_point = math.prod(len(p) for p in points)
    if len(table["states"]) != per_point * len(burnup or [None]):
        raise SystemExit(f"{name}: {len(table['states'])} states for the axes")
    return {
        "burnup": burnup,
        "names": names,
        "points": points,
        "nominal": tuple(
            p.index(float(nominal[n])) for n, p in zip(names, points, strict=True)
        ),
        "per_point": per_point,
    }


def summarize_table(
    name: str, table: dict, rodded: dict | None, nominal: dict, models: list[dict]
) -> dict:
    grid_of = table_grid(name, table, nominal)
    if rodded is not None:
        table_grid(f"{name} (rods)", rodded, nominal)
        if [(a["name"], a["points"]) for a in rodded["axes"]] != [
            (a["name"], a["points"]) for a in table["axes"]
        ]:
            raise SystemExit(
                f"{name}: the rodded table must have the axes of the table"
            )
    burnup, names, points = grid_of["burnup"], grid_of["names"], grid_of["points"]
    nominal_index = list(grid_of["nominal"])
    all_names = names + (["rods"] if rodded is not None else [])
    values = points + ([RODS] if rodded is not None else [])
    nominal_index = tuple(nominal_index + ([0] if rodded is not None else []))
    for m in models:
        for k in m["keep"]:
            unknown = [a for a in k if a not in names and a != "rods"]
            if unknown:
                raise SystemExit(
                    f"model {m['name']}: unknown axes {unknown}; the table has {names}"
                )
    used = [
        m
        for m in models
        if rodded is not None or all("rods" not in k for k in m["keep"])
    ]
    sizes = [len(v) for v in values]
    per_point = grid_of["per_point"]
    burnups = burnup if burnup is not None else [None]
    table_k = []
    out_points = []
    for b in range(len(burnups)):
        sources = [table] + ([rodded] if rodded is not None else [])
        states = {}
        for r, source in enumerate(sources):
            chunk = source["states"][b * per_point : (b + 1) * per_point]
            for flat, index in enumerate(
                itertools.product(*[range(len(p)) for p in points])
            ):
                states[index + ((r,) if rodded is not None else ())] = chunk[flat]
        grid = {i: vector(s) for i, s in states.items()}
        truth_k = {i: k_inf(s) for i, s in states.items()}
        table_k += [
            abs(1 / s["k_inf"] - 1 / truth_k[i]) * 1e5 for i, s in states.items()
        ]
        like = states[nominal_index]
        results = {}
        for m in used:
            sets = closure(m["keep"], all_names)
            worst = (0.0, nominal_index)
            squares = []
            relative = dict.fromkeys(RELATIVE, 0.0)
            for index, s in states.items():
                v = unvector(rebuild(grid, nominal_index, index, sets, all_names), like)
                rho = (1 / truth_k[index] - 1 / k_inf(v)) * 1e5
                squares.append(rho * rho)
                if abs(rho) > abs(worst[0]):
                    worst = (rho, index)
                for item, e in relative_errors(v, s).items():
                    relative[item] = max(relative[item], e)
            results[m["name"]] = {
                "keep": m["keep"],
                "states": sum(
                    math.prod(sizes[all_names.index(a)] - 1 for a in u) for u in sets
                ),
                "max_abs_pcm": abs(worst[0]),
                "worst_pcm": worst[0],
                "rms_pcm": float(np.sqrt(np.mean(squares))),
                "worst": {a: values[i][worst[1][i]] for i, a in enumerate(all_names)},
                "max_relative": relative,
            }
        out_points.append(
            {
                "burnup": burnups[b],
                "nominal_k": truth_k[nominal_index],
                "migration_area_cm2": migration_area(like),
                "models": results,
            }
        )
    return {
        "name": name,
        "axes": {a: v for a, v in zip(all_names, values, strict=True)},
        "nominal": {a: values[i][nominal_index[i]] for i, a in enumerate(all_names)},
        "grid_states": math.prod(sizes),
        "table_k_max_abs_pcm": max(table_k),
        "points": out_points,
    }


def thin(name: str, table: dict, nominal: dict, selections: list[dict]) -> dict:
    """分岐を燃焼度の点 `burnup` だけで解いたとして、公称でない各状態の、公称の状態との差（少数群定数の項目）を燃焼度について
    内挿し、選ばなかった点の公称の状態に足して、全点のテーブルの状態と比べる（Δρ は `summarize_table` と同じ向き）"""
    grid_of = table_grid(name, table, nominal)
    burnup, names, points = grid_of["burnup"], grid_of["names"], grid_of["points"]
    if burnup is None:
        raise SystemExit(f"{name}: the table has no burnup axis")
    indices = list(itertools.product(*[range(len(p)) for p in points]))
    per_point = grid_of["per_point"]
    at_nominal = indices.index(grid_of["nominal"])
    states = [
        table["states"][b * per_point : (b + 1) * per_point] for b in range(len(burnup))
    ]
    nominal_values = np.array([vector(s[at_nominal]) for s in states])
    out = []
    for sel in selections:
        chosen = [float(b) for b in sel["burnup"]]
        if (
            len(chosen) < 2
            or len(chosen) == len(burnup)
            or chosen != sorted(set(chosen))
            or any(b not in burnup for b in chosen)
            or chosen[0] != burnup[0]
            or chosen[-1] != burnup[-1]
        ):
            raise SystemExit(
                f"{name} {sel['name']}: the burnups {chosen} must be increasing points of the table from the first to the last, leaving some out"
            )
        index = [burnup.index(b) for b in chosen]
        removed = [i for i in range(len(burnup)) if i not in index]
        results = {}
        for method in THINNING_METHODS:
            worst_per_point = [0.0] * len(removed)
            worst = (0.0, None)
            for flat, grid_index in enumerate(indices):
                if flat == at_nominal:
                    continue
                difference = (
                    np.array([vector(s[flat]) for s in states]) - nominal_values
                )
                values = interpolate(
                    method,
                    np.array(chosen),
                    difference[index],
                    np.array([burnup[i] for i in removed]),
                )
                for n, i in enumerate(removed):
                    rebuilt = unvector(nominal_values[i] + values[n], states[i][flat])
                    rho = (1 / k_inf(states[i][flat]) - 1 / k_inf(rebuilt)) * 1e5
                    if abs(rho) > abs(worst_per_point[n]):
                        worst_per_point[n] = rho
                    if worst[1] is None or abs(rho) > abs(worst[0]):
                        worst = (
                            rho,
                            {"burnup": burnup[i]}
                            | {
                                a: points[k][grid_index[k]] for k, a in enumerate(names)
                            },
                        )
            results[method] = {
                "rho_pcm": worst_per_point,
                "max_abs_pcm": abs(worst[0]),
                "worst_pcm": worst[0],
                "worst": worst[1],
            }
        out.append(
            {
                "name": sel["name"],
                "burnup": chosen,
                "removed": [burnup[i] for i in removed],
                "results": results,
            }
        )
    return {"name": name, "burnup": burnup, "selections": out}


STATE_AXES = ["boron", "moderator", "fuel"]


def subgrid_table(table: dict, keep: list[list[str]]) -> dict:
    """直積のテーブル `table` から、公称の状態と 1 軸ずつの分岐と、`keep` の軸の組の部分格子の点だけを残したテーブル（モジュールの
    説明の `--subgrids`）"""
    axes = table["axes"]
    positions = [i for i, a in enumerate(axes) if a["name"] in STATE_AXES]
    if any("nominal" not in axes[i] for i in positions):
        raise SystemExit(
            "the table needs the nominal values of the boron, moderator and fuel axes"
        )
    unknown = [a for k in keep for a in k if a not in STATE_AXES]
    if unknown:
        raise SystemExit(f"unknown axes {unknown}; subgrids take {STATE_AXES}")
    sizes = [len(a["points"]) if "points" in a else len(a["values"]) for a in axes]
    product = [list(p) for p in itertools.product(*[range(n) for n in sizes])]
    if table.get("points") != product:
        raise SystemExit("the table must hold the product of its axes")
    nominal = {i: axes[i]["points"].index(axes[i]["nominal"]) for i in positions}
    kept = [{axes[i]["name"] for i in positions if axes[i]["name"] in k} for k in keep]

    def solved(point: list[int]) -> bool:
        off = {axes[i]["name"] for i in positions if point[i] != nominal[i]}
        # 1 軸ずつの分岐はどのモデルも持つ（`closure`）
        return len(off) <= 1 or any(off <= k for k in kept)

    sub = table["split"] ** 2
    points, states = [], []
    for n, point in enumerate(product):
        if solved(point):
            points.append(point)
            states += table["states"][n * sub : (n + 1) * sub]
    return {**table, "points": points, "states": states}


def summarize(tables: list[tuple], models: list[dict]) -> dict:
    """`tables` は (名前, テーブル, 制御棒を入れたテーブルか None, 公称の状態) の列"""
    return {
        "tables": [
            summarize_table(n, t, r, nominal, models) for n, t, r, nominal in tables
        ]
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    for t in summary["tables"]:
        axes = " × ".join(f"{a} {v}" for a, v in t["axes"].items())
        lines.append(
            f"{t['name']}: {axes}（直積 {t['grid_states']} 状態）、公称 {t['nominal']}、"
            f"k∞ の解き直しとテーブルの k_inf の差の最大 {t['table_k_max_abs_pcm']:.4f} pcm"
        )
        for p in t["points"]:
            where = "" if p["burnup"] is None else f"{p['burnup']:g} MWd/kg: "
            lines.append(
                f"  {where}公称の k∞ {p['nominal_k']:.6f}、M² {p['migration_area_cm2']:.2f} cm²"
            )
            for model, r in p["models"].items():
                rel = "、".join(
                    f"{k} {v * 100:.3f}%" for k, v in r["max_relative"].items()
                )
                lines.append(
                    f"    {model}（{r['states']} 状態）: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_pcm']:+.2f}、"
                    f"{r['worst']}）、rms {r['rms_pcm']:.2f} pcm、相対の差の最大 {rel}"
                )
    for t in summary.get("thinning", []):
        lines.append(f"{t['name']}: 燃焼度の点 {len(t['burnup'])}")
        for sel in t["selections"]:
            lines.append(f"  {sel['name']}（{len(sel['burnup'])} 点）: {sel['burnup']}")
            for method, r in sel["results"].items():
                lines.append(
                    f"    {method}: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_pcm']:+.2f}、{r['worst']}）"
                )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 3 and args[0] == "--subgrids":
        keep = [k.split(",") for k in args[3:]]
        out = subgrid_table(msgpack.unpackb(Path(args[1]).read_bytes()), keep)
        Path(args[2]).write_bytes(msgpack.packb(out))
    elif len(args) == 3 and args[0] == "--json":
        config_path = Path(args[2])
        config = tomllib.loads(config_path.read_text())
        files = [config_path]
        loaded = {}

        def read(path: Path) -> dict:
            """テーブルを 1 度だけ読み、`inputs` に 1 度だけ載せる"""
            if path not in loaded:
                files.append(path)
                loaded[path] = msgpack.unpackb(path.read_bytes())
            return loaded[path]

        tables = []
        for t in config["table"]:
            rodded = None
            if "rods" in t:
                rodded = read(config_path.parent / t["rods"])
            tables.append(
                (t["name"], read(config_path.parent / t["file"]), rodded, t["nominal"])
            )
        thinning = []
        for t in config.get("thinning", []):
            thinning.append(
                thin(
                    t["name"],
                    read(config_path.parent / t["file"]),
                    t["nominal"],
                    t["selections"],
                )
            )
        summary = {
            "format": "xlc-compare-branch-cross-terms",
            "thinning": thinning,
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(tables, config["model"]),
        }
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
