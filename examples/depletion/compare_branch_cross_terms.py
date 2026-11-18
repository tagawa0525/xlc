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

    [[thinning]]                           # 省略できる。燃焼度の軸を持ち、すべての状態をすべての燃焼点で解いたテーブル（直積か
    name = "2A"                            # ADR-0200 の部分格子）で、分岐を解く燃焼点を間引く
    file = "table-2A.msgpack"
    nominal = { boron = 650, moderator = 580, fuel = 900 }
    selections = [{ name = "3 点", burnup = [0, 0.5, 60] }]   # 最初と最後の点を含み、どれかの点を除く

    [thinning.timing]                      # 省略できる。テーブルを作る時間の見込み（`thinning_time`）
    total = 2349                           # すべての状態を解いた燃焼点 1 つの秒（`other` と一方だけ）
    solve = 36.1                           # 自己遮蔽を解く 1 状態の秒
    reuse = 8.7                            # 自己遮蔽を再利用する 1 状態の秒
    reuse-axis = "boron"
    # other = 14.6                         # `total` の代わりに、燃焼点 1 つの状態によらない秒
    # states = { boron = 5, ... }          # 省略できる。見込むテーブルの種類ごとの状態の数（どの種類も書く。`other` と使う）。省略するとこのテーブルの数
    # nominal-burnup = [0, 0.1, ...]       # 省略できる。見込むテーブルの公称の状態の燃焼点（このテーブルの点から。選び方の点を含む。
    #                                      # 省略するとこのテーブルの点）

    [[thinning.plans]]                     # 省略できる（timing が要る）。状態の種類ごとに選び方を決めた組
    name = "混ぜる"
    kinds = { boron = "3 点", ... }        # すべての種類に、selections の名前を書く

    [[nominal]]                            # 省略できる。`[[thinning]]` と同じ形のテーブルで、公称の状態を解く燃焼点を間引く
    name = "2O"
    file = "table-2O-all.msgpack"
    nominal = { boron = 650, moderator = 580, fuel = 900 }
    selections = [{ name = "27 点", burnup = [0, 0.1, ...] }]

    [[combined]]                           # 省略できる。`[[thinning]]` と同じ形のテーブルで、公称の状態と公称でない状態の燃焼点を
    name = "2O"                            # 組にして間引く
    file = "table-2O-all.msgpack"
    nominal = { boron = 650, moderator = 580, fuel = 900 }
    selections = [{ name = "33 点", burnup = [...] }, { name = "19 点", burnup = [...] }]
    pairs = [{ nominal = "33 点", off-nominal = "19 点" }]   # 公称でない点は公称の点の中から

    [combined.timing]                      # 省略できる。`[thinning.timing]` と同じ（`nominal-burnup` は書かない）

`[[table]]` と `[[model]]` も省略できる（間引きだけを比べる入力）。

- 間引き: 選んだ燃焼点の、公称でない各状態の公称の状態との差を燃焼度について内挿し（`THINNING_METHODS`）、除いた点で全点の
  テーブルの状態と比べる。最大は、すべての状態と、状態の種類（公称と違う軸の組。`boron`、`boron+moderator` など）ごとに書く。
  テーブルが微視的断面積を持てば、微視的断面積と数密度の微分も同じ内挿で比べる（core の `interpolate_micro` は少数群定数と
  同じ重みで足す。差の測り方は `micro_errors`）

- 公称の状態の燃焼点: 選んだ燃焼点の公称の状態を燃焼度について内挿し（`THINNING_METHODS`）、除いた点で全点のテーブルの公称の
  状態と比べる（`nominal_points`）。間引きは公称の状態をテーブルのすべての点の値で持つので、この誤差を含まない（時間の見込みの
  `nominal-burnup` の外の点でも、公称の状態は内挿でなくテーブルの値）

- 公称の状態と公称でない状態の燃焼点の組: 公称の状態を組の公称の点から、公称でない状態の公称の状態との差を組の公称でない点から
  内挿して足し（core の組み立て。`combine`）、公称でない点から除いた点で、全点のテーブルのすべての状態と比べる

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
    # 状態を軸の直積の順の番号で読むので、解いた点が直積の点を同じ順に持つこと（ADR-0200 の部分格子のテーブルは読まない）
    sizes = [len(burnup)] if burnup is not None else []
    sizes += [len(p) for p in points]
    product = [list(p) for p in itertools.product(*[range(n) for n in sizes])]
    if table.get("points") != product:
        raise SystemExit(f"{name}: the table must hold the product of its axes")
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


def state_kind(names: list[str], index: tuple, nominal: tuple) -> str:
    """公称と違う軸の名前を軸の順に `+` でつないだもの（状態の種類。`boron+moderator` など）"""
    return "+".join(n for k, n in enumerate(names) if index[k] != nominal[k])


def kind_order(kind: str, names: list[str]) -> tuple:
    """変えた軸の数、軸の順の番号の順（1 軸ずつ、2 軸、3 軸）"""
    axes = kind.split("+")
    return (len(axes), [names.index(a) for a in axes])


def burnup_states(name: str, table: dict, nominal: dict) -> dict:
    """燃焼度を最初の軸に持つテーブル（直積か、ADR-0200 の部分格子）の、燃焼度の点、状態の軸の名前と点、公称の状態の点の番号と、
    状態の点の番号（燃焼度を除く。テーブルの順）ごとの、燃焼度の点ごとの `states` の位置。間引きの真値にするので、どの状態も
    すべての燃焼度の点で解いてあること"""
    if table.get("format") != "xlc-table" or table.get("split") != 1:
        raise SystemExit(
            f"{name}: the table must be an xlc-table that does not split the assembly"
        )
    axes = table["axes"]
    if any("points" not in a for a in axes):
        raise SystemExit(f"{name}: give the rods as a second table, not as a rod axis")
    if axes[0]["name"] != "burnup":
        raise SystemExit(f"{name}: the table has no burnup axis")
    burnup = [float(b) for b in axes[0]["points"]]
    names = [a["name"] for a in axes[1:]]
    points = [[float(p) for p in a["points"]] for a in axes[1:]]
    if sorted(nominal) != sorted(names):
        raise SystemExit(
            f"{name}: give the nominal value of every axis {names}, not {sorted(nominal)}"
        )
    for a, p in zip(axes[1:], points, strict=True):
        given = float(nominal[a["name"]])
        if given not in p:
            raise SystemExit(
                f"{name}: the nominal {a['name']} {given} is not a point of the axis {p}"
            )
        if "nominal" in a and float(a["nominal"]) != given:
            raise SystemExit(
                f"{name}: the table's nominal {a['name']} {a['nominal']} is not the given {given}"
            )
    if len(table["states"]) != len(table["points"]):
        raise SystemExit(
            f"{name}: {len(table['states'])} states for {len(table['points'])} points"
        )
    rows: dict[tuple, list] = {}
    for position, point in enumerate(table["points"]):
        rows.setdefault(tuple(point[1:]), [None] * len(burnup))[point[0]] = position
    if any(None in r for r in rows.values()):
        raise SystemExit(f"{name}: solve every state at every burnup point")
    at_nominal = tuple(
        p.index(float(nominal[n])) for n, p in zip(names, points, strict=True)
    )
    if at_nominal not in rows:
        raise SystemExit(f"{name}: the table has no nominal state {nominal}")
    return {
        "burnup": burnup,
        "names": names,
        "points": points,
        "nominal": at_nominal,
        "rows": rows,
    }


# 微視的断面積（ADR-0216）と数密度の微分（ADR-0221）の、比べる項目
MICRO_ITEMS = [
    "total",
    "nu_fission",
    "fission",
    "scatter",
    "reactions",
    "thermal_fission",
]
DERIVATIVE_ITEMS = [
    "total",
    "nu_fission",
    "fission",
    "scatter",
    "rate_fission",
    "rate_reactions",
    "rate_thermal_fission",
]


def micro_parts(state: list[dict], nuclides: list[str]) -> list[tuple]:
    """状態の微視的断面積を、(`micro` か `derivative`, 項目, 核種（反応は mt を添える）, 値の列) の区切りに並べる。境界より下の
    核分裂は、割合でなく断面積（割合 × 群の和の σf）にする（core の `interpolate_micro` と同じく、その形で足して割合に戻す）"""
    out = []

    def put(kind, item, label, values):
        out.append((kind, item, label, np.ravel(np.asarray(values, dtype=float))))

    if len(state) != len(nuclides):
        raise SystemExit(f"{len(state)} micro nuclides for {len(nuclides)} names")
    for name, n in zip(nuclides, state, strict=True):
        for item in ["total", "nu_fission", "fission", "scatter"]:
            put("micro", item, name, n[item])
        for r in n["reactions"]:
            put("micro", "reactions", f"{name} mt{r['mt']}", r["channels"])
        if n.get("thermal_fission") is not None:
            put(
                "micro",
                "thermal_fission",
                name,
                [n["thermal_fission"] * sum(n["fission"])],
            )
        d = n.get("derivative")
        if d is None:
            continue
        for item in ["total", "nu_fission", "fission", "scatter"]:
            put("derivative", item, name, d[item])
        rate = d.get("rate")
        if rate is not None:
            put("derivative", "rate_fission", name, rate["fission"])
            for r in rate["reactions"]:
                put(
                    "derivative", "rate_reactions", f"{name} mt{r['mt']}", r["channels"]
                )
            put("derivative", "rate_thermal_fission", name, [rate["thermal_fission"]])
    return out


class MicroLayout:
    """`micro_parts` の区切りの並び（どの状態も同じ）と、並べた値から比べる値への変換"""

    def __init__(self, parts: list[tuple]):
        self.labels = [(kind, item, label) for kind, item, label, _ in parts]
        self.bounds = []
        start = 0
        for *_, values in parts:
            self.bounds.append((start, start + len(values)))
            start += len(values)
        self.size = start
        # 境界より下の核分裂の区切りから、同じ核種の σf の区切りへ
        self.fission_of = {}
        for k, (kind, item, label) in enumerate(self.labels):
            if kind == "micro" and item == "thermal_fission":
                self.fission_of[k] = self.labels.index(("micro", "fission", label))

    def vector(self, state: list[dict], nuclides: list[str]) -> np.ndarray:
        parts = micro_parts(state, nuclides)
        if [(kind, item, label) for kind, item, label, _ in parts] != self.labels or (
            sum(len(p[3]) for p in parts) != self.size
        ):
            raise SystemExit("the micro data differ in their items between states")
        return np.concatenate([p[3] for p in parts])

    def compared(self, v: np.ndarray) -> np.ndarray:
        """境界より下の核分裂を割合に戻した値（`v` は [点][値]）"""
        out = v.copy()
        for k, f in self.fission_of.items():
            a, _ = self.bounds[k]
            fa, fb = self.bounds[f]
            fission = v[:, fa:fb].sum(axis=1)
            out[:, a] = np.divide(
                v[:, a], fission, out=np.zeros_like(fission), where=fission != 0
            )
        return out


def micro_errors(
    layout: MicroLayout, rebuilt: np.ndarray, truth: np.ndarray
) -> np.ndarray:
    """点（行）と値ごとの差。微視的断面積は |組み立て − 真値| / |真値|、数密度の微分は区切りの真値の大きさの最大で割る（微分は
    群や反応の間で 0 を横切り、要素ごとの比は 0 の近くで大きくなるだけで量の誤差を表さない）。真値が 0 の値（微分は区切りの
    すべてが 0）は、組み立ても 0 でなければエラー"""
    r = layout.compared(rebuilt)
    t = layout.compared(truth)
    scale = np.abs(t)
    for (kind, *_), (a, b) in zip(layout.labels, layout.bounds, strict=True):
        if kind == "derivative":
            scale[:, a:b] = np.abs(t[:, a:b]).max(axis=1, keepdims=True)
    difference = np.abs(r - t)
    zero = scale == 0
    if (difference[zero] != 0).any():
        column = int(np.argwhere(zero & (difference != 0))[0][1])
        k = next(n for n, (a, b) in enumerate(layout.bounds) if a <= column < b)
        raise SystemExit(
            f"{layout.labels[k]}: the interpolation is not zero where the table is"
        )
    return np.divide(difference, scale, out=np.zeros_like(difference), where=~zero)


def item_worst(
    layout: MicroLayout, errors: np.ndarray, at: list[dict], into: dict
) -> None:
    """`errors`（[点][値]。点の位置は `at`）の項目ごとの最大を、`into`（kind → 項目 → 最大）に足す"""
    for k, ((kind, item, label), (a, b)) in enumerate(
        zip(layout.labels, layout.bounds, strict=True)
    ):
        if a == b:
            continue
        block = errors[:, a:b]
        n = int(np.argmax(block.max(axis=1)))
        e = float(block[n].max())
        worst = into.setdefault(kind, {}).setdefault(
            item, {"max_relative": 0.0, "nuclide": label, "worst": at[n]}
        )
        if e > worst["max_relative"]:
            worst.update({"max_relative": e, "nuclide": label, "worst": at[n]})


def ordered_items(found: dict) -> dict:
    """`item_worst` の結果を、`MICRO_ITEMS` と `DERIVATIVE_ITEMS` の順に並べる"""
    out = {}
    for kind, items in [("micro", MICRO_ITEMS), ("derivative", DERIVATIVE_ITEMS)]:
        if kind in found:
            out[kind] = {i: found[kind][i] for i in items if i in found[kind]}
    return out


def thin(
    name: str,
    table: dict,
    nominal: dict,
    selections: list[dict],
    timing: dict | None = None,
    plans: list[dict] = (),
) -> dict:
    """分岐を燃焼度の点 `burnup` だけで解いたとして、公称でない各状態の、公称の状態との差（少数群定数の項目と、テーブルが持てば
    微視的断面積と数密度の微分）を燃焼度について内挿し、選ばなかった点の公称の状態に足して、全点のテーブルの状態と比べる（Δρ は
    `summarize_table` と同じ向き。微視的断面積の差は `micro_errors`）。最大は、すべての状態についてのものと、状態の種類
    （`state_kind`）ごとのもの（`kinds`）を書く。テーブルは直積でも ADR-0200 の部分格子でもよい（`burnup_states`）。

    `timing`（`thinning_time`）があれば、選び方ごとにテーブルを作る時間の見込みを書く。`plans` は状態の種類ごとに選び方を
    決めた組で、誤差は種類ごとにその選び方の値を取る（状態ごとの内挿はほかの状態の燃焼点によらない）"""
    grid_of = burnup_states(name, table, nominal)
    burnup, names, points = grid_of["burnup"], grid_of["names"], grid_of["points"]
    rows = grid_of["rows"]
    at_nominal = grid_of["nominal"]
    states = table["states"]
    nominal_values = np.array([vector(states[p]) for p in rows[at_nominal]])
    micro = table.get("micro")
    if micro is not None:
        if len(micro["states"]) != len(states):
            raise SystemExit(
                f"{name}: {len(micro['states'])} micro states for {len(states)} states"
            )
        nuclides = micro["nuclides"]
        layout = MicroLayout(
            micro_parts(micro["states"][rows[at_nominal][0]], nuclides)
        )

        def micro_rows(index):
            return np.array(
                [layout.vector(micro["states"][p], nuclides) for p in rows[index]]
            )

        nominal_micro = micro_rows(at_nominal)
    kinds_of = {
        index: state_kind(names, index, at_nominal)
        for index in rows
        if index != at_nominal
    }
    order = sorted(set(kinds_of.values()), key=lambda k: kind_order(k, names))
    counts = {k: sum(1 for v in kinds_of.values() if v == k) for k in order}
    seconds = None if timing is None else thinning_time(timing, counts)
    for sel in selections:
        chosen_points(name, sel, burnup)
    nominal_points_of = nominal_burnup(name, timing, burnup, selections)
    burnups = len(nominal_points_of)
    out = []
    for sel in selections:
        chosen, index, removed = chosen_points(name, sel, burnup)
        results = {}
        for method in THINNING_METHODS:
            worst_per_point = [0.0] * len(removed)
            worst = (0.0, None)
            per_kind = {}
            micro_kind = {}
            for grid_index, kind in kinds_of.items():
                kind_worst = per_kind.setdefault(kind, [0.0, None])
                positions = rows[grid_index]
                difference = (
                    np.array([vector(states[p]) for p in positions]) - nominal_values
                )
                values = interpolate(
                    method,
                    np.array(chosen),
                    difference[index],
                    np.array([burnup[i] for i in removed]),
                )
                at = [
                    {"burnup": burnup[i]}
                    | {a: points[k][grid_index[k]] for k, a in enumerate(names)}
                    for i in removed
                ]
                for n, i in enumerate(removed):
                    truth = states[positions[i]]
                    rebuilt = unvector(nominal_values[i] + values[n], truth)
                    rho = (1 / k_inf(truth) - 1 / k_inf(rebuilt)) * 1e5
                    if abs(rho) > abs(worst_per_point[n]):
                        worst_per_point[n] = rho
                    if worst[1] is None or abs(rho) > abs(worst[0]):
                        worst = (rho, at[n])
                    if kind_worst[1] is None or abs(rho) > abs(kind_worst[0]):
                        kind_worst[:] = [rho, at[n]]
                if micro is not None:
                    truth = micro_rows(grid_index)
                    rebuilt = nominal_micro[removed] + interpolate(
                        method,
                        np.array(chosen),
                        (truth - nominal_micro)[index],
                        np.array([burnup[i] for i in removed]),
                    )
                    item_worst(
                        layout,
                        micro_errors(layout, rebuilt, truth[removed]),
                        at,
                        micro_kind.setdefault(kind, {}),
                    )
            kinds = {}
            for k in order:
                kinds[k] = {
                    "max_abs_pcm": abs(per_kind[k][0]),
                    "worst_pcm": per_kind[k][0],
                    "worst": per_kind[k][1],
                    "states": counts[k],
                }
                if seconds is not None:
                    kinds[k]["seconds_per_point"] = seconds["kinds"][k]
                    kinds[k]["timed_states"] = seconds["counts"][k]
                if micro is not None:
                    kinds[k] |= ordered_items(micro_kind[k])
            results[method] = {
                "rho_pcm": worst_per_point,
                "max_abs_pcm": abs(worst[0]),
                "worst_pcm": worst[0],
                "worst": worst[1],
                "kinds": kinds,
            }
            if micro is not None:
                overall = {}
                for k in order:
                    for kind, items in micro_kind[k].items():
                        for item, e in items.items():
                            best = overall.setdefault(kind, {}).setdefault(item, e)
                            if e["max_relative"] > best["max_relative"]:
                                overall[kind][item] = e
                results[method] |= ordered_items(overall)
        entry = {
            "name": sel["name"],
            "burnup": chosen,
            "removed": [burnup[i] for i in removed],
            "results": results,
        }
        if seconds is not None:
            entry["hours"] = table_hours(
                seconds, burnups, {k: len(chosen) for k in order}
            )
        out.append(entry)
    summary = {
        "name": name,
        "burnup": burnup,
        "nominal_burnup": nominal_points_of,
        "selections": out,
    }
    if seconds is not None:
        summary["hours"] = burnups * seconds["total"] / 3600
        summary["other_seconds"] = seconds["other"]
        summary["timing"] = dict(timing)
    if plans:
        if seconds is None:
            raise SystemExit(f"{name}: plans need the timing")
        summary["plans"] = [
            thinning_plan(name, p, out, order, timing["reuse-axis"], seconds, burnups)
            for p in plans
        ]
    return summary


def chosen_points(name: str, sel: dict, burnup: list[float]) -> tuple:
    """選び方 `sel` の燃焼点（テーブルの点から、最初から最後まで、どれかを除いて昇順に選んだもの）と、その点と除いた点の番号"""
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
    return chosen, index, removed


def nominal_burnup(
    name: str, timing: dict | None, burnup: list[float], selections: list[dict]
) -> list[float]:
    """時間の見込みで、公称の状態と状態によらない秒を数える燃焼点。`timing` の `nominal-burnup`（見込むテーブルの公称の状態の
    燃焼点。間引きの真値のテーブルより少ない点で公称の状態を解くテーブルの見込み）か、このテーブルの点。公称でない状態は公称の
    状態を解く燃焼点でしか解けない（公称の状態との差を持つ）ので、どの選び方の点もその中にあること"""
    if timing is None or "nominal-burnup" not in timing:
        return burnup
    given = timing["nominal-burnup"]
    points = (
        [float(b) for b in given]
        if isinstance(given, list)
        and all(isinstance(b, (int, float)) and not isinstance(b, bool) for b in given)
        else None
    )
    if (
        points is None
        or len(points) < 2
        or points != sorted(set(points))
        or any(b not in burnup for b in points)
        or points[0] != burnup[0]
        or points[-1] != burnup[-1]
    ):
        raise SystemExit(
            f"{name}: the timing nominal-burnup {given!r} must be increasing points of the table from the first to the last"
        )
    for sel in selections:
        outside = sorted({float(b) for b in sel["burnup"]} - set(points))
        if outside:
            raise SystemExit(
                f"{name} {sel['name']}: the burnups {outside} are not in the timing nominal-burnup"
            )
    return points


def nominal_points(
    name: str, table: dict, nominal: dict, selections: list[dict]
) -> dict:
    """公称の状態を燃焼度の点 `burnup` だけで解いたとして、その点の公称の状態（少数群定数の項目と、テーブルが持てば微視的断面積と
    数密度の微分）を燃焼度について内挿し（`THINNING_METHODS`）、選ばなかった点で全点のテーブルの公称の状態と比べる（Δρ は
    `summarize_table` と同じ向き、項目の相対の差は `relative_errors`、微視的断面積の差は `micro_errors`）。公称でない状態は使わない。
    テーブルは `thin` と同じ（`burnup_states`）"""
    grid_of = burnup_states(name, table, nominal)
    burnup = grid_of["burnup"]
    positions = grid_of["rows"][grid_of["nominal"]]
    states = [table["states"][p] for p in positions]
    values = np.array([vector(s) for s in states])
    micro = table.get("micro")
    if micro is not None:
        if len(micro["states"]) != len(table["states"]):
            raise SystemExit(
                f"{name}: {len(micro['states'])} micro states for {len(table['states'])} states"
            )
        nuclides = micro["nuclides"]
        layout = MicroLayout(micro_parts(micro["states"][positions[0]], nuclides))
        micro_values = np.array(
            [layout.vector(micro["states"][p], nuclides) for p in positions]
        )
    out = []
    for sel in selections:
        chosen, index, removed = chosen_points(name, sel, burnup)
        x = np.array(chosen)
        q = np.array([burnup[i] for i in removed])
        results = {}
        for method in THINNING_METHODS:
            interpolated = interpolate(method, x, values[index], q)
            rho = []
            relative = dict.fromkeys(RELATIVE, 0.0)
            for n, i in enumerate(removed):
                rebuilt = unvector(interpolated[n], states[i])
                rho.append((1 / k_inf(states[i]) - 1 / k_inf(rebuilt)) * 1e5)
                for item, e in relative_errors(rebuilt, states[i]).items():
                    relative[item] = max(relative[item], e)
            worst = int(np.argmax(np.abs(rho)))
            r = {
                "rho_pcm": rho,
                "max_abs_pcm": abs(rho[worst]),
                "worst_pcm": rho[worst],
                "worst_burnup": burnup[removed[worst]],
                "rms_pcm": float(np.sqrt(np.mean(np.square(rho)))),
                "max_relative": relative,
            }
            if micro is not None:
                found = {}
                item_worst(
                    layout,
                    micro_errors(
                        layout,
                        interpolate(method, x, micro_values[index], q),
                        micro_values[removed],
                    ),
                    [{"burnup": burnup[i]} for i in removed],
                    found,
                )
                r |= ordered_items(found)
            results[method] = r
        out.append(
            {
                "name": sel["name"],
                "burnup": chosen,
                "removed": q.tolist(),
                "results": results,
            }
        )
    return {"name": name, "burnup": burnup, "selections": out}


def nominal_at(
    method: str,
    values: np.ndarray,
    burnup: list[float],
    at_index: list[int],
    removed: list[int],
) -> np.ndarray:
    """公称の状態の値 `values`（[燃焼度の点][値]）の、除いた点 `removed` の値: 公称の点 `at_index` ではテーブルの値、ほかは公称の
    点から内挿した値"""
    held = values[removed].copy()
    between = [n for n, i in enumerate(removed) if i not in at_index]
    if between:
        held[between] = interpolate(
            method,
            np.array([burnup[i] for i in at_index]),
            values[at_index],
            np.array([burnup[removed[n]] for n in between]),
        )
    return held


def combine(
    name: str,
    table: dict,
    nominal: dict,
    selections: list[dict],
    pairs: list[dict],
    timing: dict | None = None,
) -> dict:
    """公称の状態を組の `nominal` の選び方の燃焼点だけで、公称でない状態を `off-nominal` の選び方の燃焼点だけで解いたとして、
    core の組み立て（state.rs の `Assembly::corner`。公称の状態を公称の点から内挿し、公称でない状態の公称の状態との差を、その状態を
    解いた点の公称の状態との差から内挿して足す。ADR-0200 の「core の内挿」）の値を、公称でない点から除いた点で、全点のテーブルの
    すべての状態（公称の状態を含む）と比べる。公称の点では公称の状態はテーブルの値（core は格子点で重みが厳密に 1）。内挿は
    `THINNING_METHODS` で、公称の状態と差に同じ方式を当てる。出す量は `thin` と同じ（全体と状態の種類ごとの最大、微視的断面積の
    差）に、公称の状態の最大（`nominal`）と、少数群定数の項目の相対の差の最大（`relative_errors`）を足したもの。

    core の燃焼度の軸の点は公称の状態を解いた点で、公称でない状態はその軸の点の番号で解いた点を持つので、公称でない点は公称の点の
    中にあること。`timing`（`thinning_time`）があれば、組ごとに、公称の状態と状態によらない秒を公称の点の数で、公称でない状態の秒を
    公称でない点の数で数えた時間を書く（`table_hours`。`nominal-burnup` は組が持つので書かない）"""
    grid_of = burnup_states(name, table, nominal)
    burnup, names, points = grid_of["burnup"], grid_of["names"], grid_of["points"]
    rows = grid_of["rows"]
    at_nominal = grid_of["nominal"]
    states = table["states"]
    nominal_values = np.array([vector(states[p]) for p in rows[at_nominal]])
    micro = table.get("micro")
    if micro is not None:
        if len(micro["states"]) != len(states):
            raise SystemExit(
                f"{name}: {len(micro['states'])} micro states for {len(states)} states"
            )
        nuclides = micro["nuclides"]
        layout = MicroLayout(
            micro_parts(micro["states"][rows[at_nominal][0]], nuclides)
        )

        def micro_rows(index):
            return np.array(
                [layout.vector(micro["states"][p], nuclides) for p in rows[index]]
            )

        nominal_micro = micro_rows(at_nominal)
    kinds_of = {
        index: state_kind(names, index, at_nominal)
        for index in rows
        if index != at_nominal
    }
    order = sorted(set(kinds_of.values()), key=lambda k: kind_order(k, names))
    counts = {k: sum(1 for v in kinds_of.values() if v == k) for k in order}
    if timing is not None and "nominal-burnup" in timing:
        raise SystemExit(
            f"{name}: the pairs give the nominal burnups; drop the timing nominal-burnup"
        )
    seconds = None if timing is None else thinning_time(timing, counts)
    by_name = {}
    for sel in selections:
        if sel["name"] in by_name:
            raise SystemExit(f"{name}: two selections named {sel['name']}")
        chosen_points(name, sel, burnup)
        by_name[sel["name"]] = sel
    out = []
    for pair in pairs:
        unknown = [
            pair[k] for k in ("nominal", "off-nominal") if pair[k] not in by_name
        ]
        if unknown:
            raise SystemExit(
                f"{name}: unknown selections {unknown}; the entry has {list(by_name)}"
            )
        nominal_sel, off_sel = by_name[pair["nominal"]], by_name[pair["off-nominal"]]
        label = f"{pair['nominal']} / {pair['off-nominal']}"
        at_points, at_index, _ = chosen_points(name, nominal_sel, burnup)
        off_points, off_index, removed = chosen_points(name, off_sel, burnup)
        outside = sorted(set(off_points) - set(at_points))
        if outside:
            raise SystemExit(
                f"{name} {label}: the off-nominal burnups {outside} are not nominal burnups "
                "(core takes the burnup axis from the nominal state, ADR-0200)"
            )
        q = np.array([burnup[i] for i in removed])
        results = {}
        for method in THINNING_METHODS:
            rebuilt_nominal = nominal_at(
                method, nominal_values, burnup, at_index, removed
            )
            worst_per_point = [0.0] * len(removed)
            worst = (0.0, None)
            relative = dict.fromkeys(RELATIVE, 0.0)
            per_kind = {}
            micro_kind = {}
            if micro is not None:
                micro_nominal = nominal_at(
                    method, nominal_micro, burnup, at_index, removed
                )
            for grid_index in rows:
                kind = kinds_of.get(grid_index, "nominal")
                kind_worst = per_kind.setdefault(kind, [0.0, None])
                positions = rows[grid_index]
                values = np.array([vector(states[p]) for p in positions])
                rebuilt_all = rebuilt_nominal.copy()
                if grid_index != at_nominal:
                    rebuilt_all += interpolate(
                        method,
                        np.array(off_points),
                        (values - nominal_values)[off_index],
                        q,
                    )
                at = [
                    {"burnup": burnup[i]}
                    | {a: points[k][grid_index[k]] for k, a in enumerate(names)}
                    for i in removed
                ]
                for n, i in enumerate(removed):
                    truth = states[positions[i]]
                    rebuilt = unvector(rebuilt_all[n], truth)
                    rho = (1 / k_inf(truth) - 1 / k_inf(rebuilt)) * 1e5
                    if abs(rho) > abs(worst_per_point[n]):
                        worst_per_point[n] = rho
                    if worst[1] is None or abs(rho) > abs(worst[0]):
                        worst = (rho, at[n])
                    if kind_worst[1] is None or abs(rho) > abs(kind_worst[0]):
                        kind_worst[:] = [rho, at[n]]
                    for item, e in relative_errors(rebuilt, truth).items():
                        relative[item] = max(relative[item], e)
                if micro is not None:
                    truth = micro_rows(grid_index)
                    rebuilt = micro_nominal.copy()
                    if grid_index != at_nominal:
                        rebuilt += interpolate(
                            method,
                            np.array(off_points),
                            (truth - nominal_micro)[off_index],
                            q,
                        )
                    item_worst(
                        layout,
                        micro_errors(layout, rebuilt, truth[removed]),
                        at,
                        micro_kind.setdefault(kind, {}),
                    )
            kinds = {}
            for k in order:
                kinds[k] = {
                    "max_abs_pcm": abs(per_kind[k][0]),
                    "worst_pcm": per_kind[k][0],
                    "worst": per_kind[k][1],
                    "states": counts[k],
                }
                if micro is not None:
                    kinds[k] |= ordered_items(micro_kind[k])
            at_worst = per_kind["nominal"]
            results[method] = {
                "rho_pcm": worst_per_point,
                "max_abs_pcm": abs(worst[0]),
                "worst_pcm": worst[0],
                "worst": worst[1],
                "nominal": {
                    "max_abs_pcm": abs(at_worst[0]),
                    "worst_pcm": at_worst[0],
                    "worst_burnup": at_worst[1]["burnup"],
                },
                "max_relative": relative,
                "kinds": kinds,
            }
            if micro is not None:
                overall = {}
                for k in ["nominal", *order]:
                    for group, items in micro_kind[k].items():
                        for item, e in items.items():
                            best = overall.setdefault(group, {}).setdefault(item, e)
                            if e["max_relative"] > best["max_relative"]:
                                overall[group][item] = e
                results[method] |= ordered_items(overall)
        entry = {
            "name": label,
            "nominal_burnup": at_points,
            "off_nominal_burnup": off_points,
            "removed": q.tolist(),
            "results": results,
        }
        if seconds is not None:
            entry["hours"] = table_hours(
                seconds, len(at_points), {k: len(off_points) for k in order}
            )
        out.append(entry)
    summary = {"name": name, "burnup": burnup, "pairs": out}
    if seconds is not None:
        summary["timing"] = dict(timing)
    return summary


def thinning_time(timing: dict, counts: dict) -> dict:
    """テーブルを作る時間の見込みの、燃焼点ごとの秒。`timing` は、すべての状態を解いた燃焼点 1 つの全体の秒 `total`、自己遮蔽を
    解く状態の 1 状態の秒 `solve`、自己遮蔽を再利用する状態の秒 `reuse`、再利用する状態の変えた軸 `reuse-axis`（ADR-0201 の
    ホウ素。その軸を公称から変えた状態は、ほかの軸の同じ状態の自己遮蔽を再利用する）。公称の状態は `solve`。状態によらない秒
    （`other`）は、`total` から状態の秒を引いた残りか、`total` の代わりに `timing` の `other`。状態の種類ごと（`counts` は種類 →
    状態の数）の秒は、状態の数 × 1 状態の秒。`timing` の `states`（種類 → 状態の数。テーブルのどの種類も書く）があれば、`counts` の
    代わりにその数で数え、状態によらない秒は `other` で与える（間引きの誤差を測ったテーブルと別の、本番の形のテーブルの見込み）。数えた状態の数 `counts` と全体の秒 `total` も返す"""
    unknown_keys = sorted(
        set(timing)
        - {
            "total",
            "other",
            "solve",
            "reuse",
            "reuse-axis",
            "states",
            "nominal-burnup",
        }
    )
    if unknown_keys:
        raise SystemExit(f"unknown timing keys {unknown_keys}")
    if ("total" in timing) == ("other" in timing):
        raise SystemExit("give the timing either the total or the other seconds")
    given = timing.get("states")
    if given is not None:
        # total は誤差を測ったテーブルの全体の秒で、別のテーブルの状態の秒を引いても状態によらない秒にならない
        if "other" not in timing:
            raise SystemExit("timing states need the other seconds, not the total")
        if sorted(given) != sorted(counts):
            raise SystemExit(
                f"the timing states {sorted(given)} must name every kind of the table {list(counts)} (0 for kinds the estimated table does not solve)"
            )
        bad = {
            k: v
            for k, v in given.items()
            if isinstance(v, bool) or not isinstance(v, int) or v < 0
        }
        if bad:
            raise SystemExit(
                f"the timing states {bad} must be counts (integers 0 or more)"
            )
        counts = {k: given[k] for k in counts}
    axis = timing["reuse-axis"]
    per_state = {
        k: float(timing["reuse"] if axis in k.split("+") else timing["solve"])
        for k in counts
    }
    kinds = {k: counts[k] * per_state[k] for k in counts}
    states = float(timing["solve"]) + sum(kinds.values())
    if "other" in timing:
        other = float(timing["other"])
        if other < 0:
            raise SystemExit(f"the other seconds {other} are negative")
        total = other + states
    else:
        total = float(timing["total"])
        other = total - states
        if other < 0:
            raise SystemExit(
                f"the states take {states} s, more than the total {timing['total']} s"
            )
    return {
        "kinds": kinds,
        "counts": counts,
        "other": other,
        "nominal": float(timing["solve"]),
        "total": total,
    }


def table_hours(seconds: dict, burnups: int, points: dict) -> float:
    """公称の状態と状態によらない秒を見込むテーブルの公称の状態の燃焼点の数（`burnups`。`nominal_burnup`）で、状態の種類ごとの秒を
    その種類を解く燃焼点の数（`points`）で数えた時間"""
    total = burnups * (seconds["other"] + seconds["nominal"])
    total += sum(seconds["kinds"][k] * n for k, n in points.items())
    return total / 3600


def thinning_plan(
    name: str,
    plan: dict,
    selections: list[dict],
    order: list[str],
    reuse_axis: str,
    seconds: dict,
    burnups: int,
) -> dict:
    """状態の種類ごとに選び方（`plan["kinds"]`。種類 → 選び方の名前）を決めた組の、誤差と時間（`thin` の説明）"""
    by_name = {s["name"]: s for s in selections}
    chosen = plan["kinds"]
    if sorted(chosen) != sorted(order):
        raise SystemExit(
            f"{name} {plan['name']}: choose a selection for every kind {order}, not {sorted(chosen)}"
        )
    unknown = sorted({v for v in chosen.values() if v not in by_name})
    if unknown:
        raise SystemExit(f"{name} {plan['name']}: unknown selections {unknown}")
    for kind in order:
        axes = kind.split("+")
        if reuse_axis not in axes or len(axes) == 1:
            continue
        # 再利用する状態は、再利用する先の状態を解いた燃焼点でしか解けない（ADR-0201）。公称の状態は、どの選び方の点でも解く（`nominal_burnup`）
        base = "+".join(a for a in axes if a != reuse_axis)
        if not set(by_name[chosen[kind]]["burnup"]) <= set(
            by_name[chosen[base]]["burnup"]
        ):
            raise SystemExit(
                f"{name} {plan['name']}: {kind} reuse the shielding of {base}, so solve them at some of its burnups"
            )
    results = {}
    for method in THINNING_METHODS:
        per_kind = {
            k: by_name[chosen[k]]["results"][method]["kinds"][k]["max_abs_pcm"]
            for k in order
        }
        worst = max(order, key=lambda k: per_kind[k])
        results[method] = {"max_abs_pcm": per_kind[worst], "worst_kind": worst}
    return {
        "name": plan["name"],
        "kinds": {k: chosen[k] for k in order},
        "results": results,
        "hours": table_hours(
            seconds, burnups, {k: len(by_name[chosen[k]]["burnup"]) for k in order}
        ),
    }


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
        if any(
            "kinds" not in r for s in t["selections"] for r in s["results"].values()
        ):
            raise SystemExit(
                f"{t['name']}: the thinning has no errors per kind of state; rerun the comparison"
            )
        hours = (
            f"、公称の状態の {len(t['nominal_burnup'])} 点のすべてで {t['hours']:.2f} 時間の見込み"
            if "hours" in t
            else ""
        )
        lines.append(f"{t['name']}: 燃焼度の点 {len(t['burnup'])}{hours}")
        for sel in t["selections"]:
            hours = f"、{sel['hours']:.2f} 時間の見込み" if "hours" in sel else ""
            lines.append(
                f"  {sel['name']}（{len(sel['burnup'])} 点）: {sel['burnup']}{hours}"
            )
            for method, r in sel["results"].items():
                lines.append(
                    f"    {method}: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_pcm']:+.2f}、{r['worst']}）"
                )
                lines += micro_lines("      ", r)
                for kind, k in r["kinds"].items():
                    lines.append(
                        f"      {kind}: |Δρ| の最大 {k['max_abs_pcm']:.2f} pcm（{k['worst_pcm']:+.2f}、{k['worst']}）"
                    )
        for plan in t.get("plans", []):
            lines.append(f"  {plan['name']}: {plan['hours']:.2f} 時間の見込み")
            for method, r in plan["results"].items():
                lines.append(
                    f"    {method}: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_kind']}）"
                )
    for t in summary.get("nominal", []):
        lines.append(f"{t['name']}: 公称の状態の燃焼度の点 {len(t['burnup'])}")
        for sel in t["selections"]:
            lines.append(
                f"  公称の状態を {sel['name']}（{len(sel['burnup'])} 点）から: {sel['burnup']}"
            )
            for method, r in sel["results"].items():
                rel = "、".join(
                    f"{k} {v * 100:.3f}%" for k, v in r["max_relative"].items()
                )
                lines.append(
                    f"    {method}: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_pcm']:+.2f}、{r['worst_burnup']:g} MWd/kg）、"
                    f"rms {r['rms_pcm']:.2f} pcm、相対の差の最大 {rel}"
                )
                lines += micro_lines("      ", r)
    for t in summary.get("combined", []):
        lines.append(f"{t['name']}: 公称の状態と公称でない状態の燃焼点を合わせた組")
        for pair in t["pairs"]:
            hours = f"、{pair['hours']:.2f} 時間の見込み" if "hours" in pair else ""
            lines.append(
                f"  {pair['name']}（公称 {len(pair['nominal_burnup'])} 点、公称でない "
                f"{len(pair['off_nominal_burnup'])} 点）{hours}"
            )
            for method, r in pair["results"].items():
                n = r["nominal"]
                lines.append(
                    f"    {method}: |Δρ| の最大 {r['max_abs_pcm']:.2f} pcm（{r['worst_pcm']:+.2f}、{r['worst']}）、"
                    f"公称の状態 {n['max_abs_pcm']:.2f} pcm（{n['worst_pcm']:+.2f}、{n['worst_burnup']:g} MWd/kg）"
                )
                lines += micro_lines("      ", r)
                for kind, k in r["kinds"].items():
                    lines.append(
                        f"      {kind}: |Δρ| の最大 {k['max_abs_pcm']:.2f} pcm（{k['worst_pcm']:+.2f}、{k['worst']}）"
                    )
    return lines


def micro_lines(indent: str, result: dict) -> list[str]:
    """微視的断面積と数密度の微分の、項目ごとの差の最大"""
    lines = []
    for kind, label in [("micro", "微視的断面積"), ("derivative", "数密度の微分")]:
        if kind in result:
            items = "、".join(
                f"{item} {e['max_relative']:.2e}（{e['nuclide']}、{e['worst']['burnup']:g}）"
                for item, e in result[kind].items()
            )
            lines.append(f"{indent}{label}: {items}")
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
        for t in config.get("table", []):
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
                    t.get("timing"),
                    t.get("plans", []),
                )
            )
        nominal_entries = [
            nominal_points(
                t["name"],
                read(config_path.parent / t["file"]),
                t["nominal"],
                t["selections"],
            )
            for t in config.get("nominal", [])
        ]
        combined = [
            combine(
                t["name"],
                read(config_path.parent / t["file"]),
                t["nominal"],
                t["selections"],
                t["pairs"],
                t.get("timing"),
            )
            for t in config.get("combined", [])
        ]
        summary = {
            "format": "xlc-compare-branch-cross-terms",
            "thinning": thinning,
            "nominal": nominal_entries,
            "combined": combined,
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(tables, config.get("model", [])),
        }
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
