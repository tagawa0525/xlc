"""#0292、#0303: 核データの版ごとに、データのある全核種を輸送するライブラリ（all_nuclides.py の `VERSIONS`）で燃焼させた同じセルの
履歴を、基準の版と比べる（ADR-0264 の版の判断の材料を、65 核種の組から全核種の組に広げる）。

使い方: python3 compare_versions.py compare --json <出力> --base <版> <版>=<履歴>,<テーブル>,<テーブルの入力>...
        python3 compare_versions.py shares --json <出力> <履歴> <テーブル> <テーブルの入力> <核種の一覧>
        python3 compare_versions.py sets --threshold <pcm> <基準の版の目安の JSON>... -- <版の目安の JSON>...
        python3 compare_versions.py --report <出力の JSON>...

テーブルは、履歴の組成から燃焼度の点 `MARKS` の基準条件の 1 状態を、すべての核種を UFG で解いたブランチ計算（ufg_subset.py の組
`all`。VIII.0 は all-ufg/<燃料>/out/all.msgpack、ほかは all-ufg-<版>/<燃料>/out/all.msgpack）。

- `compare`: 版の履歴が、燃焼度の点、燃焼の領域（物質と体積）、燃焼度 0 の組成（0 でない数密度を核種の名前で）で同じことを
  確かめてから、燃焼点ごとの k∞ と、基準の版との反応度の差 (1/k_基準 − 1/k) × 1e5 pcm を書く（../compare_depletion_versions.py と同じ
  量。燃焼度の点は試行の燃焼のもので、Problem 9 の 17 MWd/kg の点が無いので、そのツールは使わない）。`MARKS` の点では、どれかの版で
  0 でない核種の燃料のピン全体（領域の体積の重みの平均）の数密度を版ごとに（チェーンに無い版は null）書き、テーブルの差を核種に
  分ける: 核種 j の項は 1e5 (aⱼ,基準 − aⱼ,版)、aⱼ = N̄ⱼ σₐⱼ φ̄ / P（added_effect.py の項と同じ量を、版ごとにその版のテーブルの断面積、
  中性子束、生成で出したもの。テーブルに無い核種は 0）。残り = 直接（2 つのテーブルの k∞ の反応度の差）− 項の和で、2 つに分ける:
  項の外の吸収 = 1e5 Δ(A / P) − 項の和（A は 2 群の格子の吸収 Σ Σₐ,g φ̄_g。被覆管と減速材の吸収と、燃料の核種の項の近似（ピンの平均の
  数密度 × 冷却計算の 1 群の断面積 × 格子の中性子束）と格子の吸収の食い違いを含む）と、1/k と吸収 / 生成の差 = 直接 − 1e5 Δ(A / P)
  （2 群の無限格子では k = P / A なので、0 から離れればテーブルの吸収が中性子の収支と合わない）。項は生成 P の違い（ν と核分裂の
  断面積）を核種ごとの吸収の比に配るので、核分裂する核種の生成の差は、その核種の項でなく全核種の項に分かれて入る。片方の版の
  テーブルにしか無い核種（その版のチェーンか評価に無い核種）の項の和は、評価の有無の分の 1 次の目安（捕獲の先の核種の数密度の変化は、
  両方にある核種の項に入る）
- `shares`: 版の履歴とテーブルから、一覧の核種の格子の吸収に占める割合を、吸収の割合の JSON の形で書く（dilution_shielding.py が
  版ごとに自己遮蔽の目安の順位を付けるため。added_effect.py の `--shares` と同じ量）
- `sets`: 閾値を超える核種の組（ufg_subset.py の `threshold_sets`。どれかの燃料で超える核種）を、基準の版と比べて標準出力に書く
- `--report`: `compare` の JSON を文章で描く
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from added_effect import contributions, lattice_volume, shares
from compare_depletion_versions import close, initial, labelled
from ufg_subset import threshold_sets
from untransported_absorption import (
    absorption,
    history_point,
    nominal_states,
    read,
    region_sum,
    two_group_flux,
)

# 記録に書く燃焼度の点 (MWd/kg)。テーブルの燃焼度の点（ufg_subset.py の `BURNUP`）と同じ
MARKS = [0.0, 10.0, 30.0, 60.0]
# 報告に出す、項の絶対値の大きい核種の数
TOP = 12


def check_same(histories: dict, base: str) -> list[float]:
    """版の履歴が同じ点、領域、初めの組成を持つことを確かめ、燃焼度の点を返す"""
    if base not in histories:
        raise SystemExit(f"base {base} is not one of {sorted(histories)}")
    for v, h in histories.items():
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit(f"{v}: not an xlc depletion history")
    ref = histories[base]
    burnup = [p["burnup"] for p in ref["points"]]
    regions = [(r["material"], r["volume"]) for r in ref["regions"]]
    composition = initial(ref)
    for v, h in histories.items():
        other = [p["burnup"] for p in h["points"]]
        if len(other) != len(burnup) or not all(map(close, other, burnup)):
            raise SystemExit(f"{v}: burnup points {other}, {base} has {burnup}")
        mine = [(r["material"], r["volume"]) for r in h["regions"]]
        if len(mine) != len(regions) or any(
            m != n or not close(a, b) for (m, a), (n, b) in zip(mine, regions)
        ):
            raise SystemExit(f"{v}: regions {mine}, {base} has {regions}")
        for r, (a, b) in enumerate(zip(initial(h), composition)):
            if a.keys() != b.keys() or not all(close(a[n], b[n]) for n in b):
                raise SystemExit(
                    f"{v}: region {r} starts from another composition than {base}"
                )
    missing = [m for m in MARKS if not any(close(b, m) for b in burnup)]
    if missing:
        raise SystemExit(f"no burnup points at {missing} MWd/kg")
    return burnup


def pin_density(history: dict, index: int, name: str) -> float | None:
    """燃焼点 `index` の、燃料のピン全体（領域の体積の重みの平均）の数密度（チェーンに無い核種は None）"""
    if name not in history["nuclides"]:
        return None
    i = history["nuclides"].index(name)
    volumes = [r["volume"] for r in history["regions"]]
    densities = history["points"][index]["densities"]
    return sum(region[i] * v for region, v in zip(densities, volumes)) / sum(volumes)


def absorptions(history: dict, table: dict, volume: float) -> list[dict]:
    """テーブルの燃焼度の点ごとの k∞ と、核種ごとの aⱼ = N̄ⱼ σₐⱼ φ̄ / P"""
    burnups, states = nominal_states(table)
    if burnups != MARKS:
        raise SystemExit(f"the table has the burnup points {burnups}, not {MARKS}")
    cooling = table["cooling"]
    index = {n: i for i, n in enumerate(history["nuclides"])}
    out = []
    for k, (burnup, s) in enumerate(zip(burnups, states)):
        state = table["states"][s]
        flux = two_group_flux(state)
        production = sum(v * f for v, f in zip(state["nu_fission"], flux))
        lattice_absorption = sum(a * f for a, f in zip(state["absorption"], flux))
        point = history_point(history, burnup)
        terms = {}
        for name, sigma in zip(cooling["nuclides"], cooling["points"][k]):
            n = region_sum(history, point, index[name]) / volume
            terms[name] = n * absorption(sigma) * sum(flux) / production
        # 格子の体積を誤ると項がそろって伸び縮みし、残りに黙って入る（added_effect.py の `contributions` と同じ検査）
        if sum(terms.values()) * production > lattice_absorption:
            raise SystemExit(
                f"at {burnup} MWd/kg the fuel absorbs more than the lattice: check the lattice volume"
            )
        out.append(
            {
                "burnup": burnup,
                "k_inf": state["k_inf"],
                "absorption_per_production": lattice_absorption / production,
                "terms": terms,
            }
        )
    return out


def split(base: list[dict], version: list[dict]) -> list[dict]:
    """テーブルの点ごとの直接の差、核種ごとの項、残り"""
    out = []
    for i, (b, v) in enumerate(zip(base, version)):
        names = list(b["terms"]) + [n for n in v["terms"] if n not in b["terms"]]
        nuclides = {
            n: {"pcm": 1e5 * (b["terms"].get(n, 0.0) - v["terms"].get(n, 0.0))}
            for n in names
        }
        direct = 1e5 * (1.0 / b["k_inf"] - 1.0 / v["k_inf"])
        ratio = 1e5 * (b["absorption_per_production"] - v["absorption_per_production"])
        terms = sum(x["pcm"] for x in nuclides.values())
        out.append(
            {
                "burnup": b["burnup"],
                "index": i,
                "direct_pcm": direct,
                "remainder_pcm": direct - terms,
                "outside_pcm": ratio - terms,
                "balance_pcm": direct - ratio,
                # 片方のテーブルにしか無い核種（その版のチェーンか評価に無い核種）の項の和
                "only_version_pcm": sum(
                    x["pcm"] for n, x in nuclides.items() if n not in b["terms"]
                ),
                "only_base_pcm": sum(
                    x["pcm"] for n, x in nuclides.items() if n not in v["terms"]
                ),
                "nuclides": nuclides,
            }
        )
    return out


def summarize(versions: dict, base: str) -> dict:
    """`versions` は版 → (履歴, テーブル, 格子の体積)。モジュールの説明の `compare`"""
    histories = {v: h for v, (h, _, _) in versions.items()}
    burnup = check_same(histories, base)
    marks = [next(i for i, b in enumerate(burnup) if close(b, m)) for m in MARKS]
    k_ref = [p["k_inf"] for p in histories[base]["points"]]
    tables = {v: absorptions(h, t, vol) for v, (h, t, vol) in versions.items()}
    out = {}
    for v, h in histories.items():
        k = [p["k_inf"] for p in h["points"]]
        entry = {"k_inf": k, "table_k_inf": [p["k_inf"] for p in tables[v]]}
        if v != base:
            entry["reactivity_pcm"] = [(1 / a - 1 / b) * 1e5 for a, b in zip(k_ref, k)]
            entry["split"] = split(tables[base], tables[v])
        out[v] = entry
    names = []
    for h in histories.values():
        names += [n for n in h["nuclides"] if n not in names]
    densities = {}
    for n in names:
        values = {
            v: [pin_density(h, i, n) for i in marks] for v, h in histories.items()
        }
        if any(x for xs in values.values() for x in xs):
            densities[n] = values
    return {
        "base": base,
        "burnup": burnup,
        "marks": [{"burnup": burnup[i], "index": i} for i in marks],
        "versions": out,
        "densities": densities,
    }


def version_shares(history: dict, table: dict, volume: float, names: list[str]) -> dict:
    """版の組成の、核種 `names` の格子の吸収に占める割合（吸収の割合の JSON の形。added_effect.py の `shares`）"""
    return shares(contributions(history, history, table, table, volume, []), names)


def set_difference(base: list[dict], version: list[dict], threshold: float) -> dict:
    """閾値 `threshold` を超える核種の組（燃料の目安の JSON の並びから）を、基準の版と比べる"""
    (a,) = threshold_sets(base, (threshold,)).values()
    (b,) = threshold_sets(version, (threshold,)).values()
    return {
        "base": a,
        "version": b,
        "only_base": [n for n in a if n not in b],
        "only_version": [n for n in b if n not in a],
    }


def report(summary: dict) -> list[str]:
    base = summary["base"]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"基準 {base}。反応度の差は (1/k_{base} − 1/k) × 1e5 pcm",
    ]
    others = [v for v in summary["versions"] if v != base]
    for m in summary["marks"]:
        i = m["index"]
        ks = ", ".join(
            f"{v} {e['k_inf'][i]:.5f}" for v, e in summary["versions"].items()
        )
        diffs = ", ".join(
            f"{v} − {base} {summary['versions'][v]['reactivity_pcm'][i]:+.1f} pcm"
            for v in others
        )
        lines.append(f"{m['burnup']:g} MWd/kg: k∞ {ks}; 反応度の差 {diffs}")
    for v in others:
        for p in summary["versions"][v]["split"]:
            ranked = sorted(p["nuclides"].items(), key=lambda kv: -abs(kv[1]["pcm"]))
            terms = ", ".join(f"{n} {x['pcm']:+.1f}" for n, x in ranked[:TOP])
            lines.append(
                f"{v} − {base} の {p['burnup']:g} MWd/kg のテーブル: 直接 {p['direct_pcm']:+.1f}、"
                f"残り {p['remainder_pcm']:+.1f}（項の外の吸収 {p['outside_pcm']:+.1f}、"
                f"1/k と吸収 / 生成の差 {p['balance_pcm']:+.1f}）、"
                f"{v} だけの核種 {p['only_version_pcm']:+.1f}、{base} だけの核種 {p['only_base_pcm']:+.1f} pcm; 項 {terms}"
            )
    return lines


def check_inputs(history: Path, table: Path, table_input: Path) -> None:
    """テーブルの入力が、与えた履歴の組成を解いて、与えたテーブルを書くものであること（版の取り違えで、ある版の数密度とほかの版の
    断面積を組み合わせないため）"""
    data = tomllib.loads(table_input.read_text())
    place = table_input.parent
    if (place / data["history"]["file"]).resolve() != history.resolve():
        raise SystemExit(
            f"{table_input} solves {data['history']['file']}, not {history}"
        )
    if (place / data["few-group"]["output"]).resolve() != table.resolve():
        raise SystemExit(
            f"{table_input} writes {data['few-group']['output']}, not {table}"
        )
    # 版の置き場は <置き場>/library と <置き場>/<燃料>（all_nuclides.py）。チェーンは履歴を書いた燃焼の入力（履歴の隣の .toml）のもの
    library = (place / data["library"]).resolve()
    if library.parent.parent != history.resolve().parent.parent:
        raise SystemExit(
            f"{table_input} reads the library {data['library']}, not that of {history}"
        )
    burn = history.with_suffix(".toml")
    chain = burn.parent / tomllib.loads(burn.read_text())["depletion"]["chain"]
    if (place / data["history"]["chain"]).resolve() != chain.resolve():
        raise SystemExit(
            f"{table_input} reads the chain {data['history']['chain']}, not that of {burn}"
        )


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def compare_command(args: list[str]) -> None:
    if len(args) < 5 or args[0] != "--json" or args[2] != "--base":
        raise SystemExit(__doc__)
    output, base = Path(args[1]), args[3]
    versions, inputs = {}, []
    for label, value in labelled(args[4:]).items():
        files = [Path(p) for p in str(value).split(",")]
        if len(files) != 3:
            raise SystemExit(f"{label}: expected <history>,<table>,<table input>")
        history, table, table_input = files
        check_inputs(history, table, table_input)
        versions[label] = (
            read(history),
            read(table),
            lattice_volume(table_input.read_text()),
        )
        inputs += [{"version": label, "path": str(p), "md5": md5(p)} for p in files]
    summary = {
        "format": "xlc-compare-all-nuclides-versions",
        "inputs": inputs,
        **summarize(versions, base),
    }
    output.write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


def main() -> None:
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "--report":
        for p in args[1:]:
            print("\n".join(report(json.loads(Path(p).read_text()))))
        return
    if args[:1] == ["compare"]:
        compare_command(args[1:])
        return
    if args[:1] == ["shares"]:
        parser = argparse.ArgumentParser(prog="compare_versions.py shares")
        parser.add_argument("--json", type=Path, required=True)
        parser.add_argument("files", nargs=4, type=Path)
        a = parser.parse_args(args[1:])
        history, table, table_input, names = a.files
        check_inputs(history, table, table_input)
        out = version_shares(
            read(history),
            read(table),
            lattice_volume(table_input.read_text()),
            names.read_text().split(),
        )
        out["inputs"] = [{"path": str(p), "md5": md5(p)} for p in a.files]
        a.json.write_text(json.dumps(out, indent=1) + "\n")
        return
    if args[:1] == ["sets"]:
        if len(args) < 3 or args[1] != "--threshold" or "--" not in args:
            raise SystemExit(__doc__)
        k = args.index("--")
        threshold = float(args[2])
        load = [json.loads(Path(p).read_text()) for p in args[3:k]]
        other = [json.loads(Path(p).read_text()) for p in args[k + 1 :]]
        d = set_difference(load, other, threshold)
        print(
            f"閾値 {threshold:g} pcm: 基準 {len(d['base'])} 核種、版 {len(d['version'])} 核種"
        )
        print("基準だけ: " + (", ".join(d["only_base"]) or "なし"))
        print("版だけ: " + (", ".join(d["only_version"]) or "なし"))
        return
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
