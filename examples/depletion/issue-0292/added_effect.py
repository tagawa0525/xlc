"""全核種を輸送した燃焼（all/<燃料>/nominal.msgpack）と 154 核種の試行の燃焼（trial/<燃料>/nominal.msgpack）の k∞ の差を、核種ごとの
吸収の変化に分ける（#0292）。

使い方: python3 added_effect.py --json <出力> <試行の履歴> <全核種の履歴> <試行の組成のテーブル> <全核種の組成のテーブル>
          <全核種の組成のテーブルの入力> <足した核種の一覧>
        python3 added_effect.py --report <出力の JSON>...

テーブルは、2 つの履歴の組成から同じ燃焼点の基準条件の 1 状態を、すべての核種を UFG で解いたブランチ計算（ufg_subset.py の組
`all`。試行は ufg-subset/<燃料>/out/all.msgpack、全核種は all-ufg/<燃料>/out/all.msgpack）。燃焼度の点ごとに（燃焼度 0 は組成が同じなので除く）:

- 直接: 1e5 (1/k_試行 − 1/k_全核種)（pcm。2 つのテーブルの k∞）
- 核種 j の項: −1e5 ΔN̄ⱼ σₐⱼ φ̄ / P。σₐⱼ は全核種のテーブルの冷却計算の 1 群の微視的断面積の中性子を失う分（untransported_absorption.py の
  `absorption`）、φ̄ は 2 群の無限媒質の中性子束の和、P = Σ_G νΣf,G φ_G（中性子束の規格は比に入らない）。ΔN̄ⱼ は、足した核種（試行の
  燃焼が輸送しないので、その輸送に吸収が無い）では全核種の履歴の均質化した数密度そのもの、ほかは 2 つの履歴の差。均質化は燃焼の
  領域の数密度 × 体積の和を格子の体積で割る
- 種類ごとの和: 足した核種、両方が輸送するアクチニドでない核種（足した核種の捕獲の先の Sm-149 などを含む）、アクチニド（吸収だけで、
  核分裂の中性子の増減は数えない）。残り = 直接 − 3 つの和（アクチニドの生成の変化、スペクトルと自己遮蔽の変化、1 次の外の分）

1 群の断面積は全核種の組成のもので、試行の組成の断面積との違い（スペクトルの変化）は残りに入る。
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
from untransported_absorption import (
    absorption,
    history_point,
    nominal_states,
    read,
    region_sum,
    two_group_flux,
)
from verification_tables import markdown_table, signed

ACTINIDES = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"}
FUEL_NAMES = {"u": "U", "mox": "MOX", "gd": "Gd"}
# 一覧に出す核種の数（燃焼度の点ごとに、項の絶対値の大きい順）
TOP = 12


def heavy(name: str) -> bool:
    return name.split("-")[0] in ACTINIDES


def lattice_volume(text: str) -> float:
    """格子計算の入力の格子の体積 (cm²): ピッチの 2 乗 × セルの数（入力の `rows` の語の数）"""
    data = tomllib.loads(text)
    cells = sum(len(row.split()) for row in data["rows"])
    return data["pitch"] ** 2 * cells


def k_by_burnup(table: dict) -> tuple[list[float], dict[float, dict]]:
    burnups, states = nominal_states(table)
    return burnups, {b: table["states"][s] for b, s in zip(burnups, states)}


def contributions(before, after, before_table, after_table, volume, added):
    """燃焼度の点ごとの直接の差、核種ごとの項、種類ごとの和と残り"""
    for h in (before, after):
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit("the histories must be lattice depletion histories")
    gained = set(after["transportable"]) - set(before["transportable"])
    if gained != set(added):
        raise SystemExit(
            f"the added list and the gained transported nuclides differ in {sorted(gained ^ set(added))}"
        )
    burnups, after_states = k_by_burnup(after_table)
    before_burnups, before_states = k_by_burnup(before_table)
    if before_burnups != burnups:
        raise SystemExit("the two tables have other burnup points")
    cooling = after_table["cooling"]
    index_a = {n: i for i, n in enumerate(after["nuclides"])}
    index_b = {n: i for i, n in enumerate(before["nuclides"])}
    added = set(added)
    points = []
    for k, burnup in enumerate(burnups):
        if burnup == 0.0:
            continue
        state = after_states[burnup]
        flux = two_group_flux(state)
        production = sum(v * f for v, f in zip(state["nu_fission"], flux))
        lattice_absorption = sum(a * f for a, f in zip(state["absorption"], flux))
        hp_a, hp_b = history_point(after, burnup), history_point(before, burnup)
        nuclides = {}
        fuel_absorption = 0.0
        for name, sigma in zip(cooling["nuclides"], cooling["points"][k]):
            sigma_a = absorption(sigma)
            n_after = region_sum(after, hp_a, index_a[name]) / volume
            n_before = region_sum(before, hp_b, index_b[name]) / volume
            fuel_absorption += n_after * sigma_a * sum(flux)
            change = n_after if name in added else n_after - n_before
            nuclides[name] = {
                "pcm": -1e5 * change * sigma_a * sum(flux) / production,
                "added": name in added,
                "density_ratio": n_after / n_before if n_before > 0.0 else None,
            }
        if fuel_absorption > lattice_absorption:
            raise SystemExit(
                f"at {burnup} MWd/kg the fuel absorbs more than the lattice: check the lattice volume"
            )
        direct = 1e5 * (1.0 / before_states[burnup]["k_inf"] - 1.0 / state["k_inf"])
        sums = {"added_pcm": 0.0, "carried_pcm": 0.0, "heavy_pcm": 0.0}
        for name, v in nuclides.items():
            if v["added"]:
                sums["added_pcm"] += v["pcm"]
            elif heavy(name):
                sums["heavy_pcm"] += v["pcm"]
            else:
                sums["carried_pcm"] += v["pcm"]
        points.append(
            {
                "burnup": burnup,
                "direct_pcm": direct,
                **sums,
                "remainder_pcm": direct - sum(sums.values()),
                "fuel_absorption_fraction": fuel_absorption / lattice_absorption,
                "nuclides": nuclides,
            }
        )
    return points


def top(point: dict) -> list[dict]:
    """アクチニドでない核種を項の絶対値の大きい順に（アクチニドの項は吸収だけで、核分裂の生成の増減を含まないので並べない）"""
    ranked = sorted(
        ((n, v) for n, v in point["nuclides"].items() if not heavy(n)),
        key=lambda kv: -abs(kv[1]["pcm"]),
    )
    return [{"nuclide": n, **v} for n, v in ranked[:TOP]]


def terms_table(summaries: list[dict]) -> str:
    """直接の差と種類ごとの和（pcm）。アクチニドは吸収の項と残りを足して 1 列にする（吸収の項だけでは核分裂の生成の増減が抜ける）"""
    rows = [
        [
            FUEL_NAMES[d["fuel"]],
            f"{p['burnup']:g}",
            signed(p["direct_pcm"], 1),
            signed(p["added_pcm"], 1),
            signed(p["carried_pcm"], 1),
            signed(p["heavy_pcm"] + p["remainder_pcm"], 1),
        ]
        for d in summaries
        for p in d["points"]
    ]
    return markdown_table(
        [
            "燃料",
            "燃焼度 (MWd/kg)",
            "直接",
            "足した核種の吸収",
            "両方が輸送するアクチニドでない核種",
            "アクチニドと残り",
        ],
        rows,
    )


def top_table(summaries: list[dict]) -> str:
    """最後の燃焼度の点の、項の絶対値の大きいアクチニドでない核種"""
    rows = []
    for d in summaries:
        p = d["points"][-1]
        for r in top(p):
            ratio = r["density_ratio"]
            rows.append(
                [
                    f"{FUEL_NAMES[d['fuel']]} {p['burnup']:g}",
                    r["nuclide"] + ("（足した）" if r["added"] else ""),
                    signed(r["pcm"], 1),
                    "—" if ratio is None else f"{ratio:.3g}",
                ]
            )
    return markdown_table(
        ["燃料と燃焼度", "核種", "項 (pcm)", "数密度の比（全核種 / 試行）"], rows
    )


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--report", nargs="+", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("files", nargs="*", type=Path)
    args = parser.parse_args()
    if args.report:
        summaries = [json.loads(p.read_text()) for p in args.report]
        print(terms_table(summaries))
        print()
        print(top_table(summaries))
        return
    if args.json is None or len(args.files) != 6:
        parser.error("give --json and six files, or --report")
    before, after, before_table, after_table, table_input, added = args.files
    fuel = after.parent.name
    if fuel not in FUEL_NAMES:
        raise SystemExit(f"{after}: the history is not under all/<fuel>")
    points = contributions(
        read(before),
        read(after),
        read(before_table),
        read(after_table),
        lattice_volume(table_input.read_text()),
        added.read_text().split(),
    )
    summary = {
        "format": "xlc-added-effect",
        "fuel": fuel,
        "points": points,
        "inputs": [{"path": str(p), "md5": md5(p)} for p in args.files],
    }
    args.json.write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
