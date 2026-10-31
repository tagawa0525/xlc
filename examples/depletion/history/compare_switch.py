"""途中で条件を変えた燃焼（ADR-0215 の確認）: 格子計算の参照の k∞ と、履歴の補正の core と補正を切った core の k を比べる。

    python3 compare_switch.py --json switch.json     # このディレクトリの出力から JSON を書く
    python3 compare_switch.py --report switch.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、R を持つ公称のテーブル（`lattice histories residual.toml` の出力。w の H⁰ と Hᵇ、続ける点の条件 j の状態の k∞）、
格子計算の途中で条件を変えた燃焼の履歴（`switch-density-history.toml`。2A の公称の履歴の 20 MWd/kg の点から条件 j で続けた）、
その燃焼点を条件 j で解いた参照（`switch-density-branch.toml`。燃焼度ごとに 1 状態）、core の 2 つの燃焼の履歴（どちらも 20 MWd/kg
まで公称の状態、その後を条件 j）と、一生を条件 j で燃焼させた core の比較の要約（`core-history.json`）:

- 履歴の補正の core（`core-2A-residual-switch.toml`。H_m を燃焼度で重みを付けた平均で進め、R と Δσ を w で足す）
- 補正を切った core（`core-2A-plain-switch.toml`。R と Δσ を持たないテーブルで同じ燃焼）
- 反応率の Δσ を P / P_j で足した core（`core-2A-residual-switch-index.toml`。`rate-residual-weight = "production-index"`。履歴の補正の
  core と同じ燃焼で、H も同じ）

参照の燃焼度の点ごとに、core の k の参照の k∞ からの差と、それを追う核種の数密度の違いで説明する分（数密度の効き）: 参照の定数に、
core の数密度と参照の均質化した数密度 N̄ᵢ の差の線形の補正 Σᵢ cᵢ (Nᵢ − N̄ᵢ)（係数は core と同じ、R を持つ公称のテーブルの
同じ燃焼度の条件 j の状態のもの。Pu-240 は数密度の微分）を足した
2 群の k と、参照の定数の 2 群の k の差。核種の数密度の相対差の絶対値が最大の核種（N̄ᵢ が 0 の核種は除く）と、主な核種の相対差。
続ける点の参照の k∞ は公称の組成を条件 j で解いたもので、公称のテーブルの条件 j の状態の k∞ との差（`check_pcm`）は、同じ組成と
状態の 2 つの分岐の解の差なので、MOC の収束の程度で 0 のはず。参照の分岐の k∞ と格子計算の履歴の点の k∞ の差も並べる。目安として、一生を条件 j で燃焼させた
core と第 1 項のテーブルの差、公称の燃焼の core と公称のテーブルの差（`core-history.json`）を同じ点に並べる。
ρ = 1 − 1/k で、差は pcm。
"""

import argparse
import json
from itertools import pairwise
from pathlib import Path

from compare_core_history import (
    LARGEST,
    Table,
    largest_density_difference,
    same_burnup,
    single_node,
)
from compare_residual import coefficients, corrected, k_inf, md5, read, rho

FILES = {
    "residual_table": "table-residual-nominal-histories.msgpack",
    "reference": "table-switch-density-branch.msgpack",
    "history": "lattice-switch-density.msgpack",
    "core": "core-2A-residual-switch.msgpack",
    "plain": "core-2A-plain-switch.msgpack",
    "index": "core-2A-residual-switch-index.msgpack",
    # 補正を切った core のテーブル（読まず、md5 を残す）
    "plain_table": "table-residual-nominal.msgpack",
    "core_history": "core-history.json",
}
# 数密度の相対差を要約に残す主な核種（記録の本文に書く）
MAJOR = ["U-235", "Pu-239", "Pu-240", "Pu-241", "Xe-135"]
# 重みの見込みを出す ω（PARCS の ω = 2〜2.5 と、その間の値。ADR-0215 の案 V2）
OMEGAS = [1.25, 1.5, 2.0, 2.5]


def omega_weights(t, nodes, omega):
    """core の燃焼点ごとの w を、H(B₀ + ΔB) = (H(B₀) B₀ + ω x ΔB) / (B₀ + ωΔB) で進め直したもの。ステップの x の w は、core が ω = 1 で
    進めた H の差から戻す（w₁(B₁) B₁ − w₁(B₀) B₀) / ΔB）"""
    first = t.weight(nodes[0]["history"]) + 0.0
    weights = [first]
    for a, b in pairwise(nodes):
        b0, b1 = a["burnup"], b["burnup"]
        step = b1 - b0
        x = (t.weight(b["history"]) * b1 - t.weight(a["history"]) * b0) / step
        weights.append((weights[-1] * b0 + omega * x * step) / (b0 + omega * step))
    return weights


def reference_states(reference):
    """参照のテーブルの、燃焼度ごとの 1 状態（燃焼度、燃焼度の点の番号、状態の番号）。N̄ᵢ は燃焼度の点の順、定数は状態の順"""
    names = [a["name"] for a in reference["axes"]]
    burnups = reference["axes"][0]["points"]
    if names != ["burnup", "boron", "moderator", "fuel"] or any(
        len(a["points"]) != 1 for a in reference["axes"][1:]
    ):
        raise SystemExit("the reference must be one state per burnup")
    if len(reference["states"]) != len(burnups):
        raise SystemExit("the reference must solve every burnup")
    return [
        (u, b, reference["points"].index([b, 0, 0, 0])) for b, u in enumerate(burnups)
    ]


def at(nodes, burnup, label):
    """燃焼度 `burnup` の燃焼点"""
    found = [n for n in nodes if same_burnup(n["burnup"], burnup)]
    if len(found) != 1:
        raise SystemExit(f"{label}: no single point at burnup {burnup}")
    return found[0]


def densities(node, reference_densities, nuclides):
    """数密度の相対差: 絶対値が最大の核種と相対差 (%)、主な核種の相対差 (%)"""
    name, largest = largest_density_difference(
        node["densities"], reference_densities, nuclides
    )
    major = {
        n: (node["densities"][i] / reference_densities[i] - 1.0) * 100.0
        for i, n in enumerate(nuclides)
        if n in MAJOR and reference_densities[i] != 0.0
    }
    return name, largest, major


def composition(state, coeff, reference_densities, node):
    """数密度の効き: 参照の定数 `state` に Σᵢ cᵢ (Nᵢ − N̄ᵢ) を足した 2 群の k と、参照の定数の 2 群の k の差 (pcm)"""
    deltas = [x - y for x, y in zip(node["densities"], reference_densities)]
    return rho(k_inf(corrected(state, coeff, deltas))) - rho(k_inf(state))


def composition_nuclides(state, coeff, reference_densities, node, nuclides):
    """核種ごとの数密度の効き（その核種だけの差を足した効き）の、絶対値の大きい `LARGEST` 核種"""
    deltas = [x - y for x, y in zip(node["densities"], reference_densities)]
    k0 = rho(k_inf(state))
    out = []
    for i, name in enumerate(nuclides):
        only = [d if j == i else 0.0 for j, d in enumerate(deltas)]
        out.append(
            {"nuclide": name, "pcm": rho(k_inf(corrected(state, coeff, only))) - k0}
        )
    out.sort(key=lambda x: -abs(x["pcm"]))
    return out[:LARGEST]


def check_core_history(core_history, md5s):
    """一生を条件 j で燃焼させた比較の要約が、同じ R を持つテーブル、core の途中で条件を変えた燃焼、R を持たないテーブル（補正を切った
    core のテーブル。md5 `md5s`）から作ったものであること"""
    for key in ("residual_table", "core", "plain_table"):
        made_from = core_history.get("inputs", {}).get(FILES[key])
        if made_from != md5s[FILES[key]]:
            raise SystemExit(
                f"the core history summary was made from {FILES[key]} of md5 {made_from}, not {md5s[FILES[key]]}"
            )


def compare(table, reference, history, core, plain, core_history, index):
    """比べた要約（モジュールの説明）。`index` は反応率の Δσ を P / P_j で足した core"""
    t = Table(table)
    nuclides = table["micro"]["nuclides"]
    if reference["micro"]["nuclides"] != nuclides:
        raise SystemExit("the reference tracks other nuclides than the table")
    states = reference_states(reference)
    continued = history.get("continued")
    if continued is None or not same_burnup(continued["burnup"], states[0][0]):
        raise SystemExit(
            "the lattice history must continue from the first burnup of the reference"
        )
    core_nodes = single_node(core, "core", nuclides)
    plain_nodes = single_node(plain, "plain", nuclides)
    index_nodes = single_node(index, "index", nuclides)
    branch_rows = core_history["branch_rows"]
    rows = []
    for u, b, p in states:
        state = reference["states"][p]
        n_ref = reference["micro"]["densities"][b]
        n = at(core_nodes, u, "core")
        # 係数は core と同じ、R を持つ公称のテーブルのノードの状態（条件 j）の格子点のもの（微分を持つ核種は微分。ADR-0221）
        table_point, _ = t.point(u, n["moderator"], "core")
        coeff = coefficients(table["micro"]["states"][table_point])
        q = at(plain_nodes, u, "plain")
        x_index = at(index_nodes, u, "index")
        # Δσ の重みだけが違う 2 つの燃焼（同じ状態のステップと同じ H）
        if x_index["history"] != n["history"]:
            raise SystemExit(
                f"the production-index core has the history {x_index['history']} at burnup {u}, not "
                f"{n['history']} of the core"
            )
        whole = [x for x in branch_rows if same_burnup(x["burnup"], u)]
        if len(whole) != 1:
            raise SystemExit(
                f"the core history summary has no single row at burnup {u}"
            )
        lattice = at(history["points"], u, "lattice history")
        name, largest, major = densities(n, n_ref, nuclides)
        plain_name, plain_largest, plain_major = densities(q, n_ref, nuclides)
        index_name, index_largest, index_major = densities(x_index, n_ref, nuclides)
        rows.append(
            {
                "burnup": u,
                # 公称の履歴の w は −0 になりうる（H − H⁰ = 0 を負の数で割る）
                "weight": t.weight(n["history"]) + 0.0,
                "reference_k_inf": state["k_inf"],
                "k": n["k"],
                "plain_k": q["k"],
                "difference_pcm": rho(n["k"]) - rho(state["k_inf"]),
                "plain_difference_pcm": rho(q["k"]) - rho(state["k_inf"]),
                "composition_pcm": composition(state, coeff, n_ref, n),
                "composition_nuclides": composition_nuclides(
                    state, coeff, n_ref, n, nuclides
                ),
                "plain_composition_pcm": composition(state, coeff, n_ref, q),
                "density_percent": largest,
                "density_nuclide": name,
                "plain_density_percent": plain_largest,
                "plain_density_nuclide": plain_name,
                "densities_percent": major,
                "plain_densities_percent": plain_major,
                "branch_difference_pcm": whole[0]["difference_pcm"],
                "nominal_difference_pcm": whole[0]["nominal_difference_pcm"],
                "reference_history_pcm": rho(state["k_inf"]) - rho(lattice["k_inf"]),
                "index_k": x_index["k"],
                "index_difference_pcm": rho(x_index["k"]) - rho(state["k_inf"]),
                "index_composition_pcm": composition(state, coeff, n_ref, x_index),
                "index_composition_nuclides": composition_nuclides(
                    state, coeff, n_ref, x_index, nuclides
                ),
                "index_density_percent": index_largest,
                "index_density_nuclide": index_name,
                "index_densities_percent": index_major,
            }
        )
        row = rows[-1]
        fill = row["plain_difference_pcm"] - row["difference_pcm"]
        # 補正の量が w に比例すると置いて、参照に合う w
        row["needed_weight"] = (
            None
            if row["weight"] == 0.0 or fill == 0.0
            else row["weight"] * row["plain_difference_pcm"] / fill
        )
    # 続ける点の参照と、公称のテーブルの条件 j の状態（同じ組成と状態）
    u0, _, p0 = states[0]
    node = at(core_nodes, u0, "core")
    table_point, _ = t.point(u0, node["moderator"], "core")
    check = rho(reference["states"][p0]["k_inf"]) - rho(
        table["states"][table_point]["k_inf"]
    )
    omegas = []
    for omega in OMEGAS:
        weights = omega_weights(t, core_nodes, omega)
        at_rows = [
            next(
                w
                for n, w in zip(core_nodes, weights)
                if same_burnup(n["burnup"], x["burnup"])
            )
            for x in rows
        ]
        omegas.append(
            {
                "omega": omega,
                "weights": at_rows,
                # 補正の量が w に比例すると置いた差の見込み（core を再実行していない）
                "estimate_pcm": [
                    None
                    if x["weight"] == 0.0
                    else x["plain_difference_pcm"]
                    - w
                    / x["weight"]
                    * (x["plain_difference_pcm"] - x["difference_pcm"])
                    for x, w in zip(rows, at_rows)
                ],
            }
        )
    return {
        "format": "xlc-compare-history-switch",
        "variable": t.history["variable"],
        "continued_burnup": u0,
        "check_pcm": check,
        "rows": rows,
        "omegas": omegas,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""
    rows = summary["rows"]
    later = rows[1:]
    last = rows[-1]

    def span(key, items=later):
        values = [x[key] for x in items]
        return f"{min(values):+.3f}〜{max(values):+.3f}"

    def largest(key, name):
        worst = max(later, key=lambda x: abs(x[key]))
        return f"{worst[key]:+.3f}%（{worst[name]}、{worst['burnup']:.4g} MWd/kg）"

    def major(values):
        return "、".join(f"{n} {v:+.3f}%" for n, v in values.items())

    u0 = summary["continued_burnup"]

    def omega_text(x):
        values = [v for v in x["estimate_pcm"] if v is not None]
        if not values:
            return f"ω = {x['omega']:g} は見込みなし（補正の効く点が無い）"
        return f"ω = {x['omega']:g} で最後の点の w {x['weights'][-1]:.4f}、見込み {min(values):+.1f}〜{max(values):+.1f} pcm"

    needed = [x for x in later if x["needed_weight"] is not None]
    return (
        f"続ける点 {u0:.4g} MWd/kg の参照の k∞ と公称のテーブルの条件 j の状態の k∞ の差は {summary['check_pcm']:+.4f} pcm"
        f"（core − 参照は {rows[0]['difference_pcm']:+.3f} pcm、補正を切った core も {rows[0]['plain_difference_pcm']:+.3f} pcm）。"
        f"参照の分岐と格子計算の履歴の点の k∞ の差は {span('reference_history_pcm', rows)} pcm。"
        f"\n{u0:.4g} MWd/kg より後の参照の点 {len(later)} で、履歴の補正の core − 参照は {span('difference_pcm')} pcm"
        f"（最後の点 {last['burnup']:.4g} MWd/kg で {last['difference_pcm']:+.3f} pcm、w = {last['weight']:.6f}）。"
        f"そのうち数密度の効きは {span('composition_pcm')} pcm（最後の点で {last['composition_pcm']:+.3f} pcm）、"
        f"数密度の効きを除いた残り（core − 参照 − 数密度の効き）は "
        f"{min(x['difference_pcm'] - x['composition_pcm'] for x in later):+.3f}〜"
        f"{max(x['difference_pcm'] - x['composition_pcm'] for x in later):+.3f} pcm。"
        f"最後の点の数密度の効きの大きい核種: "
        + "、".join(
            f"{x['nuclide']} {x['pcm']:+.2f}" for x in last["composition_nuclides"]
        )
        + " pcm。"
        f"同じ点の目安: 一生を条件 j で燃焼させた core − 第 1 項 {span('branch_difference_pcm')} pcm、公称の燃焼の core − 公称の"
        f"テーブル {span('nominal_difference_pcm')} pcm。"
        f"\n反応率の Δσ を P / P_j で足した core − 参照は {span('index_difference_pcm')} pcm（最後の点で "
        f"{last['index_difference_pcm']:+.3f} pcm）、数密度の効きは {span('index_composition_pcm')} pcm（最後の点で "
        f"{last['index_composition_pcm']:+.3f} pcm）、定数の側は "
        f"{min(x['index_difference_pcm'] - x['index_composition_pcm'] for x in later):+.3f}〜"
        f"{max(x['index_difference_pcm'] - x['index_composition_pcm'] for x in later):+.3f} pcm。"
        f"参照の N̄ᵢ からの数密度の相対差の絶対値の最大は {largest('index_density_percent', 'index_density_nuclide')}。"
        f"最後の点の数密度の効きの大きい核種: "
        + "、".join(
            f"{x['nuclide']} {x['pcm']:+.2f}"
            for x in last["index_composition_nuclides"]
        )
        + " pcm。主な核種の相対差: "
        f"{major(last['index_densities_percent'])}"
        f"\n補正を切った core − 参照は {span('plain_difference_pcm')} pcm（最後の点で {last['plain_difference_pcm']:+.3f} pcm）、"
        f"数密度の効きは {span('plain_composition_pcm')} pcm（最後の点で {last['plain_composition_pcm']:+.3f} pcm）。"
        f"\n参照の N̄ᵢ からの数密度の相対差の絶対値の最大は、履歴の補正の core が {largest('density_percent', 'density_nuclide')}、"
        f"補正を切った core が {largest('plain_density_percent', 'plain_density_nuclide')}。最後の点の主な核種: 履歴の補正の core "
        f"{major(last['densities_percent'])}；補正を切った core {major(last['plain_densities_percent'])}"
        f"\n補正の量が w に比例すると置いた見込み: 参照に合う w は "
        + "、".join(
            f"{x['burnup']:.4g} MWd/kg で {x['needed_weight']:.4f}（ω = 1 で {x['weight']:.4f}）"
            for x in needed
        )
        + "。"
        + "；".join(omega_text(x) for x in summary["omegas"])
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", help="write the summary to this file")
    parser.add_argument("--report", help="print the report of this summary")
    args = parser.parse_args()
    if args.report:
        with open(args.report, encoding="utf-8") as f:
            print(report(json.load(f)))
        return
    here = Path(__file__).resolve().parent
    with open(here / FILES["core_history"], encoding="utf-8") as f:
        core_history = json.load(f)
    md5s = {file: md5(here / file) for file in FILES.values()}
    check_core_history(core_history, md5s)
    summary = compare(
        read(here / FILES["residual_table"]),
        read(here / FILES["reference"]),
        read(here / FILES["history"]),
        read(here / FILES["core"]),
        read(here / FILES["plain"]),
        core_history,
        read(here / FILES["index"]),
    )
    summary["inputs"] = md5s
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
