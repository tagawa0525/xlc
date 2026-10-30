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

参照の燃焼度の点ごとに、core の k の参照の k∞ からの差と、それを追う核種の数密度の違いで説明する分（数密度の効き）: 参照の定数に、
core の数密度と参照の均質化した数密度 N̄ᵢ の差の線形の補正 Σᵢ cᵢ (Nᵢ − N̄ᵢ)（係数は参照の微視的断面積。core と同じ選び方）を足した
2 群の k と、参照の定数の 2 群の k の差。核種の数密度の相対差の絶対値が最大の核種（N̄ᵢ が 0 の核種は除く）と、主な核種の相対差。
続ける点の参照の k∞ は公称の組成を条件 j で解いたもので、公称のテーブルの条件 j の状態の k∞ との差（`check_pcm`）は ADR-0187 の確認の
水準（0.06 pcm 程度）のはず。参照の分岐の k∞ と格子計算の履歴の点の k∞ の差（第 1 項）も並べる。目安として、一生を条件 j で燃焼させた
core と第 1 項のテーブルの差、公称の燃焼の core と公称のテーブルの差（`core-history.json`）を同じ点に並べる。
ρ = 1 − 1/k で、差は pcm。
"""

import argparse
import json
from pathlib import Path

from compare_core_history import (
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
    "core_history": "core-history.json",
}
# 数密度の相対差を要約に残す主な核種（記録の本文に書く）
MAJOR = ["U-235", "Pu-239", "Pu-240", "Pu-241", "Xe-135"]


def reference_states(reference):
    """参照のテーブルの、燃焼度ごとの 1 状態（燃焼度、状態の番号）"""
    names = [a["name"] for a in reference["axes"]]
    burnups = reference["axes"][0]["points"]
    if names != ["burnup", "boron", "moderator", "fuel"] or any(
        len(a["points"]) != 1 for a in reference["axes"][1:]
    ):
        raise SystemExit("the reference must be one state per burnup")
    if len(reference["states"]) != len(burnups):
        raise SystemExit("the reference must solve every burnup")
    return [(u, reference["points"].index([b, 0, 0, 0])) for b, u in enumerate(burnups)]


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


def compare(table, reference, history, core, plain, core_history):
    """比べた要約（モジュールの説明）"""
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
    branch_rows = core_history["branch_rows"]
    rows = []
    for u, p in states:
        state = reference["states"][p]
        micro = reference["micro"]["states"][p]
        n_ref = reference["micro"]["densities"][p]
        coeff = coefficients(micro)

        n = at(core_nodes, u, "core")
        q = at(plain_nodes, u, "plain")
        whole = [x for x in branch_rows if same_burnup(x["burnup"], u)]
        if len(whole) != 1:
            raise SystemExit(
                f"the core history summary has no single row at burnup {u}"
            )
        lattice = at(history["points"], u, "lattice history")
        name, largest, major = densities(n, n_ref, nuclides)
        plain_name, plain_largest, plain_major = densities(q, n_ref, nuclides)
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
                "plain_composition_pcm": composition(state, coeff, n_ref, q),
                "density_percent": largest,
                "density_nuclide": name,
                "plain_density_percent": plain_largest,
                "plain_density_nuclide": plain_name,
                "densities_percent": major,
                "plain_densities_percent": plain_major,
                "branch_difference_pcm": whole[0]["difference_pcm"],
                "nominal_difference_pcm": whole[0]["nominal_difference_pcm"],
                "first_term_pcm": rho(state["k_inf"]) - rho(lattice["k_inf"]),
            }
        )
    # 続ける点の参照と、公称のテーブルの条件 j の状態（同じ組成と状態）
    u0, p0 = states[0]
    node = at(core_nodes, u0, "core")
    table_point, _ = t.point(u0, node["moderator"], "core")
    check = rho(reference["states"][p0]["k_inf"]) - rho(
        table["states"][table_point]["k_inf"]
    )
    return {
        "format": "xlc-compare-history-switch",
        "variable": t.history["variable"],
        "continued_burnup": u0,
        "check_pcm": check,
        "rows": rows,
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
    return (
        f"続ける点 {u0:.4g} MWd/kg の参照の k∞ と公称のテーブルの条件 j の状態の k∞ の差は {summary['check_pcm']:+.4f} pcm"
        f"（core − 参照は {rows[0]['difference_pcm']:+.3f} pcm、補正を切った core も {rows[0]['plain_difference_pcm']:+.3f} pcm）。"
        f"参照の分岐と格子計算の履歴の点の k∞ の差は {span('first_term_pcm', rows)} pcm。"
        f"\n{u0:.4g} MWd/kg より後の参照の点 {len(later)} で、履歴の補正の core − 参照は {span('difference_pcm')} pcm"
        f"（最後の点 {last['burnup']:.4g} MWd/kg で {last['difference_pcm']:+.3f} pcm、w = {last['weight']:.6f}）。"
        f"そのうち数密度の効きは {span('composition_pcm')} pcm（最後の点で {last['composition_pcm']:+.3f} pcm）。"
        f"同じ点の目安: 一生を条件 j で燃焼させた core − 第 1 項 {span('branch_difference_pcm')} pcm、公称の燃焼の core − 公称の"
        f"テーブル {span('nominal_difference_pcm')} pcm。"
        f"\n補正を切った core − 参照は {span('plain_difference_pcm')} pcm（最後の点で {last['plain_difference_pcm']:+.3f} pcm）、"
        f"数密度の効きは {span('plain_composition_pcm')} pcm（最後の点で {last['plain_composition_pcm']:+.3f} pcm）。"
        f"\n参照の N̄ᵢ からの数密度の相対差の絶対値の最大は、履歴の補正の core が {largest('density_percent', 'density_nuclide')}、"
        f"補正を切った core が {largest('plain_density_percent', 'plain_density_nuclide')}。最後の点の主な核種: 履歴の補正の core "
        f"{major(last['densities_percent'])}；補正を切った core {major(last['plain_densities_percent'])}"
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
    summary = compare(
        read(here / FILES["residual_table"]),
        read(here / FILES["reference"]),
        read(here / FILES["history"]),
        read(here / FILES["core"]),
        read(here / FILES["plain"]),
        core_history,
    )
    summary["inputs"] = {file: md5(here / file) for file in FILES.values()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
