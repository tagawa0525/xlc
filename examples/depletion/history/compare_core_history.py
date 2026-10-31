"""履歴の補正の core（ADR-0215）の確認: R を足す core の集合体 1 体の燃焼を、テーブルと比べる。

    python3 compare_core_history.py --json core-history.json     # このディレクトリの出力から JSON を書く
    python3 compare_core_history.py --report core-history.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、R を持つ公称のテーブル（`lattice histories residual.toml` の出力）、R を書き足す前の公称のテーブル、第 1 項のテーブル（履歴の
分岐の燃焼を条件 j で解いた分岐）、R の比較の要約（`residual.json`。条件 j の R を公称の状態で使った誤差）と、core の 4 つの燃焼の履歴:

- 公称（`core-2A-residual-nominal.toml`。R を持つテーブル、565 K で一生。w = 0）と、同じ燃焼の R を持たないテーブル
  （`core-2A-plain-nominal.toml`）。2 つのファイルがバイトで同じか（`identical`）
- 条件 j（`core-2A-residual-branch.toml`。600 K で一生。w = 1）: R の燃焼度の点ごとに、core の k の第 1 項のテーブルの k∞ からの差と、
  公称の燃焼の core の k の、公称のテーブルの公称の状態の k∞ からの差（ミクロ燃焼の core と格子計算の組成の違いの目安）。core の
  追う核種の数密度の、第 1 項のテーブルの N̄ᵢ^j からの相対差の絶対値が最大の核種と、その相対差（N̄ᵢ^j が 0 の核種は除く）
- 第 1 項のテーブルで燃焼させた core（`core-2A-sigma-branch.toml`。条件 j の燃焼と同じで、テーブルだけが第 1 項）: 同じ点の k の第 1 項の
  k∞ からの差と、数密度の同じ相対差（履歴の補正が目指す水準）
- 途中で条件を変えた燃焼（`core-2A-residual-switch.toml`。20 MWd/kg まで 565 K、その後 600 K）: 燃焼点ごとに、ノードの減速材の密度の
  履歴 H、w = (H − H⁰) / (Hᵇ − H⁰)、R の状態への依存の目安 w × max |公称の状態で使った誤差|（R の全体を w で足すと置いた値。core は L を
  w でなく P / P_j で足すので、L の分の目安にならない）。R の燃焼度の点では、履歴の補正の効き ρ(f + w U + (P / P_j) L) − ρ(f) と、R の点ごとの
  目安 w × |公称の状態で使った誤差(u)|

f はノードの定数のうち履歴の補正の前まで（テーブルのその状態の格子点の定数 + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰)。Nᵢ は core の数密度。core と同じ係数）で、
R の燃焼度の点の、状態が格子点の燃焼点だけで組み立てる。U は w、L は熱群の νΣf の線形の補正の比 P / P_j で足す（`Table.node_k`）。
`reconstruct_pcm` は f + w U + (P / P_j) L の 2 群の k と core の k の差（core が補正を足したことの検査。1 ノードの反射境界の k は 2 群の
無限媒質の k）。ρ = 1 − 1/k で、差は pcm。
"""

import argparse
import json
from pathlib import Path

from compare_residual import (
    added,
    coefficients,
    corrected,
    k_inf,
    md5,
    read,
    residual,
    rho,
)

FILES = {
    "residual_table": "table-residual-nominal-histories.msgpack",
    "plain_table": "table-residual-nominal.msgpack",
    "branch_table": "table-residual-density-branch.msgpack",
    "nominal": "core-2A-residual-nominal.msgpack",
    "plain": "core-2A-plain-nominal.msgpack",
    "branch": "core-2A-residual-branch.msgpack",
    "switch": "core-2A-residual-switch.msgpack",
    "sigma": "core-2A-sigma-branch.msgpack",
    "residual_summary": "residual.json",
}
VARIABLE = "moderator-density"
# 核種ごとの効きを要約に残す、効きの大きい核種の数
LARGEST = 6


# core の燃焼度とテーブルの燃焼度の点の一致の許容（core は出力密度 × EFPD で燃焼度を進める。compare_core_micro.py と同じ）
BURNUP_TOLERANCE = 1e-9


def same_burnup(a, b):
    return abs(a - b) <= BURNUP_TOLERANCE * max(1.0, abs(b))


def single_node(core, label, nuclides):
    """集合体 1 体・ノード 1 つの core の履歴の、燃焼点ごとの（EFPD、燃焼度、状態の名前、減速材の温度、k、数密度、H）"""
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"{label}: not an xlc-core-history")
    if len(core["assemblies"]) != 1 or len(core["axial_mesh"]) != 2:
        raise SystemExit(f"{label}: the core must be one assembly of one axial node")
    if core.get("nuclides") != nuclides:
        raise SystemExit(f"{label}: the core tracks other nuclides than the table")
    if core.get("history_variables") != [VARIABLE]:
        raise SystemExit(f"{label}: the core holds no {VARIABLE} history")
    return [
        {
            "efpd": p["efpd"],
            "burnup": p["burnups"][0][0],
            "state": p["state"],
            "moderator": core["states"][p["state"]]["moderator"],
            "k": p["k"],
            "densities": p["densities"][0][0],
            "history": p["histories"][0][0][0],
        }
        for p in core["points"]
    ]


class Table:
    """R を持つ公称のテーブルの、格子点の定数の引き方"""

    def __init__(self, table):
        names = [a["name"] for a in table["axes"]]
        if names != ["burnup", "boron", "moderator", "fuel"]:
            raise SystemExit(f"the nominal table has the axes {names}")
        histories = table.get("histories") or []
        if len(histories) != 1 or histories[0]["variable"] != VARIABLE:
            raise SystemExit(f"the nominal table must hold one {VARIABLE} history")
        (self.history,) = histories
        rods = {x["rod"]: x for x in self.history["rods"]}
        if list(rods) != ["out"]:
            raise SystemExit("the history must hold the rod value out only")
        self.rod = rods["out"]
        self.table = table
        self.burnups = table["axes"][0]["points"]
        self.moderator = table["axes"][2]["points"]
        self.nominal = [a["points"].index(a["nominal"]) for a in table["axes"][1:]]
        self.points = {tuple(p): i for i, p in enumerate(table["points"])}

    def weight(self, h):
        """w = (H − H⁰) / (Hᵇ − H⁰)"""
        return (h - self.history["nominal"]) / (
            self.history["branch"] - self.history["nominal"]
        )

    def residual_index(self, burnup):
        """R の燃焼度の点の番号（点でなければ None）"""
        return next(
            (k for k, u in enumerate(self.rod["burnup"]) if same_burnup(burnup, u)),
            None,
        )

    def point(self, burnup, moderator, label):
        """燃焼度の点と減速材の格子点の状態（ほかの軸は公称）の、テーブルの点の番号と燃焼度の点の番号"""
        if moderator not in self.moderator:
            raise SystemExit(
                f"{label}: the moderator temperature {moderator} is not a point of the table"
            )
        b = next(i for i, u in enumerate(self.burnups) if same_burnup(burnup, u))
        index = (b, self.nominal[0], self.moderator.index(moderator), self.nominal[2])
        if index not in self.points:
            raise SystemExit(f"{label}: the table has not solved {index}")
        return self.points[index], b

    def node_k(self, node, w, label):
        """ノードの定数（格子点の定数 + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰)）の 2 群の k と、core と同じく w U + (P / P_j) L を足した 2 群の k（ADR-0215）。
        P は熱群（最後の群）の νΣf の線形の補正のノードの数密度の差のもの、P_j は同じ係数の履歴の分岐の数密度の差のもの。履歴の分岐の
        数密度が公称と同じ点は L を足さない（L は定義から 0）"""
        p, b = self.point(node["burnup"], node["moderator"], label)
        micro = self.table["micro"]
        factors = coefficients(micro["states"][p])
        deltas = [x - y for x, y in zip(node["densities"], micro["densities"][b])]
        f = corrected(self.table["states"][p], factors, deltas)
        if w == 0.0:
            return k_inf(f), k_inf(f)
        k = self.residual_index(node["burnup"])
        g = added(f, scaled(self.rod["untracked"][k], w))
        branch = [
            x - y for x, y in zip(self.rod["densities"][k], micro["densities"][b])
        ]
        if any(d != 0.0 for d in branch):
            own = production(factors, deltas)
            of_branch = production(factors, branch)
            if of_branch == 0.0:
                raise SystemExit(f"{label}: the branch index is 0 at {node['burnup']}")
            g = added(g, scaled(self.rod["linearization"][k], own / of_branch))
        return k_inf(f), k_inf(g)


def production(factors, deltas):
    """熱群（最後の群）の νΣf の線形の補正 Σᵢ cᵢ_νf,G Δᵢ（core の `production_index`）"""
    return sum(c["nu_fission"][-1] * d for c, d in zip(factors, deltas))


def scaled(r, w):
    """w 倍した残り"""
    return {
        key: [[w * x for x in row] for row in value]
        if key in ("scatter", "adf")
        else [w * x for x in value]
        for key, value in r.items()
    }


def check_summary(residual_summary, table_md5):
    """R の比較の要約が、同じ R を持つテーブル（md5 `table_md5`）から作ったものであること"""
    made_from = residual_summary.get("inputs", {}).get(FILES["residual_table"])
    if made_from != table_md5:
        raise SystemExit(
            f"the residual summary was made from the table of md5 {made_from}, not {table_md5}"
        )


def capture(sigma):
    """群ごとの捕獲の微視的断面積（MT 102 の経路の和）"""
    groups = len(sigma["fission"])
    return [
        sum(ch[g] for m in sigma["reactions"] if m["mt"] == 102 for ch in m["channels"])
        for g in range(groups)
    ]


def percent(a, b):
    """群ごとの相対差 a / b − 1 (%)。b が 0 の群は None"""
    return [(x / y - 1.0) * 100.0 if y != 0.0 else None for x, y in zip(a, b)]


def composition(t, node, first_densities, first_states):
    """条件 j の燃焼の R の点で、core の数密度と第 1 項のテーブルの数密度 N̄ᵢ^j の違いの効き: 基準は第 1 項を組み立てた定数
    （格子点の定数 + Σᵢ cᵢ (N̄ᵢ^j − N̄ᵢ⁰) + R）で、ρ(基準 + Σᵢ cᵢ (Nᵢ − N̄ᵢ^j)) − ρ(基準) と、核種ごと（その核種だけの差）の効きが大きい `LARGEST` 核種の効き、
    数密度の相対差と、第 1 項のテーブルの微視的断面積の公称のテーブルの条件 j の値に対する相対差（核分裂と捕獲、群ごと）"""
    p, b = t.point(node["burnup"], node["moderator"], "branch")
    micro = t.table["micro"]
    coeff = coefficients(micro["states"][p])
    k = t.residual_index(node["burnup"])
    base = added(
        corrected(
            t.table["states"][p],
            coeff,
            [x - y for x, y in zip(first_densities, micro["densities"][b])],
        ),
        residual(t.rod, k),
    )
    k_base = k_inf(base)
    deltas = [x - y for x, y in zip(node["densities"], first_densities)]
    nuclides = []
    for i, name in enumerate(micro["nuclides"]):
        only = [d if j == i else 0.0 for j, d in enumerate(deltas)]
        own, first = micro["states"][p][i], first_states[i]
        nuclides.append(
            {
                "nuclide": name,
                "pcm": rho(k_inf(corrected(base, coeff, only))) - rho(k_base),
                "density_percent": (node["densities"][i] / first_densities[i] - 1.0)
                * 100.0
                if first_densities[i] != 0.0
                else None,
                "fission_percent": percent(first["fission"], own["fission"]),
                "capture_percent": percent(capture(first), capture(own)),
            }
        )
    nuclides.sort(key=lambda x: -abs(x["pcm"]))
    # 要約を数十 KB に収める（記録が書くのは効きの大きい核種だけ）
    return rho(k_inf(corrected(base, coeff, deltas))) - rho(k_base), nuclides[:LARGEST]


def largest_density_difference(densities, first_densities, nuclides):
    """数密度の、第 1 項のテーブルの数密度に対する相対差 (%) の絶対値が最大の核種と、その相対差（第 1 項の数密度が 0 の核種は除く）"""
    values = [
        (name, (x / y - 1.0) * 100.0)
        for name, x, y in zip(nuclides, densities, first_densities)
        if y != 0.0
    ]
    return max(values, key=lambda v: abs(v[1]))


def compare(
    table, branch_table, nominal, branch, switch, sigma, residual_summary, identical
):
    """比べた要約（モジュールの説明）"""
    t = Table(table)
    nuclides = table["micro"]["nuclides"]
    nominal_nodes = single_node(nominal, "nominal", nuclides)
    branch_nodes = single_node(branch, "branch", nuclides)
    switch_nodes = single_node(switch, "switch", nuclides)
    sigma_nodes = single_node(sigma, "sigma", nuclides)
    first_burnups = branch_table["axes"][0]["points"]
    if [a["name"] for a in branch_table["axes"]] != [
        "burnup",
        "boron",
        "moderator",
        "fuel",
    ] or len(branch_table["states"]) != len(first_burnups):
        raise SystemExit("the first-term table must be one state per burnup")
    if not all(any(same_burnup(u, v) for v in first_burnups) for u in t.rod["burnup"]):
        raise SystemExit(
            "the first-term table does not solve the burnups of the residuals"
        )
    # 条件 j の R を公称の状態で使った誤差の大きさ（R の点ごと）
    applied = []
    for u in t.rod["burnup"]:
        found = [
            abs(x["applied_pcm"])
            for x in residual_summary["rows"]
            if same_burnup(x["burnup"], u)
        ]
        if len(found) != 1:
            raise SystemExit(f"the residual summary has no single row at burnup {u}")
        applied.extend(found)
    applied_max = max(applied)

    def at_residual_points(nodes):
        return [
            (n, k) for n in nodes if (k := t.residual_index(n["burnup"])) is not None
        ]

    branch_points = at_residual_points(branch_nodes)
    nominal_points = at_residual_points(nominal_nodes)
    sigma_points = at_residual_points(sigma_nodes)
    if not (
        [k for _, k in branch_points]
        == [k for _, k in nominal_points]
        == [k for _, k in sigma_points]
    ):
        raise SystemExit(
            "the nominal, the branch and the first-term core runs reach other burnup points of the residuals"
        )
    branch_rows = []
    for (n, k), (m, _), (q, _) in zip(branch_points, nominal_points, sigma_points):
        u = t.rod["burnup"][k]
        f = next(i for i, v in enumerate(first_burnups) if same_burnup(u, v))
        first = branch_table["states"][f]
        first_densities = branch_table["micro"]["densities"][f]
        name, density_difference = largest_density_difference(
            n["densities"], first_densities, nuclides
        )
        sigma_name, sigma_percent = largest_density_difference(
            q["densities"], first_densities, nuclides
        )
        composition_pcm, contributions = composition(
            t,
            n,
            branch_table["micro"]["densities"][f],
            branch_table["micro"]["states"][f],
        )
        p0, _ = t.point(u, m["moderator"], "nominal")
        _, with_r = t.node_k(n, t.weight(n["history"]), "branch")
        branch_rows.append(
            {
                "burnup": u,
                "k": n["k"],
                "table_k_inf": first["k_inf"],
                "difference_pcm": rho(n["k"]) - rho(first["k_inf"]),
                "nominal_difference_pcm": rho(m["k"])
                - rho(table["states"][p0]["k_inf"]),
                "reconstruct_pcm": rho(with_r) - rho(n["k"]),
                "composition_pcm": composition_pcm,
                "nuclides": contributions,
                "density_percent": density_difference,
                "density_nuclide": name,
                "sigma_k": q["k"],
                "sigma_difference_pcm": rho(q["k"]) - rho(first["k_inf"]),
                "sigma_density_percent": sigma_percent,
                "sigma_density_nuclide": sigma_name,
            }
        )
    switch_rows = []
    for n in switch_nodes:
        w = t.weight(n["history"])
        k = t.residual_index(n["burnup"])
        row = {
            "efpd": n["efpd"],
            "burnup": n["burnup"],
            "state": n["state"],
            "k": n["k"],
            "history": n["history"],
            # 公称の履歴の w は −0 になりうる（H − H⁰ = 0 を負の数で割る）
            "weight": w + 0.0,
            "bound_pcm": abs(w) * applied_max,
            "residual_pcm": None,
            "applied_bound_pcm": None,
            "reconstruct_pcm": None,
        }
        if k is not None:
            without_r, with_r = t.node_k(n, w, "switch")
            row["residual_pcm"] = rho(with_r) - rho(without_r)
            row["applied_bound_pcm"] = abs(w) * applied[k]
            row["reconstruct_pcm"] = rho(with_r) - rho(n["k"])
        switch_rows.append(row)
    return {
        "format": "xlc-compare-core-history",
        "variable": VARIABLE,
        "nominal_value": t.history["nominal"],
        "branch_value": t.history["branch"],
        "identical": identical,
        "applied_max_pcm": applied_max,
        "branch_rows": branch_rows,
        "switch_rows": switch_rows,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(rows, key, form="+.3f"):
        values = [x[key] for x in rows if x[key] is not None]
        if not values:
            return "（点なし）"
        return f"{min(values):{form}}〜{max(values):{form}}"

    branch = summary["branch_rows"]
    switch = summary["switch_rows"]

    def largest(rows, key, name):
        worst = max(rows, key=lambda x: abs(x[key]))
        return f"{worst[key]:+.3f}%（{worst[name]}、{worst['burnup']:.4g} MWd/kg）"

    later = [x for x in branch if x["burnup"] >= 15.0]
    weighted = [x for x in switch if x["weight"] != 0.0]
    end = branch[-1]

    def groups(values):
        return "、".join("–" if v is None else f"{v:+.3f}" for v in values)

    composition_text = (
        f"最後の R の点 {end['burnup']:.4g} MWd/kg で、core と第 1 項のテーブルの数密度の違いの効きは "
        f"{end['composition_pcm']:+.3f} pcm（core − 第 1 項 {end['difference_pcm']:+.3f} pcm）。効きの大きい核種: "
        + "；".join(
            f"{x['nuclide']} {x['pcm']:+.2f} pcm（数密度 {groups([x['density_percent']])}%、第 1 項の σ の公称の条件 j に対する差 "
            f"核分裂 {groups(x['fission_percent'])}%、捕獲 {groups(x['capture_percent'])}%）"
            for x in end["nuclides"]
        )
    )
    last = switch[-1]
    first_weight = (
        f"{weighted[0]['burnup']:.4g} MWd/kg から 0 でなく"
        if weighted
        else "どの点も 0"
    )
    return (
        f"公称の燃焼の R を持つテーブルと持たないテーブルの出力は"
        f"{'バイトで同じ' if summary['identical'] else '違う'}。"
        f"\n条件 j の燃焼（w = 1）の core の k と第 1 項のテーブルの k∞ の差は、R の点 {len(branch)} で "
        f"{span(branch, 'difference_pcm')} pcm（15 MWd/kg 以後 {span(later, 'difference_pcm')} pcm）。公称の燃焼の core と公称の"
        f"テーブルの k∞ の差は {span(branch, 'nominal_difference_pcm')} pcm（15 MWd/kg 以後 "
        f"{span(later, 'nominal_difference_pcm')} pcm）。core の k と、格子点の定数 + 線形の補正 + w U + (P / P_j) L の 2 群の k の差は "
        f"{span(branch, 'reconstruct_pcm', '+.1e')} pcm。"
        f"\n第 1 項のテーブルで燃焼させた core の k と第 1 項のテーブルの k∞ の差は {span(branch, 'sigma_difference_pcm')} pcm"
        f"（15 MWd/kg 以後 {span(later, 'sigma_difference_pcm')} pcm）。"
        f"core の追う核種の数密度の第 1 項の N̄ᵢ^j からの相対差の絶対値の最大は、R の点のすべてで "
        f"{largest(branch, 'density_percent', 'density_nuclide')}、第 1 項のテーブルで燃焼させた core は "
        f"{largest(branch, 'sigma_density_percent', 'sigma_density_nuclide')}。"
        f"\n{composition_text}"
        f"\n途中で条件を変えた燃焼: w は "
        f"{first_weight}、"
        f"最後の点 {last['burnup']:.4g} MWd/kg で "
        f"H = {last['history']:.6f} g/cc、w = {last['weight']:.6f}、R の状態への依存の目安 w × {summary['applied_max_pcm']:.2f} = "
        f"{last['bound_pcm']:.2f} pcm。R の点の履歴の補正の効きは {span(switch, 'residual_pcm')} pcm、点ごとの目安は "
        f"{span(switch, 'applied_bound_pcm')} pcm、core の k と組み立てた 2 群の k の差は {span(switch, 'reconstruct_pcm', '+.1e')} pcm"
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
    with open(here / FILES["residual_summary"], encoding="utf-8") as f:
        residual_summary = json.load(f)
    check_summary(residual_summary, md5(here / FILES["residual_table"]))
    summary = compare(
        read(here / FILES["residual_table"]),
        read(here / FILES["branch_table"]),
        read(here / FILES["nominal"]),
        read(here / FILES["branch"]),
        read(here / FILES["switch"]),
        read(here / FILES["sigma"]),
        residual_summary,
        (here / FILES["nominal"]).read_bytes() == (here / FILES["plain"]).read_bytes(),
    )
    summary["inputs"] = {file: md5(here / file) for file in FILES.values()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
