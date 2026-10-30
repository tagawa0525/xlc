"""履歴の補正の core（ADR-0215）の確認: R を足す core の集合体 1 体の燃焼を、テーブルと比べる。

    python3 compare_core_history.py --json core-history.json     # このディレクトリの出力から JSON を書く
    python3 compare_core_history.py --report core-history.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、R を持つ公称のテーブル（`lattice histories residual.toml` の出力）、R を書き足す前の公称のテーブル、第 1 項のテーブル（履歴の
分岐の燃焼を条件 j で解いた分岐）、R の比較の要約（`residual.json`。条件 j の R を公称の状態で使った誤差）と、core の 4 つの燃焼の履歴:

- 公称（`core-2A-residual-nominal.toml`。R を持つテーブル、565 K で一生。w = 0）と、同じ燃焼の R を持たないテーブル
  （`core-2A-plain-nominal.toml`）。2 つのファイルがバイトで同じか（`identical`）
- 条件 j（`core-2A-residual-branch.toml`。600 K で一生。w = 1）: R の燃焼度の点ごとに、core の k の第 1 項のテーブルの k∞ からの差と、
  公称の燃焼の core の k の、公称のテーブルの公称の状態の k∞ からの差（ミクロ燃焼の core と格子計算の組成の違いの目安）
- 途中で条件を変えた燃焼（`core-2A-residual-switch.toml`。20 MWd/kg まで 565 K、その後 600 K）: 燃焼点ごとに、ノードの減速材の密度の
  履歴 H、w = (H − H⁰) / (Hᵇ − H⁰)、誤差の上限 w × max |公称の状態で使った誤差|。R の燃焼度の点では、R の効き
  ρ(f + w R) − ρ(f) と、R の点ごとの誤差の上限 w × |公称の状態で使った誤差(u)|

f はノードの定数のうち R の前まで（テーブルのその状態の格子点の定数 + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰)。Nᵢ は core の数密度。core と同じ係数）で、
R の燃焼度の点の、状態が格子点の燃焼点だけで組み立てる。`reconstruct_pcm` は ρ(f + w R) の 2 群の k と core の k の差（core が R を
足したことの検査。1 ノードの反射境界の k は 2 群の無限媒質の k）。ρ = 1 − 1/k で、差は pcm。
"""

import argparse
import json
from pathlib import Path

from compare_residual import added, coefficients, corrected, k_inf, md5, read, rho

FILES = {
    "residual_table": "table-residual-nominal-histories.msgpack",
    "plain_table": "table-residual-nominal.msgpack",
    "branch_table": "table-residual-density-branch.msgpack",
    "nominal": "core-2A-residual-nominal.msgpack",
    "plain": "core-2A-plain-nominal.msgpack",
    "branch": "core-2A-residual-branch.msgpack",
    "switch": "core-2A-residual-switch.msgpack",
    "residual_summary": "residual.json",
}
VARIABLE = "moderator-density"


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
        """ノードの定数（格子点の定数 + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰)）の 2 群の k と、w R を足した 2 群の k"""
        p, b = self.point(node["burnup"], node["moderator"], label)
        micro = self.table["micro"]
        deltas = [x - y for x, y in zip(node["densities"], micro["densities"][b])]
        f = corrected(self.table["states"][p], coefficients(micro["states"][p]), deltas)
        if w == 0.0:
            return k_inf(f), k_inf(f)
        r = self.rod["residuals"][self.residual_index(node["burnup"])]
        scaled = {
            key: [[w * x for x in row] for row in value]
            if key in ("scatter", "adf")
            else [w * x for x in value]
            for key, value in r.items()
        }
        return k_inf(f), k_inf(added(f, scaled))


def compare(table, branch_table, nominal, branch, switch, residual_summary, identical):
    """比べた要約（モジュールの説明）"""
    t = Table(table)
    nuclides = table["micro"]["nuclides"]
    nominal_nodes = single_node(nominal, "nominal", nuclides)
    branch_nodes = single_node(branch, "branch", nuclides)
    switch_nodes = single_node(switch, "switch", nuclides)
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
    applied = {x["burnup"]: abs(x["applied_pcm"]) for x in residual_summary["rows"]}
    applied_max = max(applied.values())

    def at_residual_points(nodes):
        return [
            (n, k) for n in nodes if (k := t.residual_index(n["burnup"])) is not None
        ]

    branch_points = at_residual_points(branch_nodes)
    nominal_points = at_residual_points(nominal_nodes)
    if [k for _, k in branch_points] != [k for _, k in nominal_points]:
        raise SystemExit(
            "the nominal and the branch core runs reach other burnup points of the residuals"
        )
    branch_rows = []
    for (n, k), (m, _) in zip(branch_points, nominal_points):
        u = t.rod["burnup"][k]
        first = next(
            s
            for v, s in zip(first_burnups, branch_table["states"])
            if same_burnup(u, v)
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
            row["applied_bound_pcm"] = abs(w) * applied[t.rod["burnup"][k]]
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

    def span(rows, key):
        values = [x[key] for x in rows if x[key] is not None]
        return f"{min(values):+.3f}〜{max(values):+.3f}"

    branch = summary["branch_rows"]
    switch = summary["switch_rows"]
    later = [x for x in branch if x["burnup"] >= 15.0]
    weighted = [x for x in switch if x["weight"] != 0.0]
    last = switch[-1]
    return (
        f"公称の燃焼の R を持つテーブルと持たないテーブルの出力は"
        f"{'バイトで同じ' if summary['identical'] else '違う'}。"
        f"\n条件 j の燃焼（w = 1）の core の k と第 1 項のテーブルの k∞ の差は、R の点 {len(branch)} で "
        f"{span(branch, 'difference_pcm')} pcm（15 MWd/kg 以後 {span(later, 'difference_pcm')} pcm）。公称の燃焼の core と公称の"
        f"テーブルの k∞ の差は {span(branch, 'nominal_difference_pcm')} pcm（15 MWd/kg 以後 "
        f"{span(later, 'nominal_difference_pcm')} pcm）。core の k と、格子点の定数 + 線形の補正 + R の 2 群の k の差は "
        f"{span(branch, 'reconstruct_pcm')} pcm。"
        f"\n途中で条件を変えた燃焼: w は {weighted[0]['burnup']:.4g} MWd/kg から 0 でなく、最後の点 {last['burnup']:.4g} MWd/kg で "
        f"H = {last['history']:.6f} g/cc、w = {last['weight']:.6f}、誤差の上限 w × {summary['applied_max_pcm']:.2f} = "
        f"{last['bound_pcm']:.2f} pcm。R の点の R の効きは {span(switch, 'residual_pcm')} pcm、点ごとの誤差の上限は "
        f"{span(switch, 'applied_bound_pcm')} pcm、core の k と組み立てた 2 群の k の差は {span(switch, 'reconstruct_pcm')} pcm"
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
    summary = compare(
        read(here / FILES["residual_table"]),
        read(here / FILES["branch_table"]),
        read(here / FILES["nominal"]),
        read(here / FILES["branch"]),
        read(here / FILES["switch"]),
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
