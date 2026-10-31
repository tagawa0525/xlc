"""途中で条件を変えた燃焼（ADR-0215 の確認、#0190）の core − 参照の定数の側を、今の core の組成で U の外れと L の外れに分ける。

    python3 compare_split_core.py --json split-core.json     # このディレクトリの出力から JSON を書く
    python3 compare_split_core.py --report split-core.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、core が使った R を持つ公称のテーブル（`table-residual-nominal-histories.msgpack`。U(1)、L(1)、N̄ᵢ^j と、条件 j の状態の定数と
係数）、途中で条件を変えた燃焼の参照（`table-switch-density-branch.msgpack`）と、その追う核種だけを公称の履歴に入れたハイブリッド
（`table-hybrid-switch-branch.msgpack`）、履歴の補正の core の燃焼の履歴（`core-2A-residual-switch.msgpack`）と、その比較の要約（`switch.json`）。

参照の燃焼度の点 u ごとに、compare_split.py の `split` で参照の組成の U_switch = f(switch) − f(H_switch) と
L_switch = f(H_switch) − f(公称、j) − Σᵢ cᵢ ΔN̄ᵢ^switch を作り（係数 cᵢ は公称のテーブルの条件 j の状態のもの）、core が足す量からの外れを

- U の外れ δ_U = w U(1) − U_switch（w = (H − H⁰) / (Hᵇ − H⁰)、H は core の点の履歴の変数）
- L の外れ δ_L = (P / P_j) L(1) − L_switch（P = Σᵢ cᵢ_νf,G (Nᵢ − N̄ᵢ⁰)、P_j = Σᵢ cᵢ_νf,G ΔN̄ᵢ^j。Nᵢ は core の点の数密度、cᵢ は core と同じ
  公称のテーブルの条件 j の状態の格子点の係数。core の `production_ratio`）

とする。w が 0 の点は core が U も L も足さない（外れは −U_switch と −L_switch）。

core のノードの定数 f_core = f(公称、j) + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰) + w U(1) + (P / P_j) L(1) は、参照の定数
f(switch) = f(公称、j) + Σᵢ cᵢ ΔN̄ᵢ^switch + U_switch + L_switch を使うと、線形の補正が数密度の差に 1 次なので
f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch) + δ_U + δ_L。効きは参照の定数に足した Δρ（ρ = 1 − 1/k の差の pcm、k は 2 群の無限媒質の k）:

- 数密度の効き ρ(f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch)) − ρ(f(switch))（compare_switch.py の `composition`）
- U の外れの分 ρ(f(switch) + δ_U) − ρ(f(switch))、L の外れの分、両方の外れの分
- 定数の側 = core − 参照 − 数密度の効き、残り = 定数の側 − 両方の外れの分

残りは次の 4 つの和（定義から閉じる）:

- 再構成 ρ(core の k) − ρ(f_core)（core がこの通りに定数を組み立てて解いたこと。1 ノードの反射境界の k は 2 群の無限媒質の k）
- 恒等 ρ(f_core) − ρ(f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch) + δ_U + δ_L)（上の式の検算で丸めの程度）
- 非加法 ρ(f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch) + δ_U + δ_L) − ρ(f(switch)) − 数密度の効き − 両方の外れの分（反応度で足し合わない分）
- k の基準 ρ(f(switch)) − ρ(参照の k∞)（2 群の k と格子計算の k∞ の違い）

外れの分は、core のノードの定数から外れを引いた ρ(f_core) − ρ(f_core − δ) も並べる（足す定数による違い）。見込みとして、U の外れを 0
にした core − 参照 ρ(f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch) + δ_L) − ρ(参照の k∞) と、L の外れを 0 にした同じ量（core を流していない。
core の数密度は今のまま）。

成分の大きさは compare_split.py の比べ方（公称の状態の定数に参照の数密度の線形の補正を足した定数に足した効き）で、U_switch、w U(1)、
L_switch、(P / P_j) L(1) と、成分ごとに参照に合う重み 効き(U_switch) / 効き(U(1))、効き(L_switch) / 効き(L(1))。
"""

import argparse
import json
from pathlib import Path

from compare_core_history import Table, production, same_burnup, single_node
from compare_residual import (
    added,
    coefficients,
    corrected,
    difference,
    k_inf,
    md5,
    read,
    rho,
)
from compare_split import base, ratio, solved, split
from compare_switch import at
from compare_switch import composition as composition_effect
from compare_weight import effect, scaled

FILES = {
    "nominal": "table-residual-nominal-histories.msgpack",
    "switch": "table-switch-density-branch.msgpack",
    "switch_hybrid": "table-hybrid-switch-branch.msgpack",
    "core": "core-2A-residual-switch.msgpack",
}
SUMMARY = "switch.json"
# switch.json の値と、ここで同じ式で作り直した値の差の上限（同じ入力と式なので丸めの程度）
SUMMARY_TOLERANCE = 1e-9


def gap_effect(state, gap):
    """状態の定数 `state` に外れ `gap` を足した効き (pcm)"""
    return rho(k_inf(added(state, gap))) - rho(k_inf(state))


def node_gap_effect(node, gap):
    """core のノードの定数 `node` から外れ `gap` を引いた定数に対する効き ρ(node) − ρ(node − gap) (pcm)"""
    return rho(k_inf(node)) - rho(k_inf(added(node, scaled(gap, -1.0))))


def core_parts(t, point, b, node, weight):
    """core が足す w U(1) と (P / P_j) L(1) と P / P_j。w が 0 なら足さない（None）"""
    if weight == 0.0:
        return None, None, None
    k = t.residual_index(node["burnup"])
    if k is None:
        raise SystemExit(f"the burnup {node['burnup']} is not a residual point")
    micro = t.table["micro"]
    factors = coefficients(micro["states"][point])
    nominal = micro["densities"][b]
    branch = [x - y for x, y in zip(t.rod["densities"][k], nominal)]
    if all(d == 0.0 for d in branch):
        raise SystemExit(
            f"the branch densities equal the nominal ones at {node['burnup']}"
        )
    of_branch = production(factors, branch)
    if of_branch == 0.0:
        raise SystemExit(f"the branch index is 0 at {node['burnup']}")
    index = (
        production(factors, [x - y for x, y in zip(node["densities"], nominal)])
        / of_branch
    )
    return (
        scaled(t.rod["untracked"][k], weight),
        scaled(t.rod["linearization"][k], index),
        index,
    )


def check_values(row, s, u):
    """switch.json の点 `s` が、同じ入力でここで作った値と同じこと"""
    pairs = [
        ("weight", row["weight"]),
        ("difference_pcm", row["difference_pcm"]),
        ("composition_pcm", row["composition_pcm"]),
    ]
    if row["index_ratio"] is not None:
        pairs.append(("index_ratio", row["index_ratio"]))
    for key, value in pairs:
        if s[key] is None or abs(s[key] - value) > SUMMARY_TOLERANCE:
            raise SystemExit(f"the switch summary has another {key} at {u}")


def compare(tables, summary):
    """要約（モジュールの説明）"""
    nominal = tables["nominal"]
    t = Table(nominal)
    nuclides = nominal["micro"]["nuclides"]
    condition, rows = split(
        nominal, tables["switch"], tables["switch_hybrid"], "switch"
    )
    at_condition = solved(nominal, condition)
    nodes = single_node(tables["core"], "core", nuclides)
    out = []
    for x in rows:
        u = x["burnup"]
        found = [s for s in summary["rows"] if same_burnup(s["burnup"], u)]
        if len(found) != 1:
            raise SystemExit(f"the switch summary has no single point at {u}")
        n = at(nodes, u, "core")
        point, b = t.point(u, n["moderator"], "core")
        if at_condition[u] != (b, point):
            raise SystemExit(f"the core node at {u} is not at the condition j")
        coeff = coefficients(nominal["micro"]["states"][point])
        # 公称の履歴の w は −0 になりうる（H − H⁰ = 0 を負の数で割る）
        weight = t.weight(n["history"]) + 0.0
        core_untracked, core_linearization, index = core_parts(t, point, b, n, weight)
        reference = x["branch"]
        untracked_gap = (
            scaled(x["untracked"], -1.0)
            if core_untracked is None
            else difference(core_untracked, x["untracked"])
        )
        linearization_gap = (
            scaled(x["linearization"], -1.0)
            if core_linearization is None
            else difference(core_linearization, x["linearization"])
        )
        dense = corrected(
            reference,
            coeff,
            [a - c for a, c in zip(n["densities"], x["densities"])],
        )
        node = corrected(
            nominal["states"][point],
            coeff,
            [a - c for a, c in zip(n["densities"], nominal["micro"]["densities"][b])],
        )
        if core_untracked is not None:
            node = added(added(node, core_untracked), core_linearization)
        both = added(untracked_gap, linearization_gap)
        whole = added(dense, both)
        rho_reference = rho(k_inf(reference))
        row = {
            "burnup": u,
            "weight": weight,
            "index_ratio": index,
            "difference_pcm": rho(n["k"]) - rho(reference["k_inf"]),
            "composition_pcm": composition_effect(reference, coeff, x["densities"], n),
            "untracked_gap_pcm": gap_effect(reference, untracked_gap),
            "linearization_gap_pcm": gap_effect(reference, linearization_gap),
            "gap_pcm": gap_effect(reference, both),
            "node_untracked_gap_pcm": node_gap_effect(node, untracked_gap),
            "node_linearization_gap_pcm": node_gap_effect(node, linearization_gap),
            "node_gap_pcm": node_gap_effect(node, both),
            "reconstruct_pcm": rho(n["k"]) - rho(k_inf(node)),
            "identity_pcm": rho(k_inf(node)) - rho(k_inf(whole)),
            "basis_pcm": rho_reference - rho(reference["k_inf"]),
            "without_untracked_pcm": rho(k_inf(added(dense, linearization_gap)))
            - rho(reference["k_inf"]),
            "without_linearization_pcm": rho(k_inf(added(dense, untracked_gap)))
            - rho(reference["k_inf"]),
        }
        check_values(row, found[0], u)
        row["constant_pcm"] = row["difference_pcm"] - row["composition_pcm"]
        row["remainder_pcm"] = row["constant_pcm"] - row["gap_pcm"]
        row["additivity_pcm"] = (
            rho(k_inf(whole)) - rho_reference - row["composition_pcm"] - row["gap_pcm"]
        )
        state = base(nominal, u, x["densities"])
        for key, added_part in (
            ("untracked", core_untracked),
            ("linearization", core_linearization),
        ):
            k = t.residual_index(u)
            one = None if k is None else effect(state, t.rod[key][k])
            own = effect(state, x[key])
            row |= {
                f"{key}_pcm": own,
                f"scaled_{key}_pcm": 0.0
                if added_part is None
                else effect(state, added_part),
                f"one_{key}_pcm": one,
                # w = 0 の点（続ける点）は core が足さず、参照に合う重みを出さない
                f"{key}_weight": None
                if added_part is None or one is None
                else ratio(own, one),
            }
        out.append(row)
    return {
        "format": "xlc-compare-history-split-core",
        "continued_burnup": summary["continued_burnup"],
        "rows": out,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(values, form="+.3f"):
        values = [v for v in values if v is not None]
        return f"{min(values):{form}}〜{max(values):{form}}"

    def points(values, form="+.3f"):
        return "、".join(f"{v:{form}}" for v in values)

    u0 = summary["continued_burnup"]
    later = [x for x in summary["rows"] if x["burnup"] > u0]

    def column(key):
        return [x[key] for x in later]

    burnups = "、".join(f"{x['burnup']:.4g}" for x in later)
    both = [x["untracked_gap_pcm"] + x["linearization_gap_pcm"] for x in later]
    unadded = [c - b for c, b in zip(column("constant_pcm"), both)]
    node = max(
        abs(x[f"node_{key}_pcm"] - x[f"{key}_pcm"])
        for x in later
        for key in ("untracked_gap", "linearization_gap", "gap")
    )
    lines = [
        (
            f"{u0:.4g} MWd/kg より後の点 {len(later)}（{burnups} MWd/kg）: w {points(column('weight'), '.4f')}、"
            f"P / P_j {points(column('index_ratio'), '.4f')}"
        ),
        (
            f"U の外れの分 {points(column('untracked_gap_pcm'))} pcm（{span(column('untracked_gap_pcm'))}）、"
            f"L の外れの分 {points(column('linearization_gap_pcm'))} pcm（{span(column('linearization_gap_pcm'))}）、"
            f"両方の外れの分 {points(column('gap_pcm'))} pcm（参照の定数に足した効き）"
        ),
        (
            f"core のノードの定数から引いた効きでは U の外れ {span(column('node_untracked_gap_pcm'))}、"
            f"L の外れ {span(column('node_linearization_gap_pcm'))}、両方 {span(column('node_gap_pcm'))} pcm"
            f"（参照の定数に足した効きとの差の絶対値の最大 {node:.4f} pcm）"
        ),
        (
            f"定数の側（core − 参照 − 数密度の効き）{points(column('constant_pcm'))} pcm に対し、U の外れの分 + L の外れの分は "
            f"{points(both)} pcm（定数の側との差 {span(unadded)} pcm）、両方の外れの分との残り "
            f"{span(column('remainder_pcm'))} pcm（再構成 {span(column('reconstruct_pcm'), '+.1e')}、恒等 "
            f"{span(column('identity_pcm'), '+.1e')}、非加法 {span(column('additivity_pcm'))}、k の基準 "
            f"{span(column('basis_pcm'))} pcm）"
        ),
        (
            f"成分の大きさ（公称の状態の定数に参照の数密度の線形の補正を足した定数）: 参照に合う重み w_U "
            f"{points(column('untracked_weight'), '.4f')}（w {points(column('weight'), '.4f')}）、w_L "
            f"{points(column('linearization_weight'), '.4f')}（P / P_j {points(column('index_ratio'), '.4f')}）"
        ),
        (
            f"見込み（core を流していない。数密度は今の core のまま）: U の外れを 0 にした core − 参照 "
            f"{points(column('without_untracked_pcm'))} pcm、L の外れを 0 にした core − 参照 "
            f"{points(column('without_linearization_pcm'))} pcm（今の core − 参照 {points(column('difference_pcm'))} pcm）"
        ),
    ]
    # 記録の表（verification_tables.py の signed）と同じ負の符号
    return "\n".join(lines).replace("-", "−")


def check_summary(summary, md5s):
    """`switch.json` が、ここで読む公称のテーブル、参照、ハイブリッド、core の履歴から作られたこと"""
    for file in FILES.values():
        if summary.get("inputs", {}).get(file) != md5s[file]:
            raise SystemExit(f"the switch summary was not made from this {file}")


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
    md5s = {file: md5(here / file) for file in [*FILES.values(), SUMMARY]}
    with open(here / SUMMARY, encoding="utf-8") as f:
        summary = json.load(f)
    check_summary(summary, md5s)
    tables = {k: read(here / v) for k, v in FILES.items()}
    result = compare(tables, summary)
    result["inputs"] = md5s
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
