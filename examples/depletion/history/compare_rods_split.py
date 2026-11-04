"""制御棒を途中まで入れて抜く燃焼（ADR-0215 の確認の制御棒の履歴、ADR-0228、#0190）の core − 参照の定数の側を、今の core の組成で
U の外れと L の外れに分ける（compare_split_core.py の分け方を、制御棒の履歴の値 out の残りに当てたもの）。

    python3 compare_rods_split.py --json rods-split.json     # このディレクトリの出力から JSON を書く
    python3 compare_rods_split.py --report rods-split.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、core が使った 2 つの変数の残りを持つ公称のテーブル（`table-residual-rods-histories.msgpack`。制御棒の変数の値 out の U(1)、L(1)、
N̄ᵢ^j と、制御棒の値 out の公称の状態の定数と係数）、格子計算の参照の燃焼の履歴（`lattice-switch-rods.msgpack`。参照の k∞）と、その R の点の
組成を制御棒を抜いた公称の状態で解いた分岐（`table-switch-rods-branch.msgpack`）と、その追う核種だけを公称の履歴に入れたハイブリッド
（`table-hybrid-switch-rods-branch.msgpack`）、履歴の補正の core の燃焼の履歴（`core-2A-rods-switch.msgpack`）と、その比較の要約
（`rods-switch.json`。md5 と core − 参照を照らす）。

参照の R の点 u のうち切り替えの点より後ごとに、参照の組成の U_switch = f(switch) − f(H_switch) と
L_switch = f(H_switch) − f(公称、out) − Σᵢ cᵢ ΔN̄ᵢ^switch を作り（cᵢ は公称のテーブルの制御棒の値 out の公称の状態の係数。数密度の微分を
持つ核種は微分）、core が足す量からの外れを

- U の外れ δ_U = w U(1) − U_switch（w = H_r。core の点の履歴の変数）
- L の外れ δ_L = (P / P_j) L(1) − L_switch（P / P_j は core の点の数密度のもの。compare_rods_parts.py の `node`）

とする。効きは参照の定数 f(switch) に足した Δρ（ρ = 1 − 1/k の差の pcm、k は 2 群の無限媒質の k）で、項目は compare_split_core.py と同じ:

- 数密度の効き ρ(f(switch) + Σᵢ cᵢ (Nᵢ − N̄ᵢ^switch)) − ρ(f(switch))
- U の外れの分、L の外れの分、両方の外れの分
- 定数の側 = core − 参照 − 数密度の効き、残り = 定数の側 − 両方の外れの分（再構成、k の基準 ρ(f(switch)) − ρ(参照の k∞)、非加法を含む）

成分の大きさは公称の状態の定数に参照の数密度の線形の補正を足した定数に足した効きで、参照に合う重み 効き(U_switch) / 効き(U(1))、
効き(L_switch) / 効き(L(1))。

U の重みを指数で忘れる重みにした外れ: 制御棒の履歴を H_τ（compare_forgetting.py）にした重み w_τ = 1 − compare_forgetting.weight(u, B₀, τ) で
δ_U(τ) = w_τ U(1) − U_switch、τ ごとに U の外れの分と両方の外れの分（L は今のまま）と、その絶対値の最大。

切り替えの点（20 MWd/kg）の分岐は 2A-rods の組成を制御棒を抜いて解いたもので、テーブルの値 out の U(1)、L(1)、N̄ᵢ^j と同じ組成と状態なので、
その差を `switch_point` に出す（2 つの格子計算の再現の程度）。
"""

import argparse
import json
import math
import re
from pathlib import Path

import compare_forgetting as forgetting
from compare_core_history import same_burnup
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
from compare_rods_parts import (
    RECONSTRUCT_TOLERANCE,
    ROD,
    TAUS,
    history_weights,
    node,
    state_point,
)
from compare_rods_switch import VARIABLES, entry, homogenized, nodes, production, scaled
from compare_split import ratio
from compare_weight import condition_table, effect, largest_difference

FILES = {
    "table": "table-residual-rods-histories.msgpack",
    "reference": "lattice-switch-rods.msgpack",
    "switch": "table-switch-rods-branch.msgpack",
    "switch_hybrid": "table-hybrid-switch-rods-branch.msgpack",
    "core": "core-2A-rods-switch.msgpack",
}
SUMMARY = "rods-switch.json"
# rods-switch.json の core − 参照と、ここで作った値の差の上限 (pcm)（同じ入力と式なので丸めの程度）
SUMMARY_TOLERANCE = 1e-9
# 分岐のテーブルの N̄ᵢ と、参照の燃焼の履歴を均質化した数密度の相対差の上限（同じ組成を lattice と Python で体積の平均にした丸めの程度）
DENSITY_TOLERANCE = 1e-9
# 切り替えの点の分岐の U、L とテーブルの値 out の U(1)、L(1) の差の絶対値の上限 (1/cm)。同じ組成と状態を別の 2 本の格子計算（分岐と
# ハイブリッド、残りのテーブルの分岐とハイブリッド）で解いた中性子束の収束の程度（2A で 5.3e−08）。分岐かハイブリッドが別の履歴か制御棒を
# 入れた状態を解けば、残りの大きさ（1e−04 の程度）で外れる
SWITCH_POINT_TOLERANCE = 1e-6
# 切り替えの点の分岐の N̄ᵢ とテーブルの値 out の N̄ᵢ^j の相対差の上限（同じ履歴の同じ点を lattice が同じく均質化したもの）
SWITCH_DENSITY_TOLERANCE = 1e-12
# L(1) の重みを (P / P_j)^p にした外れを比べる p（1 が今の core。#0190 の候補 (j)）
POWERS = (1.0, 1.5, 2.0)
# 参照の組成の ΔN̄ᵢ / ΔN̄ᵢ^j を出す核種（L を 1 つの指標の比で足す形が、組成の差の向きに沿うかを見る。#0190 の候補 (k)）
DIRECTION_NUCLIDES = ("U-235", "Pu-239", "Pu-240", "Pu-241")


def nominal_condition(table):
    """公称のテーブルの公称の状態の条件（ホウ素、（減速材の温度、密度）、燃料）"""
    by_name = {a["name"]: a for a in table["axes"]}
    moderator = by_name["moderator"]
    m = moderator["points"].index(moderator["nominal"])
    return (
        by_name["boron"]["nominal"],
        (moderator["points"][m], moderator["density"][m]),
        by_name["fuel"]["nominal"],
    )


def switch_rows(table, branch, hybrid):
    """分岐とハイブリッドの R の点ごとの（燃焼度、f(switch)、U_switch、L_switch、N̄ᵢ^switch、線形の補正の後の公称の定数、係数、燃焼度の点の番号）"""
    condition, burnups, states = condition_table(branch, "switch branch")
    h_condition, h_burnups, h_states = condition_table(hybrid, "switch hybrid")
    if h_condition != condition or h_burnups != burnups:
        raise SystemExit(
            "the hybrid solves other conditions or burnups than the branch"
        )
    if condition != nominal_condition(table):
        raise SystemExit(f"the branch solves {condition}, not the nominal state")
    nuclides = table["micro"]["nuclides"]
    for t, name in ((branch, "branch"), (hybrid, "hybrid")):
        if t["micro"]["nuclides"] != nuclides:
            raise SystemExit(f"the {name} tracks other nuclides than the table")
    # H_X の N̄ᵢ は X と同じ（体積の平均。違えば置き換えた核種か履歴が違う）
    if hybrid["micro"]["densities"] != branch["micro"]["densities"]:
        raise SystemExit("the hybrid holds other tracked densities than the branch")
    axis = table["axes"][0]["points"]
    rows = []
    for k, (u, p, q) in enumerate(zip(burnups, states, h_states)):
        found = [i for i, v in enumerate(axis) if same_burnup(v, u)]
        if len(found) != 1:
            raise SystemExit(f"the table does not solve the burnup {u}")
        (b,) = found
        point = state_point(table, b)
        factors = coefficients(table["micro"]["states"][point])
        densities = branch["micro"]["densities"][k]
        deltas = [x - y for x, y in zip(densities, table["micro"]["densities"][b])]
        linear = corrected(table["states"][point], factors, deltas)
        fx, fh = branch["states"][p], hybrid["states"][q]
        rows.append(
            {
                "burnup": u,
                "branch": fx,
                "untracked": difference(fx, fh),
                "linearization": difference(fh, linear),
                "densities": densities,
                "linear": linear,
                "factors": factors,
                "point": b,
            }
        )
    return rows


def check_switch_point(table, row):
    """切り替えの点の U_switch、L_switch、N̄ᵢ^switch と、テーブルの値 out の U(1)、L(1)、N̄ᵢ^j の差の絶対値の最大。上限
    （`SWITCH_POINT_TOLERANCE`、`SWITCH_DENSITY_TOLERANCE`）を超えれば、分岐かハイブリッドが制御棒を抜いた参照の組成を解いていないので止める"""
    rod = entry(table, VARIABLES[1], ROD)
    found = [i for i, v in enumerate(rod["burnup"]) if same_burnup(v, row["burnup"])]
    if len(found) != 1:
        raise SystemExit(f"the table has no rods residual at {row['burnup']}")
    (k,) = found
    out = {
        "burnup": row["burnup"],
        "untracked": largest_difference(row["untracked"], rod["untracked"][k]),
        "linearization": largest_difference(
            row["linearization"], rod["linearization"][k]
        ),
        "densities": max(
            abs(x - y) / abs(y)
            for x, y in zip(row["densities"], rod["densities"][k])
            if y != 0.0
        ),
    }
    if (
        max(out["untracked"], out["linearization"]) > SWITCH_POINT_TOLERANCE
        or out["densities"] > SWITCH_DENSITY_TOLERANCE
    ):
        raise SystemExit(
            f"the branch at the switch point {row['burnup']} is not the rods-out residual of the table: {out}"
        )
    return out


def directions(table, x, listed):
    """参照の組成の数密度の差と履歴の分岐の数密度の差の比: 核種ごとの ΔN̄ᵢ / ΔN̄ᵢ^j（`DIRECTION_NUCLIDES`）と、線形の補正の指標の比
    （熱群と高速群の νΣf、熱群の全断面積、線形の補正の反応度）"""
    b = x["point"]
    nominal = table["micro"]["densities"][b]
    rod = entry(table, VARIABLES[1], ROD)
    (k,) = [i for i, v in enumerate(rod["burnup"]) if same_burnup(v, x["burnup"])]
    branch = [a - c for a, c in zip(rod["densities"][k], nominal)]
    own = [a - c for a, c in zip(x["densities"], nominal)]
    factors = x["factors"]
    state = table["states"][state_point(table, b)]

    def index(key, g, d):
        return sum(f[key][g] * v for f, v in zip(factors, d))

    def reactivity(d):
        return rho(k_inf(corrected(state, factors, d))) - rho(k_inf(state))

    names = table["micro"]["nuclides"]
    missing = [n for n in listed if n not in names]
    if missing:
        raise SystemExit(f"the table does not track {missing}")
    return {
        "nuclides": {
            n: ratio(own[names.index(n)], branch[names.index(n)]) for n in listed
        },
        "thermal_nu_fission": ratio(
            production(factors, own), production(factors, branch)
        ),
        "fast_nu_fission": ratio(
            index("nu_fission", 0, own), index("nu_fission", 0, branch)
        ),
        "thermal_total": ratio(index("total", -1, own), index("total", -1, branch)),
        "reactivity": ratio(reactivity(own), reactivity(branch)),
    }


def compare(
    table,
    reference,
    switch,
    switch_hybrid,
    core,
    taus=TAUS,
    direction_nuclides=DIRECTION_NUCLIDES,
):
    """要約（モジュールの説明）"""
    nuclides = table["micro"]["nuclides"]
    names = [h["variable"] for h in table["histories"]]
    if names != VARIABLES:
        raise SystemExit(f"the table holds the history variables {names}")
    continued = reference.get("continued")
    if continued is None:
        raise SystemExit("the reference is not a continued history")
    start = continued["burnup"]
    ref = homogenized(reference, nuclides, table["micro"]["densities"][0])
    core_nodes = nodes(core, "core", nuclides)
    split = switch_rows(table, switch, switch_hybrid)
    at_switch = [x for x in split if same_burnup(x["burnup"], start)]
    if len(at_switch) != 1:
        raise SystemExit(f"the branch does not solve the switch point {start}")
    rows = []
    bases = []
    for x in split:
        u = x["burnup"]
        # 切り替えの点は at_switch と同じ許容で除く（後の点は制御棒を抜いた core のノードと比べる）
        if u <= start or same_burnup(u, start):
            continue
        (r,) = [r for r in ref if same_burnup(r["burnup"], u)]
        for a, b in zip(x["densities"], r["densities"]):
            if abs(a - b) > DENSITY_TOLERANCE * max(abs(b), 1e-30):
                raise SystemExit(
                    f"the branch densities at {u} are not those of the reference"
                )
        found = [n for n in core_nodes if same_burnup(n["burnup"], u)]
        if len(found) != 1:
            raise SystemExit(f"core: no single point at burnup {u}")
        (n,) = found
        if n["state"] != ROD:
            raise SystemExit(f"core: the node at {u} is not in state {ROD}")
        density_weight, weight = history_weights(table, n["history"])
        if density_weight != 0.0:
            raise SystemExit(
                f"core: the moderator density history at {u} is not nominal"
            )
        f, _, untracked_one, linearization_one, index = node(
            table, u, n["densities"], weight
        )
        reconstruct = rho(n["k"]) - rho(k_inf(f))
        if abs(reconstruct) > RECONSTRUCT_TOLERANCE:
            raise SystemExit(
                f"core: the k at {u} differs by {reconstruct} pcm from its constants"
            )
        factors = x["factors"]
        reference_f = x["branch"]
        core_untracked = scaled(untracked_one, weight)
        core_linearization = scaled(linearization_one, index)
        untracked_gap = difference(core_untracked, x["untracked"])
        linearization_gap = difference(core_linearization, x["linearization"])
        both = added(untracked_gap, linearization_gap)
        dense = corrected(
            reference_f,
            factors,
            [a - c for a, c in zip(n["densities"], x["densities"])],
        )
        rho_reference = rho(k_inf(reference_f))
        difference_pcm = rho(n["k"]) - rho(r["k_inf"])
        composition = rho(k_inf(dense)) - rho_reference
        gap = effect(reference_f, both)
        whole = added(dense, both)
        row = {
            "burnup": u,
            "weight": weight,
            "index_ratio": index,
            "difference_pcm": difference_pcm,
            "composition_pcm": composition,
            "constant_pcm": difference_pcm - composition,
            "untracked_gap_pcm": effect(reference_f, untracked_gap),
            "linearization_gap_pcm": effect(reference_f, linearization_gap),
            "gap_pcm": gap,
            "remainder_pcm": difference_pcm - composition - gap,
            # 残り = 再構成 + 恒等 + 非加法 + k の基準（定義から閉じる。compare_split_core.py と同じ 4 つ）
            "reconstruct_pcm": reconstruct,
            "identity_pcm": rho(k_inf(f)) - rho(k_inf(whole)),
            "basis_pcm": rho_reference - rho(r["k_inf"]),
            "additivity_pcm": rho(k_inf(whole)) - rho_reference - composition - gap,
            "without_untracked_pcm": rho(k_inf(added(dense, linearization_gap)))
            - rho(r["k_inf"]),
            "without_linearization_pcm": rho(k_inf(added(dense, untracked_gap)))
            - rho(r["k_inf"]),
        }
        state = x["linear"]
        for key, one, added_part in (
            ("untracked", untracked_one, core_untracked),
            ("linearization", linearization_one, core_linearization),
        ):
            own = effect(state, x[key])
            whole = effect(state, one)
            row |= {
                f"{key}_pcm": own,
                f"scaled_{key}_pcm": effect(state, added_part),
                f"one_{key}_pcm": whole,
                f"{key}_weight": ratio(own, whole),
            }
        w_l = row["linearization_weight"]
        # 参照に合う L の重みを (P / P_j)^p と見た点ごとの p（重みか P / P_j が 0〜1 の外なら None）
        row["linearization_power"] = (
            math.log(w_l) / math.log(index)
            if w_l is not None and 0.0 < w_l < 1.0 and 0.0 < index < 1.0
            else None
        )
        row["directions"] = directions(table, x, direction_nuclides)
        rows.append(row)
        bases.append(
            (
                u,
                reference_f,
                untracked_one,
                x["untracked"],
                linearization_gap,
                untracked_gap,
                linearization_one,
                index,
                x["linearization"],
            )
        )
    scan = []
    for tau in taus:
        points = []
        for (
            u,
            reference_f,
            untracked_one,
            untracked_switch,
            linearization_gap,
            *_,
        ) in bases:
            w = 1.0 - forgetting.weight(u, start, tau)
            untracked_gap = difference(scaled(untracked_one, w), untracked_switch)
            points.append(
                {
                    "burnup": u,
                    "weight": w,
                    "untracked_gap_pcm": effect(reference_f, untracked_gap),
                    "gap_pcm": effect(
                        reference_f, added(untracked_gap, linearization_gap)
                    ),
                }
            )
        scan.append(
            {
                "tau": tau,
                "points": points,
                "largest_untracked_pcm": max(
                    abs(p["untracked_gap_pcm"]) for p in points
                ),
                "largest_pcm": max(abs(p["gap_pcm"]) for p in points),
            }
        )
    powers = []
    for power in POWERS:
        points = []
        for (
            u,
            reference_f,
            _,
            _,
            _,
            untracked_gap,
            linearization_one,
            index,
            linearization_switch,
        ) in bases:
            gap = difference(
                scaled(linearization_one, index**power), linearization_switch
            )
            points.append(
                {
                    "burnup": u,
                    "linearization_gap_pcm": effect(reference_f, gap),
                    "gap_pcm": effect(reference_f, added(untracked_gap, gap)),
                }
            )
        powers.append(
            {
                "power": power,
                "points": points,
                "largest_linearization_pcm": max(
                    abs(p["linearization_gap_pcm"]) for p in points
                ),
                "largest_pcm": max(abs(p["gap_pcm"]) for p in points),
            }
        )
    best = min(scan, key=lambda s: s["largest_pcm"])["tau"]
    positive = [t for t in taus if t is not None and t > 0]
    negative = [t for t in taus if t is not None and t < 0]
    edges = {min(positive, default=None), max(negative, default=None)} - {None}
    return {
        "format": "xlc-compare-rods-split",
        "switch": start,
        "switch_point": check_switch_point(table, at_switch[0]),
        "rows": rows,
        "best": best,
        "best_at_edge": best in edges,
        "forgetting": scan,
        "powers": powers,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(values, form="+.3f"):
        values = [v for v in values if v is not None]
        return f"{min(values):{form}}〜{max(values):{form}}"

    def points(values, form="+.3f"):
        return "、".join("–" if v is None else f"{v:{form}}" for v in values)

    def name(tau):
        return "∞" if tau is None else f"{tau:g}"

    rows = summary["rows"]

    def column(key):
        return [r[key] for r in rows]

    burnups = "、".join(f"{r['burnup']:g}" for r in rows)
    s = summary["switch_point"]
    scan = summary["forgetting"]
    best = next(x for x in scan if x["tau"] == summary["best"])
    infinite = next(x for x in scan if x["tau"] is None)
    lines = [
        (
            f"切り替えの点 {s['burnup']:g} MWd/kg の分岐とテーブルの値 out の差の絶対値の最大: U {s['untracked']:.2e}、L "
            f"{s['linearization']:.2e} 1/cm、N̄ᵢ の相対 {s['densities']:.2e}"
        ),
        (
            f"{summary['switch']:g} MWd/kg より後の R の点 {len(rows)}（{burnups} MWd/kg）: w {points(column('weight'), '.4f')}、"
            f"P / P_j {points(column('index_ratio'), '.4f')}"
        ),
        (
            f"core − 参照 {points(column('difference_pcm'))} pcm、数密度の効き {points(column('composition_pcm'))} pcm、"
            f"定数の側 {points(column('constant_pcm'))} pcm"
        ),
        (
            f"U の外れの分 {points(column('untracked_gap_pcm'))} pcm（{span(column('untracked_gap_pcm'))}）、"
            f"L の外れの分 {points(column('linearization_gap_pcm'))} pcm（{span(column('linearization_gap_pcm'))}）、"
            f"両方の外れの分 {points(column('gap_pcm'))} pcm。定数の側との残り {span(column('remainder_pcm'))} pcm"
            f"（再構成 {span(column('reconstruct_pcm'), '+.1e')}、恒等 {span(column('identity_pcm'), '+.1e')}、k の基準 "
            f"{span(column('basis_pcm'))}、非加法 "
            f"{span(column('additivity_pcm'))} pcm）"
        ),
        (
            f"成分の大きさ: U_switch {points(column('untracked_pcm'), '+.2f')} pcm、w U(1) {points(column('scaled_untracked_pcm'), '+.2f')} pcm、"
            f"U(1) {points(column('one_untracked_pcm'), '+.2f')} pcm、参照に合う重み w_U {points(column('untracked_weight'), '+.4f')}"
            f"（w {points(column('weight'), '.4f')}）"
        ),
        (
            f"L_switch {points(column('linearization_pcm'), '+.2f')} pcm、(P / P_j) L(1) {points(column('scaled_linearization_pcm'), '+.2f')} pcm、"
            f"L(1) {points(column('one_linearization_pcm'), '+.2f')} pcm、参照に合う重み w_L {points(column('linearization_weight'), '+.4f')}"
            f"（P / P_j {points(column('index_ratio'), '.4f')}）"
        ),
        (
            f"見込み（core を流していない。数密度は今の core のまま）: U の外れを 0 にした core − 参照 "
            f"{points(column('without_untracked_pcm'))} pcm、L の外れを 0 にした core − 参照 "
            f"{points(column('without_linearization_pcm'))} pcm"
        ),
        (
            f"U の重みを指数で忘れる w_τ にした両方の外れの分の絶対値の最大が最小の τ は {name(best['tau'])} MWd/kg で "
            f"{best['largest_pcm']:.3f} pcm（U の外れの分の最大 {best['largest_untracked_pcm']:.3f} pcm。∞ は両方 "
            f"{infinite['largest_pcm']:.3f}、U {infinite['largest_untracked_pcm']:.3f} pcm）"
            + (
                "。この τ は比べた列の 0 の側の端で、最小を挟んでいない。"
                if summary["best_at_edge"]
                else "。"
            )
            + "τ ごとの両方 / U の最大 (pcm): "
            + "、".join(
                f"{name(x['tau'])} {x['largest_pcm']:.3f} / {x['largest_untracked_pcm']:.3f}"
                for x in scan
            )
        ),
    ]
    lines.append(
        "L(1) の重みを (P / P_j)^p にした（U は w のまま）両方の外れの分の絶対値の最大 / L の外れの分の絶対値の最大 (pcm): "
        + "、".join(
            f"p = {x['power']:g} {x['largest_pcm']:.2f} / {x['largest_linearization_pcm']:.2f}"
            for x in summary["powers"]
        )
        + f"。参照に合う L の重みを (P / P_j)^p と見た点ごとの p {points(column('linearization_power'), '.2f')}"
    )
    last = rows[-1]["directions"]
    lines.append(
        f"{rows[-1]['burnup']:g} MWd/kg の参照の組成の ΔN̄ᵢ / ΔN̄ᵢ^j: "
        + "、".join(f"{n} {v:.3f}" for n, v in last["nuclides"].items())
        + f"。指標の比: 熱群の νΣf {last['thermal_nu_fission']:.4f}、高速群の νΣf {last['fast_nu_fission']:.4f}、"
        f"熱群の全断面積 {last['thermal_total']:.4f}、線形の補正の反応度 {last['reactivity']:.4f}"
        f"（参照に合う L の重み {rows[-1]['linearization_weight']:.4f}）"
    )
    return re.sub(r"(?<![A-Za-z])-", "−", "\n".join(lines))


def check_summary(summary, md5s, result):
    """`rods-switch.json` が同じ入力から作られ、core − 参照が同じであること"""
    for key in ("table", "reference", "core"):
        file = FILES[key]
        if summary.get("inputs", {}).get(file) != md5s[file]:
            raise SystemExit(f"the rods switch summary was not made from this {file}")
    for row in result["rows"]:
        found = [
            s for s in summary["switch_rows"] if same_burnup(s["burnup"], row["burnup"])
        ]
        if (
            len(found) != 1
            or abs(found[0]["core_pcm"] - row["difference_pcm"]) > SUMMARY_TOLERANCE
        ):
            raise SystemExit(
                f"the rods switch summary has another core gap at {row['burnup']}"
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
    md5s = {file: md5(here / file) for file in [*FILES.values(), SUMMARY]}
    with open(here / SUMMARY, encoding="utf-8") as f:
        summary = json.load(f)
    result = compare(**{k: read(here / v) for k, v in FILES.items()})
    check_summary(summary, md5s, result)
    result["inputs"] = md5s
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
