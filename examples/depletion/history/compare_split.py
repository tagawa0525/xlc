"""履歴の残り R（ADR-0215）を、追わない核種の分 U と線形化の誤差 L に分け、w への依存を成分ごとに比べる（#0190）。

    python3 compare_split.py --json split.json     # このディレクトリの出力から JSON を書く
    python3 compare_split.py --report split.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

履歴 X（条件を変えた燃焼の組成）と、その条件の 1 状態のテーブルごとに、公称の履歴の組成のうち R の第 3 項と同じ追う核種（分岐の入力の
`micro`）だけを X の値にしたハイブリッド H_X（make_hybrid.py の --micro）を同じ条件で解いたテーブルを読み、燃焼度の点 u ごとに

- U_X = f(X) − f(H_X)（追う核種の数密度が同じで、追わない核種だけが違う 2 つの組成の定数の差）
- L_X = f(H_X) − f(公称) − Σᵢ cᵢ ΔN̄ᵢ（追う核種だけが違う組成の定数の差から、core と同じ係数の線形の補正を引いたもの。係数は公称の組成の
  その条件の状態の値で、数密度の微分を持つ核種は微分。ΔN̄ᵢ = N̄ᵢ^X − N̄ᵢ^0 で、H_X の N̄ᵢ は X と同じ）
- R_X = U_X + L_X（compare_residual.py の R と同じ式）

を作る。N̄ᵢ は燃焼の領域の体積の平均なので、H_X の N̄ᵢ は X と同じ（違えば止める）。

効きは次の 2 つ（ρ = 1 − 1/k の差の pcm。k は少数群定数の 2 群の無限媒質の k。compare_residual.py の `k_inf`）:

- 公称の状態の定数での効き: 公称のテーブル（R を書き足した `table-residual-nominal-histories.msgpack`）の公称の状態（565 K、0.743 g/cc）に、その履歴の数密度の
  線形の補正 Σᵢ cᵢ⁰ ΔN̄ᵢ を足した定数に、成分を足した効き（compare_weight.py の比べ方。同じ定数に足すので、成分どうしを比べられる）
- 条件での効き: ρ(f(X)) − ρ(f(H_X))、ρ(f(H_X)) − ρ(f(公称) + Σᵢ cᵢ ΔN̄ᵢ)、ρ(f(X)) − ρ(f(公称) + Σᵢ cᵢ ΔN̄ᵢ)。和が閉じる

比べるもの:

1. R(1)（一生を条件 j（600 K、0.661 g/cc）で燃焼させた履歴の分岐。`residual-density-history.toml`）の U と L の大きさ
2. 途中で条件を変えた燃焼（`switch-density-history.toml`）の参照の点ごとの U と L を、core の w（ω = 1。`switch.json` の `weight`）を R(1) の
   成分に掛けた w U(1)、w L(1) と、その点の履歴の定数で比べ、成分ごとに参照に合う重み w_U = 効き(U_switch) / 効き(U(1))、w_L を出す
3. core の差（`switch.json` の `difference_pcm`）のうち、成分の比例からの外れで説明できる分: core の定数は参照の定数から
   R_switch − w R(1) を引き、数密度の違いの線形の補正を足したものなので、外れ δ_U = U_switch − w U(1) の分を ρ(参照の定数 − δ_U) − ρ(参照の
   定数) とする（δ_L、両方も同じ）。残り = core の差 − 数密度の効き（`switch.json` の `composition_pcm`）− 両方の外れの分
4. 減速材の密度だけの経路（565 K のまま、0.702 と 0.661 g/cc）の R(0.5) と R(1)/2 の違いを成分ごとに、w = 0.5 の定数（公称の状態に w = 0.5
   の数密度の線形の補正を足したもの）に足した効きで比べる（compare_weight.py の比べ方）。比 = 効き(X(0.5)) / (2 効き(X(1)/2))（比例なら 1/2）
"""

import argparse
import json
from pathlib import Path

from compare_core_history import Table, same_burnup
from compare_residual import (
    added,
    coefficients,
    corrected,
    difference,
    k_inf,
    md5,
    read,
    residual,
    rho,
)
from compare_weight import condition_table, effect, largest_difference, scaled

FILES = {
    # R を書き足した公称のテーブル（定数とミクロ燃焼のデータは R を持たないテーブルと同じ。R(1) を core が足す R と照らす）
    "nominal": "table-residual-nominal-histories.msgpack",
    "one": "table-residual-density-branch.msgpack",
    "one_hybrid": "table-hybrid-residual-density-branch.msgpack",
    "switch": "table-switch-density-branch.msgpack",
    "switch_hybrid": "table-hybrid-switch-branch.msgpack",
    "nominal_half": "table-weight-nominal-0702.msgpack",
    "half": "table-weight-branch-0702.msgpack",
    "half_hybrid": "table-hybrid-weight-0702.msgpack",
    "nominal_full": "table-weight-nominal-0661.msgpack",
    "full": "table-weight-branch-0661.msgpack",
    "full_hybrid": "table-hybrid-weight-0661.msgpack",
}
SUMMARY = "switch.json"
STATE_AXES = ("boron", "moderator", "fuel")
PARTS = ("untracked", "linearization", "residual")
# R(1) と core が足す R（`histories`）の項目の差の上限（同じ式の作り直しなので丸めの程度）
HISTORIES_TOLERANCE = 1e-12
# 密度だけの経路の密度と H⁰ + w (Hᵇ − H⁰) の差の上限 (g/cc)。入力は 3 桁で書く（0.702 = (0.743 + 0.661) / 2）
DENSITY_TOLERANCE = 1e-9


def values(table, index):
    """状態の軸の番号 `index`（ホウ素、減速材、燃料）の値（ホウ素、（減速材の温度、密度）、燃料）"""
    names = [a["name"] for a in table["axes"]]
    if names != ["burnup", *STATE_AXES]:
        raise SystemExit(f"the table has the axes {names}")
    boron, moderator, fuel = table["axes"][1:]
    return (
        boron["points"][index[0]],
        (moderator["points"][index[1]], moderator["density"][index[1]]),
        fuel["points"][index[2]],
    )


def solved(table, condition):
    """テーブルが条件 `condition` を解いた燃焼度ごとの（燃焼度の点の番号、状態の番号）"""
    burnups = table["axes"][0]["points"]
    out = {}
    for p, (b, *index) in enumerate(table["points"]):
        if values(table, index) == condition:
            out[burnups[b]] = (b, p)
    return out


def nominal_state(table, burnup):
    """公称の状態の、燃焼度 `burnup` の（燃焼度の点の番号、状態の番号）"""
    index = [a["points"].index(a["nominal"]) for a in table["axes"][1:]]
    found = solved(table, values(table, index)).get(burnup)
    if found is None:
        raise SystemExit(
            f"the nominal table does not solve the nominal state at {burnup}"
        )
    return found


def split(nominal, branch, hybrid, label):
    """条件の 1 状態の分岐とハイブリッドのテーブルと、同じ条件の状態を持つ公称の組成のテーブルから、条件と、燃焼度の点ごとの U、L、R と、
    条件の線形の補正の後の定数、ΔN̄ᵢ、分岐とハイブリッドの定数"""
    condition, burnups, states = condition_table(branch, f"{label} branch")
    h_condition, h_burnups, h_states = condition_table(hybrid, f"{label} hybrid")
    if h_condition != condition or h_burnups != burnups:
        raise SystemExit(
            f"{label}: the hybrid solves other conditions or burnups than the branch"
        )
    nuclides = nominal["micro"]["nuclides"]
    for table, name in ((branch, "branch"), (hybrid, "hybrid")):
        if table["micro"]["nuclides"] != nuclides:
            raise SystemExit(f"{label}: the {name} tracks other nuclides")
    # H_X の N̄ᵢ は X と同じ（体積の平均。違えば置き換えた核種か履歴が違う）
    if hybrid["micro"]["densities"] != branch["micro"]["densities"]:
        raise SystemExit(f"{label}: the hybrid holds other tracked densities")
    at = solved(nominal, condition)
    rows = []
    for b, (u, p, q) in enumerate(zip(burnups, states, h_states)):
        if u not in at:
            raise SystemExit(f"{label}: the nominal table does not solve {u}")
        bn, pn = at[u]
        deltas = [
            x - y
            for x, y in zip(
                branch["micro"]["densities"][b], nominal["micro"]["densities"][bn]
            )
        ]
        linear = corrected(
            nominal["states"][pn], coefficients(nominal["micro"]["states"][pn]), deltas
        )
        fx, fh = branch["states"][p], hybrid["states"][q]
        rows.append(
            {
                "burnup": u,
                "untracked": difference(fx, fh),
                "linearization": difference(fh, linear),
                "residual": difference(fx, linear),
                "linear": linear,
                "densities": branch["micro"]["densities"][b],
                "branch": fx,
                "hybrid": fh,
            }
        )
    return condition, rows


def base(nominal, burnup, densities):
    """公称の状態の定数に、数密度 `densities` の線形の補正を足したもの"""
    b, p = nominal_state(nominal, burnup)
    deltas = [x - y for x, y in zip(densities, nominal["micro"]["densities"][b])]
    return corrected(
        nominal["states"][p], coefficients(nominal["micro"]["states"][p]), deltas
    )


def one_rows(nominal, rows):
    """1（モジュールの説明）。R(1) は core が足す R（公称のテーブルの `histories`）と同じであること"""
    rod = Table(nominal).rod
    out = []
    for x in rows:
        found = [k for k, u in enumerate(rod["burnup"]) if same_burnup(u, x["burnup"])]
        if len(found) != 1:
            raise SystemExit(f"the histories have no single R at {x['burnup']}")
        check = largest_difference(x["residual"], residual(rod, found[0]))
        if check > HISTORIES_TOLERANCE:
            raise SystemExit(
                f"R(1) at {x['burnup']} differs by {check:.2e} from the R that the core adds"
            )
        state = base(nominal, x["burnup"], x["densities"])
        rho_x, rho_h = rho(k_inf(x["branch"])), rho(k_inf(x["hybrid"]))
        rho_linear = rho(k_inf(x["linear"]))
        out.append(
            {
                "burnup": x["burnup"],
                **{f"{key}_pcm": effect(state, x[key]) for key in PARTS},
                "untracked_condition_pcm": rho_x - rho_h,
                "linearization_condition_pcm": rho_h - rho_linear,
                "residual_condition_pcm": rho_x - rho_linear,
                "histories_check": check,
            }
        )
    return out


def ratio(a, b):
    return None if b == 0.0 else a / b


def switch_rows(nominal, switch, one, summary):
    """2 と 3（モジュールの説明）"""
    by_burnup = {x["burnup"]: x for x in one}
    out = []
    for x in switch:
        u = x["burnup"]
        found = [s for s in summary["rows"] if same_burnup(s["burnup"], u)]
        if len(found) != 1:
            raise SystemExit(f"the switch summary has no single point at {u}")
        (s,) = found
        if u not in by_burnup:
            raise SystemExit(f"R(1) has no point at {u}")
        o = by_burnup[u]
        weight = s["weight"]
        state = base(nominal, u, x["densities"])
        row = {"burnup": u, "weight": weight}
        reference = x["branch"]
        shifted = reference
        for key in ("untracked", "linearization"):
            own = effect(state, x[key])
            whole = effect(state, o[key])
            gap = difference(x[key], scaled(o[key], weight))
            moved = added(reference, scaled(gap, -1.0))
            shifted = added(shifted, scaled(gap, -1.0))
            row |= {
                f"{key}_pcm": own,
                f"scaled_{key}_pcm": effect(state, scaled(o[key], weight)),
                f"one_{key}_pcm": whole,
                # w = 0 の点（続ける点）は成分も w の分も 0 で、比は丸めの値
                f"{key}_weight": None if weight == 0.0 else ratio(own, whole),
                f"{key}_gap_pcm": rho(k_inf(moved)) - rho(k_inf(reference)),
            }
        row |= {
            "residual_pcm": effect(state, x["residual"]),
            "scaled_residual_pcm": effect(state, scaled(o["residual"], weight)),
            "difference_pcm": s["difference_pcm"],
            "composition_pcm": s["composition_pcm"],
            "gap_pcm": rho(k_inf(shifted)) - rho(k_inf(reference)),
        }
        row["remainder_pcm"] = (
            row["difference_pcm"] - row["composition_pcm"] - row["gap_pcm"]
        )
        out.append(row)
    return out


def density_rows(nominal, half, full):
    """4（モジュールの説明）"""
    out = []
    for h, f in zip(half, full, strict=True):
        u = h["burnup"]
        if f["burnup"] != u:
            raise SystemExit("the w = 0.5 and w = 1 tables solve other burnups")
        state = base(nominal, u, h["densities"])
        row = {"burnup": u}
        for key in PARTS:
            own = effect(state, h[key])
            share = effect(state, scaled(f[key], 0.5))
            row |= {
                f"half_{key}_pcm": own,
                f"one_half_{key}_pcm": share,
                f"{key}_difference_pcm": own - share,
                f"{key}_ratio": None if u == 0.0 else ratio(own, 2.0 * share),
            }
        out.append(row)
    return out


def check_nominal(nominal, path_nominal, label):
    """密度だけの経路の公称の組成のテーブルが、公称のテーブルと同じ履歴の組成（ΔN̄ᵢ の基準が同じ）であること"""
    _, burnups, _ = condition_table(path_nominal, f"{label} nominal")
    for b, u in enumerate(burnups):
        bn, _ = nominal_state(nominal, u)
        if path_nominal["micro"]["densities"][b] != nominal["micro"]["densities"][bn]:
            raise SystemExit(
                f"{label}: the nominal composition differs from the nominal table at {u}"
            )


def history_condition(nominal):
    """公称のテーブルの公称の状態と、`histories` の条件 j（減速材の軸だけが Hᵇ の格子点の状態）と、H⁰、Hᵇ"""
    t = Table(nominal)
    moderator = nominal["axes"][2]
    index = list(t.nominal)
    index[1] = moderator["density"].index(t.history["branch"])
    return (
        values(nominal, t.nominal),
        values(nominal, index),
        (t.history["nominal"], t.history["branch"]),
    )


def compare(tables, summary):
    """要約（モジュールの説明）"""
    nominal = tables["nominal"]
    nominal_condition, condition_j, (h0, hb) = history_condition(nominal)
    one_condition, one = split(nominal, tables["one"], tables["one_hybrid"], "R(1)")
    switch_condition, switch = split(
        nominal, tables["switch"], tables["switch_hybrid"], "switch"
    )
    if one_condition != condition_j or switch_condition != condition_j:
        raise SystemExit(
            "R(1) and the switch must be at the condition j of the histories"
        )
    check_nominal(nominal, tables["nominal_half"], "w = 0.5")
    check_nominal(nominal, tables["nominal_full"], "w = 1")
    half_condition, half = split(
        tables["nominal_half"], tables["half"], tables["half_hybrid"], "w = 0.5"
    )
    full_condition, full = split(
        tables["nominal_full"], tables["full"], tables["full_hybrid"], "w = 1"
    )
    # 密度だけの経路: 公称の状態と同じホウ素、減速材の温度、燃料で、密度が w = 0.5 と 1 の値（H⁰ + w (Hᵇ − H⁰)）
    for c, label, weight in (
        (half_condition, "w = 0.5", 0.5),
        (full_condition, "w = 1", 1.0),
    ):
        if (c[0], c[1][0], c[2]) != (
            nominal_condition[0],
            nominal_condition[1][0],
            nominal_condition[2],
        ):
            raise SystemExit(
                f"{label}: the condition changes more than the moderator density"
            )
        if abs(c[1][1] - (h0 + weight * (hb - h0))) > DENSITY_TOLERANCE:
            raise SystemExit(f"{label}: the moderator density {c[1][1]} is not at w")
    return {
        "format": "xlc-compare-history-residual-split",
        "continued_burnup": summary["continued_burnup"],
        "one": one_rows(nominal, one),
        "switch": switch_rows(nominal, switch, one, summary),
        "density": density_rows(nominal, half, full),
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(values, form="+.2f"):
        values = [v for v in values if v is not None]
        return f"{min(values):{form}}〜{max(values):{form}}"

    one = summary["one"]
    last = one[-1]
    lines = [
        (
            f"R(1): 追わない核種の分 U の効き {span([x['untracked_pcm'] for x in one])} pcm、線形化の誤差 L "
            f"{span([x['linearization_pcm'] for x in one])} pcm、R {span([x['residual_pcm'] for x in one])} pcm（公称の状態の定数）。"
            f"最後の点 {last['burnup']:.4g} MWd/kg で U {last['untracked_pcm']:+.2f}、L {last['linearization_pcm']:+.2f}、"
            f"R {last['residual_pcm']:+.2f} pcm。条件 j での効きは U {last['untracked_condition_pcm']:+.2f}、L "
            f"{last['linearization_condition_pcm']:+.2f}、R {last['residual_condition_pcm']:+.2f} pcm。"
            f"core が足す R（`histories`）との項目の差は最大 {max(x['histories_check'] for x in one):.2e}"
        )
    ]
    u0 = summary["continued_burnup"]
    later = [x for x in summary["switch"] if x["burnup"] > u0]
    s = later[-1]
    lines.append(
        f"switch（{u0:.4g} MWd/kg より後の点 {len(later)}）: w_U {span([x['untracked_weight'] for x in later], '.4f')}、"
        f"w_L {span([x['linearization_weight'] for x in later], '.4f')}（w {span([x['weight'] for x in later], '.4f')}）。"
        f"最後の点 {s['burnup']:.4g} MWd/kg で w {s['weight']:.4f}、U {s['untracked_pcm']:+.2f}（w U(1) {s['scaled_untracked_pcm']:+.2f}）、"
        f"L {s['linearization_pcm']:+.2f}（w L(1) {s['scaled_linearization_pcm']:+.2f}）pcm"
    )
    lines.append(
        f"core の差 {span([x['difference_pcm'] for x in later], '+.3f')} pcm のうち、U の外れの分 "
        f"{span([x['untracked_gap_pcm'] for x in later])}、L の外れの分 {span([x['linearization_gap_pcm'] for x in later])}、"
        f"両方 {span([x['gap_pcm'] for x in later])} pcm、数密度の効き {span([x['composition_pcm'] for x in later], '+.3f')} pcm、"
        f"残り {span([x['remainder_pcm'] for x in later], '+.3f')} pcm。最後の点で差 {s['difference_pcm']:+.3f} = 両方の外れ "
        f"{s['gap_pcm']:+.3f} + 数密度の効き {s['composition_pcm']:+.3f} + 残り {s['remainder_pcm']:+.3f} pcm（外れは反応度で"
        f"足し合わず、U の外れだけ {s['untracked_gap_pcm']:+.3f}、L の外れだけ {s['linearization_gap_pcm']:+.3f} pcm）"
    )
    density = summary["density"]
    d = density[-1]
    lines.append(
        "密度だけの経路の X(0.5) − X(1)/2 の効き（w = 0.5 の定数）: "
        + "、".join(
            f"{name} {span([x[f'{key}_difference_pcm'] for x in density])} pcm（比 "
            f"{span([x[f'{key}_ratio'] for x in density], '.4f')}）"
            for key, name in (
                ("untracked", "追わない核種の分"),
                ("linearization", "線形化の誤差"),
                ("residual", "R"),
            )
        )
        + f"。最後の点 {d['burnup']:.4g} MWd/kg で U {d['untracked_difference_pcm']:+.2f}、L "
        f"{d['linearization_difference_pcm']:+.2f}、R {d['residual_difference_pcm']:+.2f} pcm"
    )
    # 記録の表（verification_tables.py の signed）と同じ負の符号
    return "\n".join(lines).replace("-", "−")


def check_summary(summary, md5s):
    """`switch.json` が、今の参照のテーブルと、core が使った R と係数を持つ今の公称のテーブルから作られたこと"""
    for key in ("switch", "nominal"):
        file = FILES[key]
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
