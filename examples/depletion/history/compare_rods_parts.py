"""制御棒を途中まで入れて抜く燃焼（ADR-0215 の確認の制御棒の履歴、ADR-0228、#0190）の core − 参照を、追う核種の数密度の効きと定数の
側に分け、定数の側を U と L の重みと、制御棒の履歴を指数で忘れる重みで見る。

    python3 compare_rods_parts.py --json rods-parts.json     # このディレクトリの出力から JSON を書く
    python3 compare_rods_parts.py --report rods-parts.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは compare_rods_switch.py と同じ、2 つの変数の残りを持つ公称のテーブル、格子計算の参照（2A-rods の履歴の 20 MWd/kg の点から
制御棒を抜いて続けた燃焼の履歴）、履歴の補正の core の燃焼の履歴と、その比較の要約（`rods-switch.json`。md5 と core − 参照を照らす）。

参照の燃焼点のうち切り替えの点より後ごとに、core のノードの定数を core の `node_constants` と同じ式で組み立て直す:

    f_core = f(公称、out) + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰) + w U(1) + (P / P_j) L(1)

f(公称、out)、係数 cᵢ（数密度の微分を持つ核種は微分）、N̄ᵢ⁰ はテーブルの制御棒の値 out の公称の状態を、U(1)、L(1)、ΔN̄ᵢ^j は制御棒の
変数の値 out の残りを、core と同じ燃焼度の内挿（`burnup_terms`）で引く。w = (H_r − H⁰) / (Hᵇ − H⁰)、P / P_j は熱群の νΣf の
線形の補正の比（ノードの数密度 Nᵢ と履歴の分岐の数密度の差）。減速材の密度の履歴は公称（w = 0）で、core はその補正を足さない
（公称でなければ止める）。再構成 ρ(core の k) − ρ(f_core) が `RECONSTRUCT_TOLERANCE` を超えれば止める。

参照の定数（参照の組成を制御棒を抜いた状態で解いた少数群定数）は格子計算の燃焼の履歴に無いので、効きは core のノードの定数で測る
（ρ = 1 − 1/k の差の pcm、k は 2 群の無限媒質の k）:

- 数密度の効き ρ(f_core) − ρ(f_ref)。f_ref は、同じ式のノードの定数を参照の均質化した数密度 N̄ᵢ で組み立て直したもの（線形の補正と
  P / P_j の両方が N̄ᵢ のもの。P / P_j の参照の組成の値も並べる）。核種ごとの効き（その核種だけを参照の数密度にして組み立て直したもの）と
  数密度の相対差も並べる
- 定数の側 = core − 参照 − 数密度の効き（= ρ(f_ref) − ρ(参照の k∞) + 再構成）。参照の組成での R_switch と core の R の違いと、
  2 群の k と格子計算の k∞ の違いを含む
- w U(1) の効き ρ(f_ref) − ρ(f_ref − w U(1))、(P / P_j) L(1) の効きも同じ形（P / P_j は参照の組成の値）。それを引いた定数の側（見込み。
  core を流していない）
- 参照に合う重み: f_ref の U(1) の重みだけを変えて定数の側が 0 になる重み（L の重みも同じ）。`WEIGHT_RANGE` で探し、無ければ None

U_switch と L_switch（参照の組成の U と L）は、参照の組成とそのハイブリッドを制御棒を抜いた状態で解いたテーブルが要るので出さない。

指数で忘れる重み: 制御棒の履歴を H_τ（compare_forgetting.py の H_τ）にすると、燃焼度 B₀ まで制御棒を入れ、その後を抜いた履歴の重みは
1 − compare_forgetting.weight(B, B₀, τ)（τ → ∞ が今の core の w = B₀ / B）。τ ごとに、f_ref の U(1) の重みを w から w_τ にした定数の側
ρ(f_ref + (w_τ − w) U(1)) − ρ(参照の k∞) + 再構成（L と core の数密度は今のまま）と、その絶対値の最大。最大が最小の τ を `best` にし、
それが 0 の側の隣の τ（正の最小か負の最大）なら `best_at_edge`（その先を比べておらず、最小を挟んでいない）。
"""

import argparse
import json
import re
from pathlib import Path

import compare_forgetting as forgetting
from compare_core_history import same_burnup
from compare_residual import added, coefficients, corrected, k_inf, md5, read, rho
from compare_rods_switch import (
    VARIABLES,
    entry,
    homogenized,
    nodes,
    percent,
    production,
    scaled,
)

FILES = {
    "table": "table-residual-rods-histories.msgpack",
    "reference": "lattice-switch-rods.msgpack",
    "core": "core-2A-rods-switch.msgpack",
}
SUMMARY = "rods-switch.json"
# compare_forgetting.py の τ に、10 MWd/kg より小さい正の τ を足す（制御棒の経路は 10 でもまだ縮むので、最小を挟むため）
TAUS = [*forgetting.TAUS, 7, 5, 3, 2, 1]
ROD = "out"
# core の k と組み立て直した定数の 2 群の k の差の上限 (pcm)。同じ式なので丸めの程度（2A の 15 点で 4.1e−06 pcm）
RECONSTRUCT_TOLERANCE = 1e-3
# rods-switch.json の core − 参照と、ここで作った値の差の上限 (pcm)（同じ入力と式なので丸めの程度）
SUMMARY_TOLERANCE = 1e-9
# 内挿する項目: 少数群定数は 2 群の k と線形の補正、残りの足し算が使うもの（compare_residual の `corrected`、`added`）、係数は
# `corrected` が使うもの。テーブルのほかの項目（vacuum_ratio、pins、反応の一覧など）は値の無いものを含み、k に効かない
STATE_KEYS = (
    "total",
    "absorption",
    "nu_fission",
    "fission",
    "chi",
    "diffusion",
    "scatter",
    "adf",
)
FACTOR_KEYS = ("total", "nu_fission", "fission", "scatter")
# 参照に合う重みを探す範囲
WEIGHT_RANGE = (-2.0, 2.0)


def bracket(points, x):
    """線形の内挿の（点の番号、重み）。格子点ではその点だけ。範囲の外はエラー"""
    for i in range(len(points) - 1):
        a, b = points[i], points[i + 1]
        if a <= x <= b:
            t = (x - a) / (b - a)
            if t == 0.0:
                return [(i, 1.0)]
            if t == 1.0:
                return [(i + 1, 1.0)]
            return [(i, 1.0 - t), (i + 1, t)]
    raise SystemExit(f"the burnup {x} is outside the points {points}")


def slope(points, k):
    """点 `k` の傾きの、3 点の値に掛ける係数（点の番号と係数。core の `slope`）"""
    s = min(max(k - 1, 0), len(points) - 3)
    p = points[s : s + 3]
    x = points[k]

    def coefficient(m):
        a, b = (p[j] for j in range(3) if j != m)
        return (2.0 * x - a - b) / ((p[m] - a) * (p[m] - b))

    return [(s + m, coefficient(m)) for m in range(3)]


def burnup_terms(points, x):
    """燃焼度 `x` の内挿の（点の番号、重み）。core の `burnup_terms`（ADR-0192）の写し: 3 点の放物線の傾きの 3 次 Hermite で、燃焼度 0 の
    点は曲線に使わず、0 から次の点までと、0 を除いて 3 点に満たない軸は線形"""
    linear = bracket(points, x)
    first = 1 if points[0] == 0.0 else 0
    if len(linear) != 2:
        return linear
    (i, _), (_, t) = linear
    if i < first or len(points) - first < 3:
        return linear
    smooth = points[first:]
    j = i - first
    h = smooth[j + 1] - smooth[j]
    h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t)
    h01 = t * t * (3.0 - 2.0 * t)
    h10 = t * (1.0 - t) * (1.0 - t)
    h11 = -t * t * (1.0 - t)
    base = max(j - 1, 0)
    weights = [0.0] * 4
    weights[j - base] += h00
    weights[j + 1 - base] += h01
    for k, hk in ((j, h10), (j + 1, h11)):
        for p, coeff in slope(smooth, k):
            weights[p - base] += h * hk * coeff
    return [(first + base + m, w) for m, w in enumerate(weights) if w != 0.0]


def combined(items, keys=None):
    """（重み、値）の組の和 Σ w x。値は数か、数の入れ子のリスト。`keys` を渡すと値は dict で、その項目だけを足した dict"""
    if keys is not None:
        return {key: combined([(w, x[key]) for w, x in items]) for key in keys}
    _, first = items[0]
    if isinstance(first, list):
        return [combined([(w, x[i]) for w, x in items]) for i in range(len(first))]
    return sum(w * x for w, x in items)


def state_point(table, b):
    """燃焼度の点の番号 `b` の、制御棒の値 out の公称の状態（連続の軸は公称）の状態の番号"""
    point = []
    for a in table["axes"]:
        if a["name"] == "burnup":
            point.append(b)
        elif a["name"] == "rod":
            point.append(a["values"].index(ROD))
        else:
            point.append(a["points"].index(a["nominal"]))
    return table["points"].index(point)


def history_weights(table, history):
    """core の点の履歴の変数の値（`VARIABLES` の順）の、変数ごとの w = (H − H⁰) / (Hᵇ − H⁰)"""
    out = []
    for name, h in zip(VARIABLES, history, strict=True):
        (v,) = [x for x in table["histories"] if x["variable"] == name]
        out.append((h - v["nominal"]) / (v["branch"] - v["nominal"]) + 0.0)
    return out


def node(table, burnup, densities, weight):
    """core のノードの定数 f_core と、係数、U(1)、L(1)、P / P_j（モジュールの説明）"""
    axis = table["axes"][0]["points"]
    micro = table["micro"]
    terms = burnup_terms(axis, burnup)
    state = combined(
        [(w, table["states"][state_point(table, b)]) for b, w in terms], STATE_KEYS
    )
    by_point = [
        (w, coefficients(micro["states"][state_point(table, b)])) for b, w in terms
    ]
    factors = [
        combined([(w, c[i]) for w, c in by_point], FACTOR_KEYS)
        for i in range(len(micro["nuclides"]))
    ]
    nominal = combined([(w, micro["densities"][b]) for b, w in terms])
    rod = entry(table, VARIABLES[1], ROD)
    r_terms = burnup_terms(rod["burnup"], burnup)
    untracked = combined([(w, rod["untracked"][k]) for k, w in r_terms], STATE_KEYS)
    linearization = combined(
        [(w, rod["linearization"][k]) for k, w in r_terms], STATE_KEYS
    )
    branch = combined(
        [
            (
                w,
                [
                    x - y
                    for x, y in zip(
                        rod["densities"][k],
                        micro["densities"][axis.index(rod["burnup"][k])],
                    )
                ],
            )
            for k, w in r_terms
        ]
    )
    deltas = [x - y for x, y in zip(densities, nominal)]
    of_branch = production(factors, branch)
    if of_branch == 0.0:
        raise SystemExit(f"the branch index is 0 at {burnup}")
    index = production(factors, deltas) / of_branch
    f = corrected(state, factors, deltas)
    f = added(added(f, scaled(untracked, weight)), scaled(linearization, index))
    return f, factors, untracked, linearization, index


def removed_effect(f, part):
    """定数 `f` から `part` を引いた効き ρ(f) − ρ(f − part) (pcm)"""
    return rho(k_inf(f)) - rho(k_inf(added(f, scaled(part, -1.0))))


def reference_weight(base, part, now, target):
    """`base` の `part` の重みを `now` から変えて ρ が `target` に合う重み（`WEIGHT_RANGE` の二分法。無ければ None）"""

    def gap(x):
        return rho(k_inf(added(base, scaled(part, x - now)))) - target

    lo, hi = WEIGHT_RANGE
    g_lo = gap(lo)
    if g_lo * gap(hi) > 0.0:
        return None
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        g_mid = gap(mid)
        if g_lo * g_mid <= 0.0:
            hi = mid
        else:
            lo, g_lo = mid, g_mid
    return 0.5 * (lo + hi)


def compare(table, reference, core, taus=TAUS):
    """要約（モジュールの説明）"""
    nuclides = table["micro"]["nuclides"]
    names = [h["variable"] for h in table["histories"]]
    if names != VARIABLES:
        raise SystemExit(f"the table holds the history variables {names}")
    continued = reference.get("continued")
    if continued is None:
        raise SystemExit("the reference is not a continued history")
    switch = continued["burnup"]
    ref = homogenized(reference, nuclides, table["micro"]["densities"][0])
    core_nodes = nodes(core, "core", nuclides)
    rows = []
    bases = []
    for r in ref:
        if r["burnup"] <= switch:
            continue
        found = [n for n in core_nodes if same_burnup(n["burnup"], r["burnup"])]
        if len(found) != 1:
            raise SystemExit(f"core: no single point at burnup {r['burnup']}")
        (n,) = found
        if n["state"] != ROD:
            raise SystemExit(f"core: the node at {r['burnup']} is not in state {ROD}")
        density_weight, weight = history_weights(table, n["history"])
        if density_weight != 0.0:
            raise SystemExit(
                f"core: the moderator density history at {r['burnup']} is not nominal"
            )
        f, _, untracked, linearization, index = node(
            table, n["burnup"], n["densities"], weight
        )
        reconstruct = rho(n["k"]) - rho(k_inf(f))
        if abs(reconstruct) > RECONSTRUCT_TOLERANCE:
            raise SystemExit(
                f"core: the k at {r['burnup']} differs by {reconstruct} pcm from its constants"
            )
        # P / P_j もノードの数密度で決まるので、参照の数密度でノードの定数を組み立て直す（線形の補正だけを動かすと、P / P_j の変わる分が
        # 定数の側に入る）
        at_reference, _, _, _, reference_index = node(
            table, n["burnup"], r["densities"], weight
        )
        u = scaled(untracked, weight)
        lin = scaled(linearization, reference_index)
        target = rho(r["k_inf"])
        difference = rho(n["k"]) - target
        composition = rho(k_inf(f)) - rho(k_inf(at_reference))
        # 定数の側を 0 にする ρ（定数の側 = ρ(f_ref) − ρ(参照の k∞) + 再構成）
        zero = target - reconstruct
        per_nuclide = []
        for i, name in enumerate(nuclides):
            only = list(n["densities"])
            only[i] = r["densities"][i]
            g, _, _, _, _ = node(table, n["burnup"], only, weight)
            per_nuclide.append(
                {
                    "nuclide": name,
                    "percent": None
                    if r["densities"][i] == 0.0
                    else percent(n["densities"][i], r["densities"][i]),
                    "pcm": rho(k_inf(f)) - rho(k_inf(g)),
                }
            )
        rows.append(
            {
                "burnup": r["burnup"],
                "weight": weight,
                "index_ratio": index,
                "reference_index_ratio": reference_index,
                "reference_k": r["k_inf"],
                "difference_pcm": difference,
                "reconstruct_pcm": reconstruct,
                "composition_pcm": composition,
                "constant_pcm": difference - composition,
                "untracked_pcm": removed_effect(at_reference, u),
                "linearization_pcm": removed_effect(at_reference, lin),
                "without_untracked_pcm": rho(
                    k_inf(added(at_reference, scaled(u, -1.0)))
                )
                - zero,
                "without_linearization_pcm": rho(
                    k_inf(added(at_reference, scaled(lin, -1.0)))
                )
                - zero,
                "untracked_weight": reference_weight(
                    at_reference, untracked, weight, zero
                ),
                "linearization_weight": reference_weight(
                    at_reference, linearization, reference_index, zero
                ),
                "nuclides": per_nuclide,
            }
        )
        bases.append((r["burnup"], weight, at_reference, untracked, zero))
    scan = []
    for tau in taus:
        points = []
        for burnup, weight, at_reference, untracked, zero in bases:
            w = 1.0 - forgetting.weight(burnup, switch, tau)
            moved = added(at_reference, scaled(untracked, w - weight))
            points.append(
                {
                    "burnup": burnup,
                    "weight": w,
                    "constant_pcm": rho(k_inf(moved)) - zero,
                }
            )
        scan.append(
            {
                "tau": tau,
                "points": points,
                "largest_pcm": max(
                    (abs(p["constant_pcm"]) for p in points), default=0.0
                ),
            }
        )
    best = min(scan, key=lambda x: x["largest_pcm"])["tau"]
    # 0 の側の隣の τ（正の最小、負の最大）が最小なら、その先を比べておらず最小を挟んでいない
    positive = [t for t in taus if t is not None and t > 0]
    negative = [t for t in taus if t is not None and t < 0]
    edges = {min(positive, default=None), max(negative, default=None)} - {None}
    return {
        "format": "xlc-compare-rods-parts",
        "switch": switch,
        "rows": rows,
        "best": best,
        "best_at_edge": best in edges,
        "forgetting": scan,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(rows, key, form="+.3f"):
        values = [r[key] for r in rows if r[key] is not None]
        return f"{min(values):{form}}〜{max(values):{form}}"

    def name(tau):
        return "∞" if tau is None else f"{tau:g}"

    rows = summary["rows"]
    largest = max(
        (x for r in rows for x in r["nuclides"] if x["percent"] is not None),
        key=lambda x: abs(x["percent"]),
    )
    at = next(r for r in rows if largest in r["nuclides"])
    lines = [
        (
            f"切り替えの後（{summary['switch']:g} MWd/kg より後の {len(rows)} 点）の core − 参照 {span(rows, 'difference_pcm')} pcm のうち、"
            f"数密度の効き {span(rows, 'composition_pcm')} pcm、定数の側 {span(rows, 'constant_pcm')} pcm"
            f"（再構成 {span(rows, 'reconstruct_pcm', '+.1e')} pcm）"
        ),
        (
            f"数密度の相対差の絶対値の最大は {largest['nuclide']} の {at['burnup']:g} MWd/kg の {largest['percent']:+.3f}%。"
            "点ごとの核種の効きの絶対値の最大: "
            + "、".join(
                f"{r['burnup']:g} {max(r['nuclides'], key=lambda x: abs(x['pcm']))['nuclide']} "
                f"{max(r['nuclides'], key=lambda x: abs(x['pcm']))['pcm']:+.2f}"
                for r in rows
            )
            + " pcm"
        ),
        (
            f"w U(1) の効き {span(rows, 'untracked_pcm', '+.2f')} pcm、(P / P_j) L(1) の効き {span(rows, 'linearization_pcm', '+.2f')} pcm。"
            f"U を引いた定数の側 {span(rows, 'without_untracked_pcm')} pcm、L を引いた定数の側 {span(rows, 'without_linearization_pcm')} pcm"
        ),
        (
            f"参照に合う U の重み {span(rows, 'untracked_weight', '+.4f')}（w {span(rows, 'weight', '.4f')}）、"
            f"L の重み {span(rows, 'linearization_weight', '+.4f')}（参照の組成の P / P_j {span(rows, 'reference_index_ratio', '.4f')}、"
            f"core の P / P_j {span(rows, 'index_ratio', '.4f')}）"
        ),
    ]
    scan = summary["forgetting"]
    best = next(x for x in scan if x["tau"] == summary["best"])
    infinite = next(x for x in scan if x["tau"] is None)
    values = [p["constant_pcm"] for p in best["points"]]
    lines.append(
        f"U の重みを指数で忘れる w_τ にした定数の側の絶対値の最大が最小の τ は {name(best['tau'])} MWd/kg で "
        f"{min(values):+.3f}〜{max(values):+.3f} pcm（最大 {best['largest_pcm']:.3f} pcm。∞ は {infinite['largest_pcm']:.3f} pcm）"
        + (
            "。この τ は比べた列の 0 の側の端で、最小を挟んでいない。"
            if summary["best_at_edge"]
            else "。"
        )
        + "τ ごとの最大 (pcm): "
        + "、".join(f"{name(x['tau'])} {x['largest_pcm']:.3f}" for x in scan)
    )
    # 核種の名前の - は残す
    return re.sub(r"(?<![A-Za-z])-", "−", "\n".join(lines))


def check_summary(summary, md5s, result):
    """`rods-switch.json` が同じ入力から作られ、core − 参照が同じであること"""
    for file in FILES.values():
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
