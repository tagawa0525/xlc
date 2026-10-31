"""履歴の残り R の線形化の誤差 L と追わない核種の分 U を、組成の指標で表せるかを比べる（ADR-0215、#0190）。

    python3 compare_index.py --json index.json     # このディレクトリの出力から JSON を書く
    python3 compare_index.py --report index.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

compare_split.py と同じテーブルを読み、同じ式で履歴 X ごとに U_X と L_X を作る（4 本: 一生を条件 j の R(1)、途中で条件を変えた燃焼の
switch、密度だけの経路の w = 0.5 と 1）。組成の指標は 2 つ:

- E_X（線形の補正の効き）: 公称のテーブルの公称の状態（565 K、0.743 g/cc）の定数に、その履歴の数密度の線形の補正 Σᵢ cᵢ⁰ ΔN̄ᵢ^X を
  足した定数の ρ − 公称の状態の定数の ρ（pcm）。U と L の効きは、compare_split.py と同じく、その線形の補正を足した定数に足した効き。
  4 本で同じ定数（公称の状態）から測るので、同じ燃焼度で並べられる
- 主な核種の相対差 ΔN̄ᵢ^X / N̄ᵢ⁰（N̄ᵢ⁰ が 0 の核種は出さない）

燃焼度の点ごとに、L/E（1 次の比）、L/E²（2 次の比）、U/E、U/w と、それぞれの R(1) の同じ燃焼度の値との比を出す（L がその次数で指標に
比例するなら、比は履歴によらず 1）。E が 0 の点（公称の組成）は比を出さない。

switch の参照の点ごとに、R(1) の成分を何倍して足すかで外れ δ = X_switch − a X(1) を作り、その分を ρ(参照の定数 − δ) − ρ(参照の定数) と
する（compare_split.py の外れの分と同じ式。a = w で compare_split.py の値と一致することを確かめる）。a は次の 4 つ:

- w（core の今の置き方。`switch.json` の `weight`）
- 指標の比 E_switch / E(1)（1 次）と、その 2 乗（2 次）
- 条件 j での指標の比: ρ(f(公称, j) + Σᵢ cᵢ(j) ΔN̄ᵢ) − ρ(f(公称, j)) の switch と R(1) の比（L のほうだけ。指標を測る定数の選び方の確認）

ほかに、U を w、L を指標の比で足したときの両方の外れの分（compare_split.py の両方の外れと同じ式）と、主な核種の ΔN̄ᵢ の比
ΔN̄ᵢ^switch / ΔN̄ᵢ^(1) を出す。

E は燃焼度の前半で符号を変えるので、その近くでは比を指標に使えない。符号を変えない候補として、群ごとの線形の補正の量の指標を 2 つ
比べる（公称の状態の係数 cᵢ⁰ で、熱群（最後の群）の値。1/cm）:

- P（熱群の νΣf の補正）: Σᵢ cᵢ⁰_νf,G ΔN̄ᵢ
- A（熱群の吸収の補正）: Σᵢ (cᵢ⁰_t,G − Σ_g' cᵢ⁰_s,G→g') ΔN̄ᵢ

switch の参照の点ごとに、L を P と A の比で R(1) から内挿した外れの分と、U を w、L を P の比で足した両方の外れの分を出す。密度だけの
経路の 2 本（w = 0.5 と 1）は、燃焼度の点ごとに、E、P、A、L の比（0.702 / 0.661）を並べる（L は 2 本でほぼ同じ比に乗るので、
指標の比がそこから外れる点では、その指標で L を内挿できない）。
"""

import argparse
import json
import re
from pathlib import Path

import compare_split as split
from compare_core_history import same_burnup
from compare_residual import added, coefficients, difference, k_inf, md5, read, rho
from compare_weight import effect, scaled

SUMMARIES = ("switch.json", "split.json")
INPUTS = (*split.FILES.values(), *SUMMARIES)
NUCLIDES = ("U-235", "Pu-239", "Pu-240", "Pu-241")
# 履歴の名前（JSON の鍵）と、R(1) の w に対する重み（switch は点ごとに `switch.json` の w）
HISTORIES = ("one", "switch", "half", "full")
WEIGHTS = {"one": 1.0, "half": 0.5, "full": 1.0}
# 外れの分を w で作り直した値と split.json の値の差の上限 (pcm)。同じ式と同じテーブルなので丸めの程度
SPLIT_TOLERANCE = 1e-9


def ratio(a, b):
    return None if a is None or b is None or b == 0.0 else a / b


def index_effect(nominal, burnup, densities):
    """指標 E: 公称の状態の定数に数密度 `densities` の線形の補正を足した効き (pcm)"""
    _, p = split.nominal_state(nominal, burnup)
    return rho(k_inf(split.base(nominal, burnup, densities))) - rho(
        k_inf(nominal["states"][p])
    )


def group_indices(nominal, burnup, densities):
    """熱群（最後の群）の線形の補正の量 (P, A)（モジュールの説明）。係数は公称の状態の値"""
    b, p = split.nominal_state(nominal, burnup)
    deltas = [x - y for x, y in zip(densities, nominal["micro"]["densities"][b])]
    factors = coefficients(nominal["micro"]["states"][p])
    production = sum(c["nu_fission"][-1] * d for c, d in zip(factors, deltas))
    absorption = sum(
        (c["total"][-1] - sum(c["scatter"][-1])) * d for c, d in zip(factors, deltas)
    )
    return production, absorption


def positions(nominal, names):
    nuclides = nominal["micro"]["nuclides"]
    missing = [n for n in names if n not in nuclides]
    if missing:
        raise SystemExit(f"the tables do not track {missing}")
    return {n: nuclides.index(n) for n in names}


def relative_densities(nominal, burnup, densities, at):
    """核種ごとの ΔN̄ᵢ / N̄ᵢ⁰（N̄ᵢ⁰ が 0 なら None）"""
    b, _ = split.nominal_state(nominal, burnup)
    n0 = nominal["micro"]["densities"][b]
    return {n: ratio(densities[i] - n0[i], n0[i]) for n, i in at.items()}


def history_rows(nominal, rows, weights, at):
    """履歴の燃焼度の点ごとの指標と成分と比。`weights` は点ごとの w"""
    out = []
    for x, w in zip(rows, weights, strict=True):
        u = x["burnup"]
        state = split.base(nominal, u, x["densities"])
        e = index_effect(nominal, u, x["densities"])
        linearization = effect(state, x["linearization"])
        untracked = effect(state, x["untracked"])
        production, absorption = group_indices(nominal, u, x["densities"])
        out.append(
            {
                "burnup": u,
                "weight": w,
                "index_pcm": e,
                "production": production,
                "thermal_absorption": absorption,
                "linearization_pcm": linearization,
                "untracked_pcm": untracked,
                "linearization_first": ratio(linearization, e),
                "linearization_second": ratio(
                    linearization, None if e == 0.0 else e * e
                ),
                "untracked_first": ratio(untracked, e),
                "untracked_weight": ratio(untracked, w),
                "densities": relative_densities(nominal, u, x["densities"], at),
            }
        )
    return out


RATIOS = (
    "linearization_first",
    "linearization_second",
    "untracked_first",
    "untracked_weight",
)


def to_one(histories):
    """各比の、R(1) の同じ燃焼度の値との比を書き足す"""
    one = histories["one"]
    for rows in histories.values():
        for x in rows:
            found = [o for o in one if same_burnup(o["burnup"], x["burnup"])]
            if len(found) != 1:
                raise SystemExit(f"R(1) has no single point at {x['burnup']}")
            for key in RATIOS:
                x[f"{key}_to_one"] = ratio(x[key], found[0][key])


def removed(reference, delta):
    """参照の定数から外れ `delta` を引いた効き (pcm)"""
    return rho(k_inf(added(reference, scaled(delta, -1.0)))) - rho(k_inf(reference))


def gap(reference, own, whole, factor):
    """外れ δ = own − factor × whole を参照の定数から引いた効き (pcm)。factor が無ければ None"""
    if factor is None:
        return None
    return removed(reference, difference(own, scaled(whole, factor)))


def switch_rows(nominal, condition, switch, one, summary, split_summary, at):
    """switch の参照の点ごとの、指標の比と、成分を w と指標の比で R(1) から内挿した外れの分"""
    by_burnup = {x["burnup"]: x for x in one}
    states = split.solved(nominal, condition)
    out = []
    for x in switch:
        u = x["burnup"]
        rows = [s for s in summary["rows"] if same_burnup(s["burnup"], u)]
        checks = [s for s in split_summary["switch"] if same_burnup(s["burnup"], u)]
        if len(rows) != 1 or len(checks) != 1 or u not in by_burnup:
            raise SystemExit(f"the summaries or R(1) have no single point at {u}")
        (s,), (check,), o = rows, checks, by_burnup[u]
        w = s["weight"]
        own_index = index_effect(nominal, u, x["densities"])
        one_index = index_effect(nominal, u, o["densities"])
        # R(1) が公称の組成の点（燃焼度 0）では比を出さない
        index = ratio(own_index, one_index)
        own_groups = group_indices(nominal, u, x["densities"])
        one_groups = group_indices(nominal, u, o["densities"])
        production = ratio(own_groups[0], one_groups[0])
        absorption = ratio(own_groups[1], one_groups[1])
        j = rho(k_inf(nominal["states"][states[u][1]]))
        condition_index = ratio(
            rho(k_inf(x["linear"])) - j, rho(k_inf(o["linear"])) - j
        )
        reference = x["branch"]
        row = {
            "burnup": u,
            "weight": w,
            "difference_pcm": s["difference_pcm"],
            "index_ratio": index,
            "condition_index_ratio": condition_index,
            "production_ratio": production,
            "absorption_ratio": absorption,
            "linearization_weight": check["linearization_weight"],
            "untracked_weight": check["untracked_weight"],
            "density_ratios": {
                n: ratio(
                    x["densities"][i] - nominal_density(nominal, u, i),
                    o["densities"][i] - nominal_density(nominal, u, i),
                )
                for n, i in at.items()
            },
        }
        factors = {
            "weight": w,
            "index": index,
            "index_squared": None if index is None else index * index,
        }
        for key in ("untracked", "linearization"):
            for name, factor in factors.items():
                row[f"{key}_gap_{name}_pcm"] = gap(reference, x[key], o[key], factor)
            worst = abs(row[f"{key}_gap_weight_pcm"] - check[f"{key}_gap_pcm"])
            if worst > SPLIT_TOLERANCE:
                raise SystemExit(
                    f"the {key} gap at w differs from split.json by {worst:.2e} at {u}"
                )
        row["linearization_gap_condition_index_pcm"] = gap(
            reference, x["linearization"], o["linearization"], condition_index
        )
        for name, factor in (("production", production), ("absorption", absorption)):
            row[f"linearization_gap_{name}_pcm"] = gap(
                reference, x["linearization"], o["linearization"], factor
            )
        row["gap_production_pcm"] = mixed(reference, x, o, w, production)
        row["gap_mixed_pcm"] = mixed(reference, x, o, w, index)
        row["gap_weight_pcm"] = check["gap_pcm"]
        out.append(row)
    return out


def mixed(reference, own, one, w, factor):
    """U を w、L を `factor` で R(1) から内挿した両方の外れの分 (pcm)。factor が無ければ None"""
    if factor is None:
        return None
    return removed(
        reference,
        added(
            difference(own["untracked"], scaled(one["untracked"], w)),
            difference(own["linearization"], scaled(one["linearization"], factor)),
        ),
    )


def density_pairs(half, full):
    """密度だけの経路の燃焼度の点ごとの、E、P、A、L の比（0.702 / 0.661）。公称の組成の点は比を出さない"""
    out = []
    for h, f in zip(half, full, strict=True):
        if not same_burnup(h["burnup"], f["burnup"]):
            raise SystemExit("the density paths have other burnup points")
        out.append(
            {
                "burnup": h["burnup"],
                **{
                    name: ratio(h[key], f[key])
                    for key, name in (
                        ("index_pcm", "index_ratio"),
                        ("production", "production_ratio"),
                        ("thermal_absorption", "absorption_ratio"),
                        ("linearization_pcm", "linearization_ratio"),
                    )
                },
            }
        )
    return out


def nominal_density(nominal, burnup, i):
    b, _ = split.nominal_state(nominal, burnup)
    return nominal["micro"]["densities"][b][i]


def compare(tables, summary, split_summary, names=NUCLIDES):
    """要約（モジュールの説明）。条件と組成の検査は compare_split.py の compare"""
    split.compare(tables, summary)
    nominal = tables["nominal"]
    at = positions(nominal, names)
    condition, one = split.split(nominal, tables["one"], tables["one_hybrid"], "R(1)")
    _, switch = split.split(
        nominal, tables["switch"], tables["switch_hybrid"], "switch"
    )
    _, half = split.split(
        tables["nominal_half"], tables["half"], tables["half_hybrid"], "w = 0.5"
    )
    _, full = split.split(
        tables["nominal_full"], tables["full"], tables["full_hybrid"], "w = 1"
    )
    weights = []
    for x in switch:
        found = [s for s in summary["rows"] if same_burnup(s["burnup"], x["burnup"])]
        if len(found) != 1:
            raise SystemExit(f"the switch summary has no single point at {x['burnup']}")
        weights.append(found[0]["weight"])
    rows = {"one": one, "switch": switch, "half": half, "full": full}
    histories = {
        name: history_rows(
            nominal,
            rows[name],
            weights if name == "switch" else [WEIGHTS[name]] * len(rows[name]),
            at,
        )
        for name in HISTORIES
    }
    to_one(histories)
    return {
        "format": "xlc-compare-history-linearization-index",
        "continued_burnup": summary["continued_burnup"],
        "nuclides": list(names),
        "histories": histories,
        "density_pairs": density_pairs(histories["half"], histories["full"]),
        "switch": switch_rows(
            nominal, condition, switch, one, summary, split_summary, at
        ),
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def span(values, form="+.2f"):
        values = [v for v in values if v is not None]
        return f"{min(values):{form}}〜{max(values):{form}}"

    u0 = summary["continued_burnup"]
    histories = summary["histories"]
    names = {"one": "R(1)", "switch": "switch", "half": "0.702", "full": "0.661"}
    later = [x["burnup"] for x in histories["switch"] if x["burnup"] > u0]
    last = later[-1]

    def at(name, u):
        return next(x for x in histories[name] if same_burnup(x["burnup"], u))

    lines = []
    for key, label, form in (
        ("linearization_first", "L/E", ".4f"),
        ("linearization_second", "L/E² (10⁻⁶/pcm)", ".3f"),
        ("untracked_first", "U/E", ".4f"),
        ("untracked_weight", "U/w (pcm)", "+.2f"),
    ):
        scale = 1e6 if key == "linearization_second" else 1.0
        values = "、".join(
            f"{names[n]} {scale * at(n, last)[key]:{form}}" for n in HISTORIES
        )
        spans = "、".join(
            f"{names[n]} {span([at(n, u)[f'{key}_to_one'] for u in later], '.4f')}"
            for n in HISTORIES[1:]
        )
        lines.append(
            f"{label}: {last:.4g} MWd/kg で {values}。R(1) の同じ燃焼度の値との比は"
            f" {later[0]:.4g}〜{last:.4g} MWd/kg の {len(later)} 点で {spans}"
        )
    pairs = [(at("half", u), at("full", u)) for u in later]
    lines.append(
        f"密度だけの経路（それぞれの履歴の定数）: L の効きの比 L(0.702)/L(0.661) "
        f"{span([ratio(h['linearization_pcm'], f['linearization_pcm']) for h, f in pairs], '.4f')}、"
        f"E の比 {span([ratio(h['index_pcm'], f['index_pcm']) for h, f in pairs], '.4f')}、U の比 "
        f"{span([ratio(h['untracked_pcm'], f['untracked_pcm']) for h, f in pairs], '.4f')}"
        f"（{later[0]:.4g}〜{last:.4g} MWd/kg の {len(later)} 点）。最後の点で L の比 "
        f"{pairs[-1][0]['linearization_pcm'] / pairs[-1][1]['linearization_pcm']:.4f}、E の比 "
        f"{pairs[-1][0]['index_pcm'] / pairs[-1][1]['index_pcm']:.4f}、U の比 "
        f"{pairs[-1][0]['untracked_pcm'] / pairs[-1][1]['untracked_pcm']:.4f}"
    )
    rows = [x for x in summary["switch"] if x["burnup"] > u0]
    s = rows[-1]
    lines.append(
        f"switch の {len(rows)} 点: 指標の比 E_switch/E(1) {span([x['index_ratio'] for x in rows], '.4f')}"
        f"（条件 j の指標の比 {span([x['condition_index_ratio'] for x in rows], '.4f')}、w "
        f"{span([x['weight'] for x in rows], '.4f')}、w_L {span([x['linearization_weight'] for x in rows], '.4f')}、"
        f"w_U {span([x['untracked_weight'] for x in rows], '.4f')}）。ΔN̄ の比 "
        + "、".join(
            f"{n} {span([x['density_ratios'][n] for x in rows], '.4f')}"
            for n in summary["nuclides"]
        )
    )
    lines.append(
        f"L の外れの分: w で {span([x['linearization_gap_weight_pcm'] for x in rows], '+.3f')}、指標の比で "
        f"{span([x['linearization_gap_index_pcm'] for x in rows], '+.3f')}、指標の比の 2 乗で "
        f"{span([x['linearization_gap_index_squared_pcm'] for x in rows], '+.3f')}、条件 j の指標の比で "
        f"{span([x['linearization_gap_condition_index_pcm'] for x in rows], '+.3f')} pcm。最後の点 {s['burnup']:.4g} MWd/kg で "
        f"{s['linearization_gap_weight_pcm']:+.3f}、{s['linearization_gap_index_pcm']:+.3f}、"
        f"{s['linearization_gap_index_squared_pcm']:+.3f}、{s['linearization_gap_condition_index_pcm']:+.3f} pcm"
    )
    lines.append(
        f"U の外れの分: w で {span([x['untracked_gap_weight_pcm'] for x in rows], '+.3f')}、指標の比で "
        f"{span([x['untracked_gap_index_pcm'] for x in rows], '+.3f')}、指標の比の 2 乗で "
        f"{span([x['untracked_gap_index_squared_pcm'] for x in rows], '+.3f')} pcm。両方の外れの分は、U を w、L を指標の比で "
        f"{span([x['gap_mixed_pcm'] for x in rows], '+.3f')} pcm（どちらも w で {span([x['gap_weight_pcm'] for x in rows], '+.3f')} pcm、"
        f"core − 参照 {span([x['difference_pcm'] for x in rows], '+.3f')} pcm）"
    )
    lines.append(
        f"熱群の νΣf の補正の比と熱群の吸収の補正の比（switch の {len(rows)} 点）: "
        f"{span([x['production_ratio'] for x in rows], '.4f')}、{span([x['absorption_ratio'] for x in rows], '.4f')}。"
        f"L の外れの分は {span([x['linearization_gap_production_pcm'] for x in rows], '+.3f')}、"
        f"{span([x['linearization_gap_absorption_pcm'] for x in rows], '+.3f')} pcm。U を w、L を熱群の νΣf の補正の比で足した"
        f"両方の外れの分は {span([x['gap_production_pcm'] for x in rows], '+.3f')} pcm"
    )
    pairs = [x for x in summary["density_pairs"] if x["index_ratio"] is not None]
    lines.append(
        f"密度だけの経路の比（0.702 / 0.661、燃焼度 0 を除く {len(pairs)} 点）: E {span([x['index_ratio'] for x in pairs], '.4f')}、"
        f"熱群の νΣf の補正 {span([x['production_ratio'] for x in pairs], '.4f')}、熱群の吸収の補正 "
        f"{span([x['absorption_ratio'] for x in pairs], '.4f')}、L {span([x['linearization_ratio'] for x in pairs], '.4f')}。"
        f"R(1) の E の符号は {', '.join(f'{x['burnup']:.4g}: {x['index_pcm']:+.2f}' for x in histories['one'] if 0 < x['burnup'] <= 6)} pcm、"
        f"熱群の νΣf の補正は {span([x['production'] for x in histories['one'] if x['burnup'] > 0], '.4e')} 1/cm"
    )
    # 記録の表（verification_tables.py の signed）と同じ負の符号。核種の名前（U-235）のハイフンは残す
    return re.sub(r"(?<![A-Za-z])-", "−", "\n".join(lines))


def check_summaries(summaries, md5s):
    """`switch.json` と `split.json` が今のテーブルから作られたこと"""
    split.check_summary(summaries["switch.json"], md5s)
    made = summaries["split.json"].get("inputs", {})
    for file in (*split.FILES.values(), "switch.json"):
        if made.get(file) != md5s[file]:
            raise SystemExit(f"split.json was not made from this {file}")


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
    md5s = {file: md5(here / file) for file in INPUTS}
    summaries = {}
    for file in SUMMARIES:
        with open(here / file, encoding="utf-8") as f:
            summaries[file] = json.load(f)
    check_summaries(summaries, md5s)
    tables = {k: read(here / v) for k, v in split.FILES.items()}
    result = compare(tables, summaries["switch.json"], summaries["split.json"])
    result["inputs"] = md5s
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
