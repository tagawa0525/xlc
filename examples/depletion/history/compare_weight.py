"""履歴の重みと R の w への依存（ADR-0215 の確認、#0190）: 途中で条件を変えた燃焼の、履歴の変数の重み ω を変えた core の差と、
減速材の密度だけの経路の R(0.5) と R(1)/2 の効きを比べる。

    python3 compare_weight.py --json weight.json     # このディレクトリの出力から JSON を書く
    python3 compare_weight.py --report weight.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

1. ω の実測: H_m を PARCS の形 H(B₀ + ΔB) = (H(B₀) B₀ + ω x ΔB) / (B₀ + ω ΔB) で進めた core の途中で条件を変えた燃焼
   （`core-2A-residual-switch-omega*.toml`）の k と、参照（`switch-density-branch.toml`）の k∞ の差を、参照の点ごとに出す。
   `switch.json` の、補正の量が w に比例すると置いた見込み（`omegas`）と並べる
2. R の w への依存: 2A を 565 K のまま減速材の密度を 0.702 g/cc（w = 0.5）と 0.661 g/cc（w = 1）にして燃焼させた履歴の分岐の、
   自分の条件の 1 状態のテーブル（第 1 項）と、公称の履歴を同じ条件で解いた 1 状態のテーブル（第 2 項と係数）から、燃焼度の点 u
   ごとに compare_residual.py と同じ式で R(w) = f(分岐の組成、条件) − f(公称の組成、条件) − Σᵢ cᵢ(条件) ΔN̄ᵢ を作る
   （`lattice histories` を通さない。テーブルの減速材の軸は同じ温度の 2 点を置けないため）。core は R を公称の状態で w R(1) と
   足すので、R(0.5) と R(1)/2 を、同じ公称の状態の定数（`table-residual-nominal.msgpack` の 565 K、0.743 g/cc の格子点）に
   w = 0.5 の数密度の線形の補正 Σᵢ cᵢ⁰ ΔN̄ᵢ(0.5) を足した定数 f⁰ + Σᵢ cᵢ⁰ ΔN̄ᵢ(0.5) に足した効き Δρ で比べる:
   差 = ρ(f⁰ + Σᵢ cᵢ⁰ ΔN̄ᵢ(0.5) + R(0.5)) − ρ(f⁰ + Σᵢ cᵢ⁰ ΔN̄ᵢ(0.5) + R(1)/2)。反応率の微視的断面積の差 Δσ(w) = σ(分岐の組成、
   条件) − σ(公称の組成、条件) も、主な反応の Δσ(0.5) / Δσ(1) の比で比べる（比例なら 0.5）
3. 検算: 600 K、0.661 g/cc の経路（温度の履歴を含む）の R(1) を、公称のテーブルの条件 j の状態と第 1 項
   （`table-residual-density-branch.msgpack`）から同じ式で作り直し、`lattice histories` が書いた R（`histories`）との差を出す。
   その条件 j での効きを、密度だけの経路の R(1) の自分の条件での効きと並べる（差が温度の履歴の分）
4. 案 A の見込み（core を流していない。確かめていない）: R(0) = 0、R(0.5)、R(1) を通る w の 2 次の式は
   R(w) = w R(1) + 4 w (1 − w) (R(0.5) − R(1)/2)。core の w R(1) を置き換えた差の変化を、途中で条件を変えた燃焼の参照の点
   ごとに 4 (1 − w) (r − 1/2) E と見込む。r = 効き(R(0.5)) / (2 効き(R(1)/2))（密度だけの経路の、上の比べ方の効きの比。比例なら
   1/2）、E は ω = 1 の core がノードに足した w R(1) の効き（core-history.json の `switch_rows` の `residual_pcm`）、w はその w。
   密度だけの経路の非線形の割合を、温度の履歴を含む経路に当てはめる（ADR-0215 の、ノードの密度と温度の履歴が格子点の組と
   同じ関係で動くという置き方のもとで。確かめていない）。Δσ の非線形は入れていない。残りの差が、案 B（重みの定義）が埋める分の見込み

ρ = 1 − 1/k で、差は pcm。2 群の k は少数群定数の無限媒質の k（compare_residual.py の `k_inf`）。
"""

import argparse
import json
from pathlib import Path

from compare_core_history import Table, same_burnup, single_node
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
from compare_switch import at, reference_states

FILES = {
    "nominal": "table-residual-nominal.msgpack",
    "histories": "table-residual-nominal-histories.msgpack",
    "branch": "table-residual-density-branch.msgpack",
    "nominal_half": "table-weight-nominal-0702.msgpack",
    "branch_half": "table-weight-branch-0702.msgpack",
    "nominal_one": "table-weight-nominal-0661.msgpack",
    "branch_one": "table-weight-branch-0661.msgpack",
    "reference": "table-switch-density-branch.msgpack",
    "switch": "switch.json",
    "core_history": "core-history.json",
    # ω = 1 の core の途中で条件を変えた燃焼（読まず、switch.json と core-history.json が同じ燃焼から作られたことを md5 で照らす）
    "core": "core-2A-residual-switch.msgpack",
}
CORES = {
    1.25: "core-2A-residual-switch-omega125.msgpack",
    1.5: "core-2A-residual-switch-omega150.msgpack",
    2.0: "core-2A-residual-switch-omega200.msgpack",
    2.5: "core-2A-residual-switch-omega250.msgpack",
}
# Δσ の比を出す反応（核種、反応）。反応は fission か MT の番号（経路の和）
REACTIONS = [
    ("U-235", "fission"),
    ("Pu-239", "fission"),
    ("Pu-240", 102),
    ("Xe-135", 102),
]
STATE_AXES = ("boron", "moderator", "fuel")


def condition_table(table, label):
    """燃焼度の軸のほかは 1 点のテーブルの、条件（ホウ素、減速材の温度と密度、燃料）と、燃焼度の点ごとの状態の番号"""
    names = [a["name"] for a in table["axes"]]
    if names != ["burnup", *STATE_AXES] or any(
        len(a["points"]) != 1 for a in table["axes"][1:]
    ):
        raise SystemExit(f"{label}: the table must vary only the burnup axis")
    if table.get("micro") is None:
        raise SystemExit(f"{label}: the table has no micro data")
    moderator = table["axes"][2]
    condition = (
        table["axes"][1]["points"][0],
        (moderator["points"][0], moderator["density"][0]),
        table["axes"][3]["points"][0],
    )
    burnups = table["axes"][0]["points"]
    return condition, burnups, [p for _, _, p in reference_states(table)]


def scaled(residual, factor):
    """残り R の項目に `factor` を掛けたもの"""
    out = {
        key: [factor * x for x in residual[key]]
        for key in ("total", "absorption", "nu_fission", "fission", "chi", "diffusion")
    }
    for key in ("scatter", "adf"):
        out[key] = [[factor * x for x in row] for row in residual[key]]
    return out


def reaction(sigma, which):
    """微視的断面積の反応の群ごとの値（fission か、MT の経路の和）"""
    if which == "fission":
        return sigma["fission"]
    found = [r for r in sigma.get("reactions", []) if r["mt"] == which]
    if len(found) != 1:
        raise SystemExit(f"no single reaction MT {which}")
    return [sum(c[g] for c in found[0]["channels"]) for g in range(len(sigma["total"]))]


def residuals(nominal_condition, branch_condition, label):
    """条件の公称の組成と分岐の組成の 1 状態のテーブルから、燃焼度の点ごとの（R、条件の線形の補正の後の定数、ΔN̄ᵢ、分岐の微視的
    断面積、公称の微視的断面積）。R の式は compare_residual.py と同じ（係数は公称の組成の条件の状態の値。微分を持つ核種は微分）"""
    c0, burnups, states0 = condition_table(nominal_condition, f"{label} nominal")
    c1, burnups1, states1 = condition_table(branch_condition, f"{label} branch")
    if c0 != c1:
        raise SystemExit(
            f"{label}: the nominal and branch tables solve other conditions"
        )
    if burnups != burnups1:
        raise SystemExit(f"{label}: the nominal and branch tables solve other burnups")
    if nominal_condition["micro"]["nuclides"] != branch_condition["micro"]["nuclides"]:
        raise SystemExit(f"{label}: the nominal and branch tables track other nuclides")
    out = []
    for b, (p0, p1) in enumerate(zip(states0, states1)):
        deltas = [
            x - y
            for x, y in zip(
                branch_condition["micro"]["densities"][b],
                nominal_condition["micro"]["densities"][b],
            )
        ]
        linear = corrected(
            nominal_condition["states"][p0],
            coefficients(nominal_condition["micro"]["states"][p0]),
            deltas,
        )
        out.append(
            {
                "burnup": burnups[b],
                "residual": difference(branch_condition["states"][p1], linear),
                "linear": linear,
                "deltas": deltas,
                "branch_sigma": branch_condition["micro"]["states"][p1],
                "nominal_sigma": nominal_condition["micro"]["states"][p0],
                "nominal_densities": nominal_condition["micro"]["densities"][b],
            }
        )
    return c0, out


def effect(state, residual):
    """状態 `state` に残り `residual` を足した効き (pcm)"""
    return rho(k_inf(added(state, residual))) - rho(k_inf(state))


def sigma_ratios(half, one, nuclides):
    """主な反応の Δσ(0.5) / Δσ(1)（群ごと。Δσ(1) が 0 の群は None）"""
    out = {}
    for name, which in REACTIONS:
        i = nuclides.index(name)
        d_half = [
            a - b
            for a, b in zip(
                reaction(half["branch_sigma"][i], which),
                reaction(half["nominal_sigma"][i], which),
            )
        ]
        d_one = [
            a - b
            for a, b in zip(
                reaction(one["branch_sigma"][i], which),
                reaction(one["nominal_sigma"][i], which),
            )
        ]
        label = f"{name} {'fission' if which == 'fission' else f'MT {which}'}"
        out[label] = [None if y == 0.0 else x / y for x, y in zip(d_half, d_one)]
    return out


def largest_difference(a, b):
    """2 つの残りの項目の差の絶対値の最大"""
    worst = 0.0
    for key in ("total", "absorption", "nu_fission", "fission", "chi", "diffusion"):
        worst = max([worst, *(abs(x - y) for x, y in zip(a[key], b[key]))])
    for key in ("scatter", "adf"):
        worst = max(
            [worst, *(abs(x - y) for p, q in zip(a[key], b[key]) for x, y in zip(p, q))]
        )
    return worst


def compare_residuals(nominal, histories, branch, half_tables, one_tables):
    """2 と 3（モジュールの説明）。燃焼度の点ごとの行"""
    t = Table(histories)
    nuclides = nominal["micro"]["nuclides"]
    burnups = nominal["axes"][0]["points"]
    nominal_index = [a["points"].index(a["nominal"]) for a in nominal["axes"][1:]]
    moderator = nominal["axes"][2]
    condition_index = list(nominal_index)
    condition_index[1] = moderator["density"].index(t.history["branch"])
    points = {tuple(p): i for i, p in enumerate(nominal["points"])}
    nominal_condition = (
        nominal["axes"][1]["points"][nominal_index[0]],
        (moderator["points"][nominal_index[1]], moderator["density"][nominal_index[1]]),
        nominal["axes"][3]["points"][nominal_index[2]],
    )
    for table, label in ((half_tables[0], "w = 0.5"), (one_tables[0], "w = 1")):
        if table["micro"]["nuclides"] != nuclides:
            raise SystemExit(f"{label}: other micro nuclides than the nominal table")
    half_condition, half = residuals(*half_tables, "w = 0.5")
    one_condition, one = residuals(*one_tables, "w = 1")
    if [x["burnup"] for x in half] != [x["burnup"] for x in one]:
        raise SystemExit("the w = 0.5 and w = 1 tables solve other burnups")
    # 密度だけの経路: 公称と同じホウ素、減速材の温度、燃料
    for c, label in ((half_condition, "w = 0.5"), (one_condition, "w = 1")):
        if (c[0], c[1][0], c[2]) != (
            nominal_condition[0],
            nominal_condition[1][0],
            nominal_condition[2],
        ):
            raise SystemExit(
                f"{label}: the condition changes more than the moderator density"
            )
    rod = t.history["rods"][0]
    rows = []
    for h, o in zip(half, one):
        u = h["burnup"]
        b = burnups.index(u)
        p0 = points[(b, *nominal_index)]
        pj = points[(b, *condition_index)]
        f0 = nominal["states"][p0]
        c0 = coefficients(nominal["micro"]["states"][p0])
        n0 = nominal["micro"]["densities"][b]
        for x, label in ((h, "w = 0.5"), (o, "w = 1")):
            # 公称の組成を条件で解いたテーブルは、公称のテーブルと同じ履歴の組成（ΔN̄ᵢ の基準が同じ）
            if x["nominal_densities"] != n0:
                raise SystemExit(
                    f"{label}: the nominal composition differs from the nominal table at {u}"
                )
        base_half = corrected(f0, c0, h["deltas"])
        base_one = corrected(f0, c0, o["deltas"])
        half_pcm = effect(base_half, h["residual"])
        one_half_pcm = effect(base_half, scaled(o["residual"], 0.5))
        # 600 K の経路の R(1)（作り直し）
        k = rod["burnup"].index(u)
        deltas_600 = [x - y for x, y in zip(branch["micro"]["densities"][k], n0)]
        linear_600 = corrected(
            nominal["states"][pj],
            coefficients(nominal["micro"]["states"][pj]),
            deltas_600,
        )
        residual_600 = difference(branch["states"][k], linear_600)
        rows.append(
            {
                "burnup": u,
                "half_pcm": half_pcm,
                "one_half_pcm": one_half_pcm,
                "difference_pcm": half_pcm - one_half_pcm,
                "one_pcm": effect(base_one, o["residual"]),
                "half_condition_pcm": effect(h["linear"], h["residual"]),
                "one_condition_pcm": effect(o["linear"], o["residual"]),
                "temperature_path_pcm": effect(linear_600, residual_600),
                "temperature_path_rebuilt": largest_difference(
                    residual_600, rod["residuals"][k]
                ),
                "sigma_ratios": sigma_ratios(h, o, nuclides),
            }
        )
    return {
        "half_condition": half_condition,
        "one_condition": one_condition,
        "nominal_condition": nominal_condition,
        "rows": rows,
    }


def compare_omegas(histories, reference, switch, cores):
    """1（モジュールの説明）。ω ごとに、参照の点の（w、core − 参照、見込み）"""
    t = Table(histories)
    nuclides = histories["micro"]["nuclides"]
    states = reference_states(reference)
    estimates = {x["omega"]: x for x in switch["omegas"]}
    out = []
    for omega, core in cores.items():
        if omega not in estimates:
            raise SystemExit(f"the switch summary has no estimate of omega {omega}")
        nodes = single_node(core, f"omega {omega}", nuclides)
        rows = []
        for (u, _, p), estimate, row in zip(
            states, estimates[omega]["estimate_pcm"], switch["rows"]
        ):
            if not same_burnup(row["burnup"], u):
                raise SystemExit(
                    "the switch summary has other burnups than the reference"
                )
            n = at(nodes, u, f"omega {omega}")
            rows.append(
                {
                    "burnup": u,
                    "weight": t.weight(n["history"]) + 0.0,
                    "difference_pcm": rho(n["k"])
                    - rho(reference["states"][p]["k_inf"]),
                    "estimate_pcm": estimate,
                }
            )
        out.append({"omega": omega, "rows": rows})
    return out


def quadratic(switch, residual_rows, core_switch):
    """4（モジュールの説明）。途中で条件を変えた燃焼の参照の点（続ける点の後）ごとの見込み。`core_switch` は core-history.json の
    `switch_rows`（ω = 1 の core の燃焼点ごとの w と、ノードに足した w R(1) の効き `residual_pcm`）"""
    by_burnup = {x["burnup"]: x for x in residual_rows}
    out = []
    for row in switch["rows"][1:]:
        u, w = row["burnup"], row["weight"]
        x = by_burnup.get(u)
        if x is None:
            raise SystemExit(f"no residual at burnup {u}")
        found = [y for y in core_switch if same_burnup(y["burnup"], u)]
        if len(found) != 1 or found[0]["weight"] != w:
            raise SystemExit(
                f"the core history summary has no single point of w {w} at burnup {u}"
            )
        ratio = x["half_pcm"] / (2.0 * x["one_half_pcm"])
        # w R(1) を w R(1) + 4 w (1 − w) (R(0.5) − R(1)/2) に置き換える。R(0.5) − R(1)/2 を (2 r − 1) R(1)/2 と置き、効きは R に
        # 1 次とする: 変化 = 4 (1 − w) (r − 1/2) × (core がノードに足した w R(1) の効き)。R(1) の効きは足す定数で変わる（同じ R(1) でも
        # 60 MWd/kg で w = 0.5 と 1 の公称の状態で 14% 違う）ので、条件 j の格子点の効きでなく、core のノードの効きを使う
        change = 4.0 * (1.0 - w) * (ratio - 0.5) * found[0]["residual_pcm"]
        out.append(
            {
                "burnup": u,
                "weight": w,
                "needed_weight": row["needed_weight"],
                "difference_pcm": row["difference_pcm"],
                "ratio": ratio,
                "residual_pcm": found[0]["residual_pcm"],
                "change_pcm": change,
                "estimate_pcm": row["difference_pcm"] + change,
            }
        )
    return out


def check_switch(switch, md5s, core_history=None):
    """`switch.json` が、同じ R を持つテーブル、参照、core の途中で条件を変えた燃焼から作ったものであること。`core_history`
    （core-history.json）があれば、それも同じ core の燃焼から作ったものであること"""
    for key in ("histories", "reference", "core"):
        made_from = switch.get("inputs", {}).get(FILES[key])
        if made_from != md5s[FILES[key]]:
            raise SystemExit(
                f"the switch summary was made from {FILES[key]} of md5 {made_from}, not {md5s[FILES[key]]}"
            )
    if core_history is not None:
        made_from = core_history.get("inputs", {}).get(FILES["core"])
        if made_from != md5s[FILES["core"]]:
            raise SystemExit(
                f"the core history summary was made from {FILES['core']} of md5 {made_from}, not {md5s[FILES['core']]}"
            )


def compare(tables, switch, core_history, cores):
    """要約（モジュールの説明）"""
    residual = compare_residuals(
        tables["nominal"],
        tables["histories"],
        tables["branch"],
        (tables["nominal_half"], tables["branch_half"]),
        (tables["nominal_one"], tables["branch_one"]),
    )
    return {
        "format": "xlc-compare-history-weight",
        "continued_burnup": switch["continued_burnup"],
        "omegas": compare_omegas(
            tables["histories"], tables["reference"], switch, cores
        ),
        "omega_one": [
            {
                "burnup": x["burnup"],
                "weight": x["weight"],
                "difference_pcm": x["difference_pcm"],
            }
            for x in switch["rows"]
        ],
        "residual": residual,
        "quadratic": quadratic(switch, residual["rows"], core_history["switch_rows"]),
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""
    u0 = summary["continued_burnup"]

    def span(values, form="+.3f"):
        return f"{min(values):{form}}〜{max(values):{form}}"

    later = [x for x in summary["omega_one"] if x["burnup"] > u0]
    lines = [
        (
            f"ω = 1: {u0:.4g} MWd/kg より後の参照の点 {len(later)} で core − 参照は "
            f"{span([x['difference_pcm'] for x in later])} pcm"
        )
    ]
    for o in summary["omegas"]:
        rows = [x for x in o["rows"] if x["burnup"] > u0]
        gap = [x["difference_pcm"] - x["estimate_pcm"] for x in rows]
        last = rows[-1]
        lines.append(
            f"ω = {o['omega']:g}: 実測 {span([x['difference_pcm'] for x in rows])} pcm（見込み "
            f"{span([x['estimate_pcm'] for x in rows], '+.1f')} pcm、実測 − 見込み {span(gap)} pcm）。"
            f"最後の点 {last['burnup']:.4g} MWd/kg で w {last['weight']:.4f}、実測 {last['difference_pcm']:+.3f} pcm、"
            f"見込み {last['estimate_pcm']:+.1f} pcm"
        )
    r = summary["residual"]
    rows = r["rows"]
    last = rows[-1]
    worst = max(rows, key=lambda x: abs(x["difference_pcm"]))
    lines.append(
        f"R(0.5) − R(1)/2 の効き（公称の状態 + w = 0.5 の線形の補正に足した効きの差）は {span([x['difference_pcm'] for x in rows])} pcm"
        f"（絶対値の最大 {abs(worst['difference_pcm']):.3f} pcm、{worst['burnup']:.4g} MWd/kg）。最後の点 {last['burnup']:.4g} MWd/kg で "
        f"R(0.5) の効き {last['half_pcm']:+.3f}、R(1)/2 の効き {last['one_half_pcm']:+.3f}、差 {last['difference_pcm']:+.3f} pcm。"
        f"R(1) の効き（公称の状態）は {span([x['one_pcm'] for x in rows], '+.2f')} pcm"
    )
    lines.append(
        f"検算: 600 K の経路の R(1) の作り直しと `histories` の R の差は最大 "
        f"{max(x['temperature_path_rebuilt'] for x in rows):.2e}。条件での R(1) の効きは、600 K の経路で "
        f"{span([x['temperature_path_pcm'] for x in rows], '+.2f')} pcm、565 K の密度だけの経路で "
        f"{span([x['one_condition_pcm'] for x in rows], '+.2f')} pcm（最後の点で {last['temperature_path_pcm']:+.2f} と "
        f"{last['one_condition_pcm']:+.2f} pcm）"
    )
    ratios = {}
    for x in rows[1:]:
        for name, values in x["sigma_ratios"].items():
            ratios.setdefault(name, []).extend(v for v in values if v is not None)
    lines.append(
        "Δσ(0.5) / Δσ(1)（燃焼度 0 を除く点と群）: "
        + "、".join(f"{n} {span(v, '.3f')}" for n, v in ratios.items())
        + "；最後の点: "
        + "、".join(
            f"{n} " + "/".join("–" if v is None else f"{v:.3f}" for v in values)
            for n, values in last["sigma_ratios"].items()
        )
    )
    q = summary["quadratic"]
    lines.append(
        f"案 A（2 次の内挿）の見込み: 差の変化 {span([x['change_pcm'] for x in q], '+.2f')} pcm、差 "
        f"{span([x['estimate_pcm'] for x in q], '+.2f')} pcm（ω = 1 は {span([x['difference_pcm'] for x in q], '+.3f')} pcm）。"
        f"効きの比 r = {span([x['ratio'] for x in q], '.4f')}"
    )
    return "\n".join(lines)


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
    md5s = {file: md5(here / file) for file in [*FILES.values(), *CORES.values()]}
    with open(here / FILES["switch"], encoding="utf-8") as f:
        switch = json.load(f)
    with open(here / FILES["core_history"], encoding="utf-8") as f:
        core_history = json.load(f)
    check_switch(switch, md5s, core_history)
    tables = {
        k: read(here / v)
        for k, v in FILES.items()
        if k not in ("switch", "core_history", "core")
    }
    cores = {omega: read(here / file) for omega, file in CORES.items()}
    summary = compare(tables, switch, core_history, cores)
    summary["inputs"] = md5s
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
