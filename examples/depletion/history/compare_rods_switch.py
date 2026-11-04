"""制御棒を途中まで入れて抜く燃焼（ADR-0215 の確認の制御棒の履歴、ADR-0228、#0190）: 格子計算の参照の k∞ と、履歴の補正の core と
補正を切った core の k を比べる。

    python3 compare_rods_switch.py --json rods-switch.json     # このディレクトリの出力から JSON を書く
    python3 compare_rods_switch.py --report rods-switch.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、減速材の密度と制御棒の履歴の残りを持つ公称のテーブル（`lattice histories residual-rods.toml` の出力）、格子計算の参照（2A-rods の
履歴の 20 MWd/kg の点から制御棒を抜いて続けた燃焼の履歴。`switch-rods-history.toml`）、core の 3 つの燃焼の履歴（どれも集合体 1 体、
ノード 1 つ）:

- 履歴の補正の core（`core-2A-rods-switch.toml`。20 MWd/kg まで制御棒を入れ、その後を抜く。H_r を燃焼度で重みを付けた平均で進める）
- 補正を切った core（`core-2A-rods-plain-switch.toml`。履歴の残りを持たない公称のテーブルで同じ燃焼）
- 一生を制御棒を入れて燃焼させた core（`core-2A-rods-branch.toml`。w = 1）

と、制御棒の履歴の分岐の自分の条件のテーブル（`residual-rods-branch.toml`。値 aic の状態が第 1 項）。

- 切り替え: 参照の燃焼点のうち切り替えの点（参照の履歴の `continued`）から後ごとに、core の k の参照の k∞ からの差（補正あり、
  補正なし）と、追う核種の数密度の相対差（主な核種）。参照の数密度は燃焼の領域の数密度を体積で重みを付けて格子の面積で割った
  均質化した数密度（lattice の `[history]` の N̄ᵢ と同じ。面積は燃焼度 0 の U-238 で公称のテーブルの N̄ᵢ⁰ に合わせ、U-235 で照らす）
- 一生: 履歴の分岐の自分の条件のテーブルの値 aic の燃焼度の点ごとに、一生を制御棒を入れた core の k の、その k∞ からの差
- 分け方（ADR-0228）: 2 つの変数の R の点ごとに、ノードの組成の差が見込み w_m ΔN̄ᵢ^m + w_r ΔN̄ᵢ^r（w_m = w_r = `SPLIT_WEIGHT`）の、公称の
  状態（制御棒の値 out）の定数に線形の補正と w_m U_m + w_r U_r を足した定数に、L を ADR-0228 の a_j（このノードでは a_j = w_j）で
  足した k と、変数ごとに P / P_j で重ねて足した k（ADR-0228 の案 A）の差。U と L は制御棒の値 out の補正

ρ = 1 − 1/k で、差は pcm（a − b は ρ(a) − ρ(b)）。
"""

import argparse
import json
from pathlib import Path

from compare_core_history import same_burnup
from compare_residual import added, coefficients, corrected, k_inf, md5, read, rho

FILES = {
    "table": "table-residual-rods-histories.msgpack",
    "branch": "table-residual-rods-branch.msgpack",
    "reference": "lattice-switch-rods.msgpack",
    "core": "core-2A-rods-switch.msgpack",
    "plain": "core-2A-rods-plain-switch.msgpack",
    "life": "core-2A-rods-branch.msgpack",
    # 補正を切った core のテーブル（読まず、md5 を残す）
    "plain_table": "table-residual-rods-nominal.msgpack",
}
VARIABLES = ["moderator-density", "rods"]
INSERTED = "aic"
# 数密度の相対差を要約に残す主な核種
MAJOR = ["U-235", "Pu-239", "Pu-240", "Pu-241", "Xe-135"]
# 分け方の比べ方の、2 つの変数の重み
SPLIT_WEIGHT = 0.5


def nodes(core, label, nuclides):
    """集合体 1 体・ノード 1 つの core の履歴の、燃焼点ごとの（燃焼度、状態の名前、k、数密度、履歴の変数）"""
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"{label}: not an xlc-core-history")
    if len(core["assemblies"]) != 1 or len(core["axial_mesh"]) != 2:
        raise SystemExit(f"{label}: the core must be one assembly of one axial node")
    if core.get("nuclides") != nuclides:
        raise SystemExit(f"{label}: the core tracks other nuclides than the table")
    if core.get("history_variables") != VARIABLES:
        raise SystemExit(
            f"{label}: the core holds the history variables {core.get('history_variables')}"
        )
    return [
        {
            "burnup": p["burnups"][0][0],
            "state": p["state"],
            "k": p["k"],
            "densities": p["densities"][0][0],
            "history": p["histories"][0][0],
        }
        for p in core["points"]
    ]


def at_burnup(rows, burnup, label):
    found = [r for r in rows if same_burnup(r["burnup"], burnup)]
    if len(found) != 1:
        raise SystemExit(f"{label}: no single point at burnup {burnup}")
    return found[0]


def homogenized(history, nuclides, nominal_first):
    """格子計算の燃焼の履歴の燃焼点ごとの（燃焼度、k∞、追う核種の均質化した数密度）。格子の面積は、燃焼度 0 の U-238 の Σ_r V_r N_r を
    公称のテーブルの N̄ᵢ⁰（`nominal_first`。追う核種の順）で割ったもので、U-235 が相対 1e-9 で合うこと（初めの燃料は同じ）"""
    index = {n: i for i, n in enumerate(history["nuclides"])}
    volumes = [r["volume"] for r in history["regions"]]

    def summed(point, name):
        return sum(v * d[index[name]] for v, d in zip(volumes, point["densities"]))

    first = history["points"][0]
    if first["burnup"] != 0.0:
        raise SystemExit("the reference history must start at burnup 0")
    area = summed(first, "U-238") / nominal_first[nuclides.index("U-238")]
    check = summed(first, "U-235") / area
    expected = nominal_first[nuclides.index("U-235")]
    if abs(check - expected) > 1e-9 * expected:
        raise SystemExit(
            f"the reference U-235 at burnup 0 is {check}, the table {expected}"
        )
    return [
        {
            "burnup": p["burnup"],
            "k_inf": p["k_inf"],
            "densities": [summed(p, n) / area for n in nuclides],
        }
        for p in history["points"]
    ]


def percent(a, b):
    return 100.0 * (a - b) / b


def state_index(table, burnup, rod):
    """燃焼度の点 `burnup` の公称の状態（連続の軸は公称、制御棒の値 `rod`）の状態の番号"""
    point = []
    for a in table["axes"]:
        if a["name"] == "burnup":
            point.append(a["points"].index(burnup))
        elif a["name"] == "rod":
            point.append(a["values"].index(rod))
        else:
            point.append(a["points"].index(a["nominal"]))
    return table["points"].index(point)


def entry(table, variable, rod):
    (history,) = [h for h in table["histories"] if h["variable"] == variable]
    (found,) = [r for r in history["rods"] if r["rod"] == rod]
    return found


def production(coeffs, deltas):
    """熱群（最後の群）の νΣf の線形の補正（ADR-0215 の P）"""
    return sum(c["nu_fission"][-1] * d for c, d in zip(coeffs, deltas))


def scaled(r, w):
    out = {}
    for key, value in r.items():
        out[key] = (
            [[w * x for x in row] for row in value]
            if isinstance(value[0], list)
            else [w * x for x in value]
        )
    return out


def split_rows(table):
    """分け方の比べ方（モジュールの説明）。2 つの変数がどちらも R を持つ燃焼度の点のうち、履歴の分岐の数密度が公称と違う点"""
    micro = table["micro"]
    density, rods = entry(table, VARIABLES[0], "out"), entry(table, VARIABLES[1], "out")
    burnups = table["axes"][0]["points"]
    rows = []
    for b in density["burnup"]:
        if b not in rods["burnup"]:
            continue
        u = burnups.index(b)
        km, kr = density["burnup"].index(b), rods["burnup"].index(b)
        nominal = micro["densities"][u]
        dm = [x - y for x, y in zip(density["densities"][km], nominal)]
        dr = [x - y for x, y in zip(rods["densities"][kr], nominal)]
        if all(x == 0.0 for x in dm) or all(x == 0.0 for x in dr):
            continue
        p = state_index(table, b, "out")
        coeffs = coefficients(micro["states"][p])
        w = SPLIT_WEIGHT
        node = [w * x + w * y for x, y in zip(dm, dr)]
        base = corrected(table["states"][p], coeffs, node)
        base = added(base, scaled(density["untracked"][km], w))
        base = added(base, scaled(rods["untracked"][kr], w))
        index, pm, pr = (
            production(coeffs, node),
            production(coeffs, dm),
            production(coeffs, dr),
        )
        expected = w * pm + w * pr
        shares = [index / pm * (w * pm / expected), index / pr * (w * pr / expected)]
        independent = [index / pm, index / pr]
        lm, lr = density["linearization"][km], rods["linearization"][kr]

        def with_l(a, base=base, lm=lm, lr=lr):
            return added(added(base, scaled(lm, a[0])), scaled(lr, a[1]))

        rows.append(
            {
                "burnup": b,
                "shares": shares,
                "independent": independent,
                "density_l_pcm": rho(k_inf(added(base, scaled(lm, w))))
                - rho(k_inf(base)),
                "rods_l_pcm": rho(k_inf(added(base, scaled(lr, w)))) - rho(k_inf(base)),
                "overcount_pcm": rho(k_inf(with_l(independent)))
                - rho(k_inf(with_l(shares))),
            }
        )
    return rows


def compare(table, branch, reference, core, plain, life):
    nuclides = table["micro"]["nuclides"]
    names = [h["variable"] for h in table["histories"]]
    if names != VARIABLES:
        raise SystemExit(f"the table holds the history variables {names}")
    switch = reference.get("continued")
    if switch is None:
        raise SystemExit("the reference is not a continued history")
    ref = homogenized(reference, nuclides, table["micro"]["densities"][0])
    core_nodes = nodes(core, "core", nuclides)
    plain_nodes = nodes(plain, "plain", nuclides)
    major = [nuclides.index(n) for n in MAJOR]
    switch_rows = []
    for r in ref:
        if r["burnup"] < switch:
            continue
        c = at_burnup(core_nodes, r["burnup"], "core")
        q = at_burnup(plain_nodes, r["burnup"], "plain")
        switch_rows.append(
            {
                "burnup": r["burnup"],
                "weight": c["history"][1],
                "density_weight": c["history"][0],
                "reference_k": r["k_inf"],
                "core_pcm": rho(c["k"]) - rho(r["k_inf"]),
                "plain_pcm": rho(q["k"]) - rho(r["k_inf"]),
                "core_densities": {
                    nuclides[i]: percent(c["densities"][i], r["densities"][i])
                    for i in major
                },
                "plain_densities": {
                    nuclides[i]: percent(q["densities"][i], r["densities"][i])
                    for i in major
                },
            }
        )
    life_nodes = nodes(life, "life", nuclides)
    burnups = branch["axes"][0]["points"]
    life_rows = []
    for b in burnups:
        k = branch["states"][state_index(branch, b, INSERTED)]["k_inf"]
        c = at_burnup(life_nodes, b, "life")
        life_rows.append(
            {"burnup": b, "branch_k_inf": k, "life_pcm": rho(c["k"]) - rho(k)}
        )
    return {
        "format": "xlc-compare-rods-switch",
        "switch": switch,
        "switch_rows": switch_rows,
        "life_rows": life_rows,
        "split_weight": SPLIT_WEIGHT,
        "split_rows": split_rows(table),
    }


def report(summary):
    def span(rows, key):
        values = [r[key] for r in rows]
        return f"{min(values):+.3f}〜{max(values):+.3f}"

    after = [r for r in summary["switch_rows"] if r["burnup"] > summary["switch"]]
    lines = [
        (
            f"切り替えの後（{summary['switch']} MWd/kg より後の {len(after)} 点）の core − 参照: "
            f"補正あり {span(after, 'core_pcm')} pcm、補正なし {span(after, 'plain_pcm')} pcm"
        ),
        "絶対値の最大: 補正あり {:.3f} pcm、補正なし {:.3f} pcm".format(
            max(abs(r["core_pcm"]) for r in after),
            max(abs(r["plain_pcm"]) for r in after),
        ),
        f"一生を制御棒を入れた core − 第 1 項: {span(summary['life_rows'], 'life_pcm')} pcm",
        (
            f"2 つの変数の L を P / P_j で重ねた形 − ADR-0228 の形（w = {summary['split_weight']}）: "
            f"{span(summary['split_rows'], 'overcount_pcm')} pcm"
        ),
    ]
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
    read_files = {k: read(here / v) for k, v in FILES.items() if k != "plain_table"}
    result = compare(**read_files)
    result["inputs"] = {v: md5(here / v) for v in FILES.values()}
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
