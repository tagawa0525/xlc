"""履歴の残り R（ADR-0215）の大きさと、瞬時の状態への依存を、燃焼度の点ごとに出す。

    python3 compare_residual.py --json residual.json     # このディレクトリの出力から JSON を書く
    python3 compare_residual.py --report residual.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、R を書き足した公称のテーブル（`lattice histories residual.toml` の出力。減速材の密度の履歴の補正を 1 つ、制御棒の値 out で
持つ）、履歴の分岐の自分の条件のテーブル（R の第 1 項）、履歴の分岐の組成を公称の条件で解いた逆向きのテーブル、履歴の分岐の燃焼の
履歴。R を持つ燃焼度の点 u ごとに、条件 j（減速材の軸だけが Hᵇ の格子点の状態）で:

- 効き: ρ(分岐の定数) − ρ(公称の条件 j の状態の定数)。履歴の効き（組成の違いの、条件 j での効き）
- 線形: ρ(公称の条件 j の状態 + Σᵢ cᵢ ΔN̄ᵢ) − ρ(公称の条件 j の状態)。core の線形の補正の効き（係数 cᵢ は数密度の微分を持つ核種は
  その微分、ほかは σᵢ。公称のテーブルの条件 j の状態の値。ΔN̄ᵢ = N̄ᵢ^j − N̄ᵢ^0）
- R: ρ(線形の補正の後 + R) − ρ(線形の補正の後)。効き = 線形 + R
- 再現: ρ(線形の補正の後 + R) − ρ(分岐のテーブルの k∞)（R の定義の検査。2 群の k とテーブルの k∞ の差を含む）と、分岐の定数の 2 群の
  k との差（丸めだけ）
- 逆向きの R: 公称の条件の状態で同じ式で作った残り R⁰ = f(分岐の組成、公称の条件) − (公称の状態 + Σᵢ cᵢ⁰ ΔN̄ᵢ) の効き
  ρ(逆向きの定数) − ρ(公称の状態 + Σᵢ cᵢ⁰ ΔN̄ᵢ)。依存 = R の効き − 逆向きの R の効き
- 公称の状態での R: ρ(公称の状態 + Σᵢ cᵢ⁰ ΔN̄ᵢ + R) − ρ(公称の状態 + Σᵢ cᵢ⁰ ΔN̄ᵢ)。条件 j の R の、公称の状態での効き
- 公称の状態で使った誤差: ρ(公称の状態 + Σᵢ cᵢ⁰ ΔN̄ᵢ + R) − ρ(逆向きの定数) = 公称の状態での R − 逆向きの R。core が条件 j の R を
  公称の状態のノードに足したときの誤差（R は瞬時の状態に依らないと置く。ADR-0215）
- 差の全体を公称の状態で使った誤差: ρ(公称の状態 + (分岐の定数 − 公称の条件 j の状態)) − ρ(逆向きの定数)。線形の補正を引かない
  定数の差の全体を、線形の補正なしに公称の状態のノードに足したときの誤差（ADR-0215 の案 S1 の R'。ノードの数密度が履歴の分岐の値
  なら案 S1 の線形の補正は 0）
- 第 1 項: ρ(分岐のテーブルの k∞) − ρ(履歴の分岐の燃焼の点の k∞)（ADR-0187 の確認と同じ比べ方）

ρ は反応度 1 − 1/k で、差は pcm。2 群の k は少数群定数（全断面積、散乱行列、νΣf、χ）の無限媒質の固有値で、テーブルの k∞ との差
（`self_check_pcm`。3 つのテーブルのすべての状態の最大）を出す。
"""

import argparse
import hashlib
import json
from pathlib import Path

import msgpack

FILES = {
    "nominal": "table-residual-nominal-histories.msgpack",
    "branch": "table-residual-density-branch.msgpack",
    "reverse": "table-residual-density-reverse.msgpack",
    "history": "lattice-residual-density.msgpack",
}
VARIABLE = "moderator-density"
STATE_AXES = ("boron", "moderator", "fuel")


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def solve(matrix, vector):
    """小さな連立 1 次方程式（部分ピボットの Gauss の消去）"""
    n = len(vector)
    a = [list(row) + [b] for row, b in zip(matrix, vector)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(a[r][c]))
        a[c], a[p] = a[p], a[c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            for j in range(c, n + 1):
                a[r][j] -= f * a[c][j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - sum(a[r][j] * x[j] for j in range(r + 1, n))) / a[r][r]
    return x


def k_inf(state):
    """少数群定数の無限媒質の k: (diag(Σt) − Sᵀ) φ = χ を解き、k = νΣf · φ（χ に比例する核分裂源は 1 つの固有値だけを持つ）"""
    total, scatter = state["total"], state["scatter"]
    groups = len(total)
    removal = [
        [(total[g] if g == h else 0.0) - scatter[h][g] for h in range(groups)]
        for g in range(groups)
    ]
    flux = solve(removal, state["chi"])
    return sum(n * f for n, f in zip(state["nu_fission"], flux))


def rho(k):
    """反応度 1 − 1/k (pcm)"""
    return 1e5 * (1.0 - 1.0 / k)


def coefficients(nuclides):
    """線形の補正の係数: 数密度の微分（ADR-0221）を持つ核種はそれ、ほかは微視的断面積（core の補正と同じ）"""
    return [n.get("derivative") or n for n in nuclides]


def corrected(state, coefficients, deltas):
    """状態の少数群定数に、核種ごとの係数 × 数密度の差を足したもの（全断面積、νΣf、核分裂、散乱行列。吸収は全断面積の係数から散乱の
    行の和を引いた係数。χ と拡散係数はそのまま）"""
    out = dict(state)
    for key in ("total", "nu_fission", "fission"):
        out[key] = [
            value + sum(c[key][g] * d for c, d in zip(coefficients, deltas))
            for g, value in enumerate(state[key])
        ]
    out["absorption"] = [
        value
        + sum(
            (c["total"][g] - sum(c["scatter"][g])) * d
            for c, d in zip(coefficients, deltas)
        )
        for g, value in enumerate(state["absorption"])
    ]
    out["scatter"] = [
        [
            value + sum(c["scatter"][g][h] * d for c, d in zip(coefficients, deltas))
            for h, value in enumerate(row)
        ]
        for g, row in enumerate(state["scatter"])
    ]
    return out


def difference(a, b):
    """少数群定数の差 a − b（R の項目だけ）"""
    out = {
        key: [x - y for x, y in zip(a[key], b[key])]
        for key in ("total", "absorption", "nu_fission", "fission", "chi", "diffusion")
    }
    for key in ("scatter", "adf"):
        out[key] = [[x - y for x, y in zip(p, q)] for p, q in zip(a[key], b[key])]
    return out


def added(state, residual):
    """状態の少数群定数に残り R を足したもの（R の項目だけ）"""
    out = dict(state)
    for key in ("total", "absorption", "nu_fission", "fission", "chi", "diffusion"):
        out[key] = [a + b for a, b in zip(state[key], residual[key])]
    for key in ("scatter", "adf"):
        out[key] = [
            [a + b for a, b in zip(p, q)] for p, q in zip(state[key], residual[key])
        ]
    return out


def residual(rod, k):
    """制御棒の値 `rod` の R の点 `k` の残り R = U + L（テーブルは追わない核種の分 U と線形化の誤差 L に分けて持つ。ADR-0215）"""
    return added(rod["untracked"][k], rod["linearization"][k])


def axis(table, name):
    found = [a for a in table["axes"] if a["name"] == name]
    if len(found) != 1:
        raise SystemExit(f"the table has no {name} axis")
    return found[0]


def single_state(table, label):
    """燃焼度の軸のほかは 1 点のテーブルの、ホウ素、減速材（温度と密度）、燃料の値と燃焼度の点"""
    names = [a["name"] for a in table["axes"]]
    if names != ["burnup", *STATE_AXES] or any(
        len(a["points"]) != 1 for a in table["axes"][1:]
    ):
        raise SystemExit(f"{label}: the table must vary only the burnup axis")
    burnups = table["axes"][0]["points"]
    if len(table["states"]) != len(burnups) or table.get("micro") is None:
        raise SystemExit(f"{label}: one state per burnup and micro data are needed")
    moderator = axis(table, "moderator")
    condition = (
        axis(table, "boron")["points"][0],
        (moderator["points"][0], moderator["density"][0]),
        axis(table, "fuel")["points"][0],
    )
    return condition, burnups


def compare(nominal, branch, reverse, history):
    histories = nominal.get("histories") or []
    if len(histories) != 1 or histories[0]["variable"] != VARIABLE:
        raise SystemExit(f"the nominal table must hold one {VARIABLE} history")
    (h,) = histories
    rods = {r["rod"]: r for r in h["rods"]}
    if list(rods) != ["out"]:
        raise SystemExit("the history must hold the rod value out only")
    rod = rods["out"]
    names = [a["name"] for a in nominal["axes"]]
    if names != ["burnup", *STATE_AXES]:
        raise SystemExit(f"the nominal table has the axes {names}")
    burnups = nominal["axes"][0]["points"]
    moderator = axis(nominal, "moderator")
    nominal_index = [
        a["points"].index(a["nominal"]) for a in nominal["axes"][1:]
    ]  # ホウ素、減速材、燃料
    density = moderator["density"]
    condition_index = list(nominal_index)
    condition_index[1] = density.index(h["branch"])
    if density[nominal_index[1]] != h["nominal"]:
        raise SystemExit(
            "the history nominal value is not the nominal moderator density"
        )

    def values(index):
        return (
            axis(nominal, "boron")["points"][index[0]],
            (moderator["points"][index[1]], density[index[1]]),
            axis(nominal, "fuel")["points"][index[2]],
        )

    branch_condition, branch_burnups = single_state(branch, "branch")
    reverse_condition, reverse_burnups = single_state(reverse, "reverse")
    if branch_condition != values(condition_index):
        raise SystemExit("the branch table is not at the condition of the history")
    if reverse_condition != values(nominal_index):
        raise SystemExit("the reverse table is not at the nominal state")
    if branch_burnups != rod["burnup"] or reverse_burnups != rod["burnup"]:
        raise SystemExit("the branch tables do not solve the burnups of the residuals")
    micro = nominal["micro"]
    for table, label in ((branch, "branch"), (reverse, "reverse")):
        if table["micro"]["nuclides"] != micro["nuclides"]:
            raise SystemExit(f"{label}: other micro nuclides")
    # 逆向きは同じ履歴の組成を解いたもの
    if reverse["micro"]["densities"] != branch["micro"]["densities"]:
        raise SystemExit(
            "the reverse table holds other densities than the branch table"
        )
    history_k = {p["burnup"]: p["k_inf"] for p in history["points"]}
    points = {tuple(p): i for i, p in enumerate(nominal["points"])}
    self_check = max(
        abs(rho(k_inf(s)) - rho(s["k_inf"]))
        for t in (nominal, branch, reverse)
        for s in t["states"]
    )
    rows = []
    for k, u in enumerate(rod["burnup"]):
        b = burnups.index(u)
        pj = points[(b, *condition_index)]
        p0 = points[(b, *nominal_index)]
        deltas = [
            x - y
            for x, y in zip(branch["micro"]["densities"][k], micro["densities"][b])
        ]
        r = residual(rod, k)
        fj, f0 = nominal["states"][pj], nominal["states"][p0]
        fb, fr = branch["states"][k], reverse["states"][k]
        linear_j = corrected(fj, coefficients(micro["states"][pj]), deltas)
        linear_0 = corrected(f0, coefficients(micro["states"][p0]), deltas)
        k_linear_j = k_inf(linear_j)
        k_with_r = k_inf(added(linear_j, r))
        k_linear_0 = k_inf(linear_0)
        k_with_r_0 = k_inf(added(linear_0, r))
        residual_pcm = rho(k_with_r) - rho(k_linear_j)
        reverse_pcm = rho(k_inf(fr)) - rho(k_linear_0)
        if u not in history_k:
            raise SystemExit(f"the history has no burnup {u}")
        rows.append(
            {
                "burnup": u,
                "effect_pcm": rho(k_inf(fb)) - rho(k_inf(fj)),
                "linear_pcm": rho(k_linear_j) - rho(k_inf(fj)),
                "residual_pcm": residual_pcm,
                "reproduce_pcm": rho(k_with_r) - rho(fb["k_inf"]),
                "reproduce_constants_pcm": rho(k_with_r) - rho(k_inf(fb)),
                "reverse_residual_pcm": reverse_pcm,
                "dependence_pcm": residual_pcm - reverse_pcm,
                "residual_nominal_pcm": rho(k_with_r_0) - rho(k_linear_0),
                "applied_pcm": rho(k_with_r_0) - rho(k_inf(fr)),
                "applied_whole_pcm": rho(k_inf(added(f0, difference(fb, fj))))
                - rho(k_inf(fr)),
                "first_term_pcm": rho(fb["k_inf"]) - rho(history_k[u]),
            }
        )
    return {
        "format": "xlc-compare-history-residual",
        "variable": VARIABLE,
        "nominal_value": h["nominal"],
        "branch_value": h["branch"],
        "nuclides": len(micro["nuclides"]),
        "derivative_nuclides": [
            n
            for n, s in zip(micro["nuclides"], micro["states"][0])
            if s.get("derivative")
        ],
        "self_check_pcm": self_check,
        "rows": rows,
    }


def report(summary):
    rows = summary["rows"]

    def span(key):
        values = [r[key] for r in rows]
        return f"{min(values):+.2f}〜{max(values):+.2f}"

    def largest(key):
        worst = max(rows, key=lambda r: abs(r[key]))
        return f"{abs(worst[key]):.2g} pcm（{worst['burnup']:g} MWd/kg）"

    last = rows[-1]
    return (
        f"減速材の密度の履歴（H⁰ = {summary['nominal_value']:g}、Hᵇ = {summary['branch_value']:g} g/cc）、追う核種 "
        f"{summary['nuclides']}（微分は {'、'.join(summary['derivative_nuclides']) or 'なし'}）、燃焼度の点 {len(rows)}。"
        f"R の効きは {span('residual_pcm')} pcm（最後の点 {last['burnup']:g} MWd/kg で {last['residual_pcm']:+.2f} pcm。"
        f"そのときの効き {last['effect_pcm']:+.1f} pcm、線形の補正 {last['linear_pcm']:+.1f} pcm）。"
        f"R を足した定数の 2 群の k は、分岐のテーブルの k∞ を最大 {largest('reproduce_pcm')}、分岐の定数の 2 群の k を最大 "
        f"{largest('reproduce_constants_pcm')} の差で再現する。"
        f"\n状態への依存（R の効き − 逆向きの R の効き）は {span('dependence_pcm')} pcm（絶対値の最大 {largest('dependence_pcm')}）。"
        f"条件の R を公称の状態で使った誤差は {span('applied_pcm')} pcm（絶対値の最大 {largest('applied_pcm')}。"
        f"最後の点の公称の状態での R の効きは {last['residual_nominal_pcm']:+.2f} pcm、逆向きの R は "
        f"{last['reverse_residual_pcm']:+.2f} pcm）。"
        f"線形の補正を引かない定数の差の全体を公称の状態で使った誤差は {span('applied_whole_pcm')} pcm（絶対値の最大 "
        f"{largest('applied_whole_pcm')}）。"
        f"\n第 1 項の分岐は履歴の分岐の燃焼の点の k∞ を最大 {largest('first_term_pcm')} の差で再現する。2 群の k とテーブルの k∞ の"
        f"差は最大 {summary['self_check_pcm']:.2e} pcm"
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
    tables = {name: read(here / file) for name, file in FILES.items()}
    summary = compare(
        tables["nominal"], tables["branch"], tables["reverse"], tables["history"]
    )
    summary["inputs"] = {file: md5(here / file) for file in FILES.values()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
