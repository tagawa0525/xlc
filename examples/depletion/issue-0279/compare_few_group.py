"""少数群定数の数密度の微分（lattice の `[sensitivity.few-group]`。ADR-0247 の段 4）を、同じ入力の直接の摂動の中心差分と
核種ごと、定数ごとに比べ、自己遮蔽の暗黙の効きの大きさを書く（#0279）。

- 自己遮蔽を固定した微分（案 I0）: 一般化随伴の値 G（`relative`）と、数密度を (1 ± ε) 倍と (1 ± ε/2) 倍した直接の摂動の
  Richardson の外挿 D の差を、許容 max(2 e, 1e-4 |D|, 3 τ R / ε) と比べる。e = |D(ε/2) − D(ε)| / 3 は刻みの誤差の見積もり、
  3 τ R / ε は中性子束の判定 τ の相対の誤差が R の外挿に入る大きさ（compare_k_sensitivity.py の k の判定の項と同じ形）
- 自己遮蔽を解き直した核種（`shielded`）: 暗黙の効きは、自己遮蔽を解き直した直接の摂動 D_s と D の差（案 I3 の量）。
  自己遮蔽した物質の差を一般化随伴に掛けた値 I1（`shielded.first_order`）を D_s と、両方の刻みの誤差を含む同じ形の許容で比べる。
  k の I1（随伴に掛けた値）も、自己遮蔽を解き直した k の直接の摂動と比べる
- 反応度: 格子全体の 2 群の定数の無限媒質の k∞ = νΣfᵀ A⁻¹ χ（A = diag Σt − Σsᵀ。χ はテーブルの値で固定）の、定数による
  微分で、定数の変化を k の変化に直し、Δρ = Δk∞ / k² を pcm で書く。鏡面反射の格子では k∞ は MOC の k と一致する（反応率と
  中性子束を保存する均質化）。自己遮蔽を解き直した k の直接の摂動の差（Δk / k²）も書く。どちらも数密度を 1 割変えたときの値
  （ADR-0221 の摂動の割合 0.1。1 次の範囲で割合に比例する）
- 隣り合う燃焼度の点の組の自己遮蔽の分の見込み（`--pair`）: 核種ごとの暗黙の効き（k の直接の摂動の差、数密度の相対の変化
  あたり）に、組の数密度の相対の変化（inputs.py の pair-56-60.toml）を掛けて足す。ADR-0221 の背景の「自己遮蔽の分」（線形の
  補正が表さない、組成の差による自己遮蔽の変化の効き）の 1 次の見込み。組の変化はリングごとに違うが、均質化した数密度の
  相対の変化で、すべてのリングを同じ割合で変える向き（微分の向き）に置き換える

使い方: python3 compare_few_group.py [--pair <組の TOML>] <結果の TOML>...（I0 か I1（少数群定数と k）のどれかが許容を超えれば終了コード 1）
"""

import sys
from pathlib import Path

import numpy as np
import tomllib

RELATIVE = 1.0e-4
# 2 群のテーブルの定数の名前と順（lattice の `table_responses`）
TWO_GROUP = [
    "total-1", "total-2", "nu-fission-1", "nu-fission-2", "fission-1", "fission-2",
    "scatter-1-1", "scatter-1-2", "scatter-2-1", "scatter-2-2",
]  # fmt: skip
# 暗黙の効きを反応度に直す数密度の変化の割合（ADR-0221 の perturbation）
FRACTION = 0.1


def allowed(direct: dict, i: int, value: float, tau: float, epsilon: float) -> float:
    """直接の摂動の Richardson の外挿 `direct`（`Steps`）の応答 i の許容"""
    coarse, fine = direct["central"][0][i], direct["central"][1][i]
    reference = direct["richardson"][i]
    return max(
        2.0 * abs(fine - coarse) / 3.0,
        RELATIVE * abs(reference),
        3.0 * tau * abs(value) / epsilon,
    )


def k_allowed(shielded: dict, tau: float, epsilon: float) -> float:
    """自己遮蔽を解き直した k の 1 次の変化と直接の摂動の比較の許容（両方の刻みの誤差、相対 1e-4、k の判定 τ が外挿に入る
    3 τ / ε の大きい方。compare_k_sensitivity.py と同じ形）"""

    def step(s: dict) -> float:
        coarse, fine = s["central"]
        return 2.0 * abs(fine - coarse) / 3.0

    return max(
        step(shielded["k_first_order"]),
        step(shielded["k_direct"]),
        RELATIVE * abs(shielded["k_direct"]["richardson"]),
        3.0 * tau / epsilon,
    )


def relative(value: float, reference: float) -> float:
    """`value` の `reference` に対する相対の差（両方 0 なら 0、参照だけ 0 なら無限大）"""
    if reference != 0.0:
        return (value - reference) / reference
    return 0.0 if value == 0.0 else float("inf")


def two_group_k(values: dict[str, float], chi: list[float]) -> float:
    """2 群の定数の無限媒質の k∞（`values` は応答の名前から値）"""
    a = np.array(
        [
            [
                values["total-1"] - values["scatter-1-1"],
                -values["scatter-2-1"],
            ],
            [
                -values["scatter-1-2"],
                values["total-2"] - values["scatter-2-2"],
            ],
        ]
    )
    nu = np.array([values["nu-fission-1"], values["nu-fission-2"]])
    return float(nu @ np.linalg.solve(a, np.array(chi)))


def k_gradient(values: dict[str, float], chi: list[float]) -> dict[str, float]:
    """k∞ の定数による微分（中心差分。k∞ は定数の有理式なので、相対 1e-6 の刻みで丸めより十分正確）"""
    gradient = {}
    for name, value in values.items():
        h = 1e-6 * abs(value) if value != 0.0 else 1e-12
        up, down = dict(values), dict(values)
        up[name], down[name] = value + h, value - h
        gradient[name] = (two_group_k(up, chi) - two_group_k(down, chi)) / (2.0 * h)
    return gradient


def compare(report: dict) -> dict:
    """I0 の比較の行、I1 と暗黙の効きの行"""
    few = report["few_group"]
    names = few["responses"]
    # 反応度に直す k∞ は 2 群の式で、比べる参照は直接の摂動
    if names != TWO_GROUP:
        sys.exit(f"the comparison takes the two-group table constants, got {names}")
    if "epsilon" not in report:
        sys.exit("the comparison needs the direct perturbations ([sensitivity] direct)")
    tau = report["flux_tolerance"]
    epsilon = report["epsilon"]
    k = report["k"]
    values = dict(zip(names, few["values"]))
    gradient = k_gradient(values, few["chi"])
    densities = {n["name"]: n["density"] for n in report["nuclides"]}
    k_direct = {n["name"]: n["direct"]["richardson"] for n in report["nuclides"]}
    explicit, implicit = [], []
    for n in few["nuclides"]:
        direct = n["direct"]
        for i, response in enumerate(names):
            g, d = n["relative"][i], direct["richardson"][i]
            limit = allowed(direct, i, few["values"][i], tau, epsilon)
            explicit.append(
                {
                    "name": n["name"],
                    "response": response,
                    "value": few["values"][i],
                    "adjoint": g,
                    "direct": d,
                    "relative": relative(g, d),
                    "allowed": limit,
                    "ok": abs(g - d) <= limit,
                }
            )
        shielded = n.get("shielded")
        if shielded is None:
            continue
        reshielded = shielded["direct"]
        first = shielded["first_order"]
        change = {
            response: reshielded["richardson"][i] - direct["richardson"][i]
            for i, response in enumerate(names)
        }
        rows = []
        for i, response in enumerate(names):
            s = reshielded["richardson"][i]
            # I1 も刻み ε と ε/2 の中心差分で、自己遮蔽の数密度への非線形の依存による刻みの誤差を持つ
            limit = max(
                allowed(reshielded, i, few["values"][i], tau, epsilon),
                allowed(first, i, few["values"][i], tau, epsilon),
            )
            rows.append(
                {
                    "response": response,
                    "implicit": change[response],
                    # 暗黙の分の、定数に対する割合と、陽の分に対する割合
                    "per_value": change[response] / few["values"][i]
                    if few["values"][i] != 0.0
                    else float("inf"),
                    "per_explicit": change[response] / direct["richardson"][i]
                    if direct["richardson"][i] != 0.0
                    else float("inf"),
                    "first_order": first["richardson"][i],
                    "reshielded": s,
                    "ok": abs(first["richardson"][i] - s) <= limit,
                    "allowed": limit,
                    "pcm": 1e5
                    * FRACTION
                    * gradient[response]
                    * change[response]
                    / k**2,
                }
            )
        implicit.append(
            {
                "name": n["name"],
                "density": densities[n["name"]],
                "rows": rows,
                "pcm": sum(r["pcm"] for r in rows),
                "k_pcm": 1e5
                * FRACTION
                * (shielded["k_direct"]["richardson"] - k_direct[n["name"]])
                / k**2,
                "k_ok": abs(
                    shielded["k_first_order"]["richardson"]
                    - shielded["k_direct"]["richardson"]
                )
                <= k_allowed(shielded, report["k_tolerance"], epsilon),
                "k_first_order_pcm": 1e5
                * FRACTION
                * (shielded["k_first_order"]["richardson"] - k_direct[n["name"]])
                / k**2,
            }
        )
    return {
        "explicit": explicit,
        "implicit": implicit,
        "k_infinity": two_group_k(values, few["chi"]),
    }


def pair_estimate(
    implicit: list[dict], relative: dict[str, float]
) -> list[tuple[str, float]]:
    """組の数密度の相対の変化 `relative` による、核種ごとの自己遮蔽の分の見込み (pcm)（効きの大きい順）"""
    terms = [(n["name"], n["k_pcm"] / FRACTION * relative[n["name"]]) for n in implicit]
    return sorted(terms, key=lambda t: -abs(t[1]))


def main(arguments: list[str]) -> int:
    relative = None
    if arguments[:1] == ["--pair"]:
        relative = tomllib.loads(Path(arguments[1]).read_text())["relative"]
        arguments = arguments[2:]
    failed = False
    for path in arguments:
        report = tomllib.loads(Path(path).read_text())
        few = report["few_group"]
        result = compare(report)
        print(
            f"{path}: k {report['k']:.10f}, 2-group k-infinity {result['k_infinity']:.10f}, "
            f"eps {report['epsilon']}, table difference {few['table_difference']:.1e}, "
            f"GMRES tolerance {few['tolerance']:.0e}, sweeps {sum(few['sweeps'])}"
        )
        rows = result["explicit"]
        outside = [r for r in rows if not r["ok"]]
        print(f"I0 (self-shielding fixed): {len(rows)} pairs, {len(outside)} outside")
        for r in outside:
            print(
                f"  OUTSIDE {r['name']:8s} {r['response']:13s} G {r['adjoint']:.6e} D {r['direct']:.6e} "
                f"rel {r['relative']:.2e} allowed {r['allowed'] / abs(r['direct']):.2e}"
            )
        significant = [r for r in rows if abs(r["direct"]) > 1e-6 * abs(r["value"])]
        if not significant:
            sys.exit(f"{path}: no derivative is larger than 1e-6 of its constant")
        worst = max(significant, key=lambda r: abs(r["relative"]))
        print(
            f"  largest relative difference where |D| > 1e-6 R: {worst['name']} {worst['response']} "
            f"{worst['relative']:.2e} ({len(significant)} pairs)"
        )
        worst = max(rows, key=lambda r: abs(r["relative"]))
        print(
            f"  largest relative difference: {worst['name']} {worst['response']} {worst['relative']:.2e} "
            f"(D / R {worst['direct'] / worst['value']:.1e})"
        )
        if result["implicit"]:
            print(
                "implicit (self-shielding solved again minus fixed), per 10 % density change; "
                "largest |implicit / R| response"
            )
            print(
                f"{'nuclide':8s} {'pcm (2g)':>9s} {'pcm (k)':>9s} {'pcm (I1 k)':>10s} {'response':13s} "
                f"{'impl/R':>10s} {'impl/expl':>10s} {'I1 vs I3':>9s}"
            )
            for n in sorted(result["implicit"], key=lambda n: -abs(n["k_pcm"])):
                top = max(n["rows"], key=lambda r: abs(r["per_value"]))
                bad = [r["response"] for r in n["rows"] if not r["ok"]]
                if not n["k_ok"]:
                    bad.append("k")
                failed |= bool(bad)
                print(
                    f"{n['name']:8s} {n['pcm']:9.3f} {n['k_pcm']:9.3f} {n['k_first_order_pcm']:10.3f} "
                    f"{top['response']:13s} {top['per_value']:10.2e} {top['per_explicit']:10.2e} "
                    f"{'ok' if not bad else 'OUTSIDE ' + ','.join(bad):>9s}"
                )
            if relative is not None:
                terms = pair_estimate(result["implicit"], relative)
                total = sum(t for _, t in terms)
                listed = ", ".join(f"{n} {t:+.2f}" for n, t in terms[:5])
                print(
                    f"self-shielding part of the pair (first order): {total:+.2f} pcm; largest {listed}"
                )
        failed |= bool(outside)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
