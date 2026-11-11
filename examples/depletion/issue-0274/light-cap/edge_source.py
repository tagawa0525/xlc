"""自由ガスの核を入射エネルギー E_c で静止した標的に切り替えたときの、減速の散乱源の形（ADR-0246、#0274）。

使い方: python3 edge_source.py

断面積一定の自由ガスの核（lattice の `free_gas::kernel` と同じ閉じた式の、下方散乱の側）と静止した標的の核で、衝突密度が 1/E' の
一様な減速（1 つの核種）の散乱源 S(E) = ∫ K(E'→E) / E' dE' を数値積分し、次を出す:

- O-16（燃料 900 K）で、上限 E_c の核（E_c より下の入射は自由ガス、上は静止した標的）の S(E) と、上限のない核の S(E) の比 − 1 を、
  出射 E = α E_c の前後で
- H-1（減速材 580 K）の 1 回の散乱で 0.625 eV より下に落ちる割合の、自由ガスと静止した標的の比 − 1

積分は入射 E' ≥ E に限る（上方散乱の側の核はこの式では扱わない）。上限が E と E/α の間の出射では、上方散乱の側の入射は
上限より下でどちらの核も自由ガスなので、比 − 1 は下端の側の差だけで決まる。
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfcx

BOLTZMANN = 8.617333262e-5  # eV/K
O16 = 15.857510  # ENDF/B-VIII.0 の O-16 の AWR
H1 = 0.999167


def kernel(e: float, e_out: float, awr: float, kt: float) -> float:
    """断面積一定の自由ガスの核 S(E→E')（自由原子の断面積あたり、1/eV）。E' ≤ E の側だけ"""
    if e_out > e:
        raise ValueError("only the down-scattering side")
    theta = (awr + 1) / (2 * np.sqrt(awr))
    zeta = (awr - 1) / (2 * np.sqrt(awr))
    s, s_out = np.sqrt(e / kt), np.sqrt(e_out / kt)
    x1, x2 = theta * s + zeta * s_out, theta * s - zeta * s_out
    u, v = zeta * s - theta * s_out, zeta * s + theta * s_out
    if u > 0:
        # 遠い下方散乱: 同じ指数因子の項をまとめて打ち消し合いを避ける（lattice と同じ）
        bracket = np.exp(-u * u) * (erfcx(u) - erfcx(x2)) + np.exp(-v * v) * (
            erfcx(x1) - erfcx(v)
        )
    else:
        bracket = (
            erf(-u) + erf(v) + np.exp(-v * v) * erfcx(x1) - np.exp(-u * u) * erfcx(x2)
        )
    return max(theta * theta / (2 * e) * bracket, 0.0)


def static(e: float, e_out: float, awr: float) -> float:
    alpha = ((awr - 1) / (awr + 1)) ** 2
    return 1 / ((1 - alpha) * e) if alpha * e <= e_out <= e else 0.0


def source(e_out: float, cap: float, awr: float, kt: float) -> float:
    """入射 E_c = `cap` で核を切り替えた散乱源（衝突密度 1/E'）"""
    alpha = ((awr - 1) / (awr + 1)) ** 2
    lo, hi = e_out, e_out / alpha * 1.15
    edge = e_out / alpha

    def part(a: float, b: float, free: bool) -> float:
        if b <= a:
            return 0.0
        f = (
            (lambda e: kernel(e, e_out, awr, kt) / e)
            if free
            else (lambda e: static(e, e_out, awr) / e)
        )
        return quad(f, a, b, points=[edge] if a < edge < b else None, limit=400)[0]

    c = min(max(cap, lo), hi)
    return part(lo, c, True) + part(c, hi, False)


def thermal_fraction(e: float, threshold: float, awr: float, kt: float) -> float:
    """入射 `e` の自由ガスの核で `threshold` より下に落ちる割合"""
    return quad(
        lambda e_out: kernel(e, e_out, awr, kt),
        1e-8,
        threshold,
        points=[kt, 5 * kt],
        limit=400,
    )[0]


def main() -> None:
    kt = BOLTZMANN * 900.0
    alpha = ((O16 - 1) / (O16 + 1)) ** 2
    cap = 27.14
    print(
        f"O-16、900 K、E_c = {cap} eV（α E_c = {alpha * cap:.3f} eV）: S(E)/S_上限なし(E) − 1"
    )
    for e in np.arange(alpha * cap - 2.0, alpha * cap + 2.01, 0.25):
        r = source(e, cap, O16, kt) / source(e, np.inf, O16, kt) - 1
        print(f"  {e:7.3f} eV {r * 100:+7.3f} %")
    kt = BOLTZMANN * 580.0
    print("H-1、580 K: 0.625 eV より下に落ちる割合の、自由ガス / 静止した標的 − 1")
    for e in [20.0, 50.0, 200.0]:
        f = thermal_fraction(e, 0.625, H1, kt)
        print(f"  {e:6.1f} eV {(f / (0.625 / e) - 1) * 100:+7.3f} %")


if __name__ == "__main__":
    main()
