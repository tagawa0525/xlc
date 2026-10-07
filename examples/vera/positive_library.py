"""多群ライブラリの散乱を、P_L で打ち切っても負にならない形にする（OpenMC の多群モードとの比較用）。

使い方: python3 positive_library.py <ライブラリ (msgpack)> <出力 (msgpack)> <次数 L>

打ち切ったルジャンドル展開 f(μ) = Σ_l (2l+1)/2 σ_l P_l(μ) / σ_0 は、水素などで負になる。OpenMC の多群
モードは負の部分を 0 に切って規格化し直す（src/scattdata.cpp）ので、モーメントの違う散乱核を解くことになる。
そこで核種・温度・群の移行ごとに、P1〜P_L のモーメントに同じ係数 λ（0〜1）を掛け、f が [-1, 1] で負に
ならない最大の λ にする（P_L より上は捨てる。f の最小は、端と導関数の根で正確に求める）。物質の散乱核は核種の正の重みの和なので、これも負にならず、
xlc と OpenMC が同じ散乱核を解く。物理的な散乱核ではなく、輸送の解き方の検証のためだけに使う。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
from numpy.polynomial import legendre


def minimum(coefficients: np.ndarray) -> float:
    """ルジャンドル級数 Σ_l c_l P_l(μ) の [-1, 1] での最小値（端と、導関数の実の根で求める）。"""
    roots = (
        legendre.legroots(legendre.legder(coefficients))
        if len(coefficients) > 2
        else []
    )
    points = [-1.0, 1.0] + [
        r.real
        for r in np.atleast_1d(roots)
        if abs(r.imag) < 1e-12 and -1 <= r.real <= 1
    ]
    return float(np.min(legendre.legval(np.array(points), coefficients)))


def positive(scatter: list, order: int) -> list:
    """[l][g][g'] の散乱行列を、P_order までの負にならない形にする。"""
    s = np.array(scatter[: order + 1])
    if s.shape[0] < order + 1:
        raise SystemExit(f"the library has scattering moments up to P{s.shape[0] - 1}")
    # 負の P0 は物質の散乱核を負にしうる（モーメントを縮めても直らない）
    if (s[0] < 0).any():
        raise SystemExit("the library has negative P0 transfers")
    # P0 が 0 の移行に高次のモーメントだけがあると、必ずどこかで負になる
    s[1:, s[0] == 0] = 0.0
    weights = (2 * np.arange(order + 1) + 1) / 2
    for g, h in zip(*np.nonzero(s[0] > 0)):
        # f(μ) = 1/2 + λ h(μ)、h は P1〜P_L の項
        coefficients = weights * s[:, g, h] / s[0, g, h]
        coefficients[0] = 0.0
        lowest = minimum(coefficients)
        if lowest < -0.5:
            s[1:, g, h] *= 0.5 / -lowest
    return s.tolist()


def main() -> None:
    source, output, order = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    library = msgpack.unpackb(source.read_bytes())
    if library.get("format") != "xlc-multigroup":
        raise SystemExit(f"{source}: not an xlc-multigroup library")
    for nuclide in library["nuclides"]:
        for table in nuclide["temperatures"]:
            table["scatter"] = positive(table["scatter"], order)
    output.write_bytes(msgpack.packb(library))


if __name__ == "__main__":
    main()
