"""多群ライブラリの散乱を、P_L で打ち切っても負にならない形にする（OpenMC の多群モードとの比較用）。

使い方: python3 positive_library.py <ライブラリ (msgpack)> <出力 (msgpack)> <次数 L>

打ち切ったルジャンドル展開 f(μ) = Σ_l (2l+1)/2 σ_l P_l(μ) / σ_0 は、水素などで負になる。OpenMC の多群
モードは負の部分を 0 に切って規格化し直す（src/scattdata.cpp）ので、モーメントの違う散乱核を解くことになる。
そこで核種・温度・群の移行ごとに、P1〜P_L のモーメントに同じ係数 λ（0〜1）を掛け、f が [-1, 1] で負に
ならない最大の λ にする（P_L より上は捨てる）。物質の散乱核は核種の正の重みの和なので、これも負にならず、
xlc と OpenMC が同じ散乱核を解く。物理的な散乱核ではなく、輸送の解き方の検証のためだけに使う。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
from numpy.polynomial import legendre

# 負でないことを確かめる μ の点（λ はこの点で 0 以上になるように決め、点の間の小さな負は許す）
MU = np.linspace(-1.0, 1.0, 4001)


def positive(scatter: list, order: int) -> list:
    """[l][g][g'] の散乱行列を、P_order までの負にならない形にする。"""
    s = np.array(scatter[: order + 1])
    if s.shape[0] < order + 1:
        raise SystemExit(f"the library has scattering moments up to P{s.shape[0] - 1}")
    # [l][μ] の (2l+1)/2 P_l(μ)
    basis = np.array(
        [(2 * l + 1) / 2 * legendre.legval(MU, [0] * l + [1]) for l in range(order + 1)]
    )
    # P0 が 0 の移行に高次のモーメントだけがあると、必ずどこかで負になる
    s[1:, s[0] <= 0] = 0.0
    for g, h in zip(*np.nonzero(s[0] > 0)):
        # f(μ) = 1/2 + λ h(μ)
        higher = (s[1:, g, h] / s[0, g, h]) @ basis[1:]
        negative = higher < 0
        if negative.any():
            factor = min(1.0, float(np.min(0.5 / -higher[negative])))
            s[1:, g, h] *= factor
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
