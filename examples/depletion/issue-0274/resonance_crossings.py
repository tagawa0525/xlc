"""分離共鳴を持つ核種の自由ガスの上限（kT の倍数。ADR-0204）が、燃料の温度の軸（ADR-0239）の 560〜1350 K で越える捕獲の極大を
並べる（#0274）。

使い方: python3 resonance_crossings.py <核種の出力> <kT の倍数>...

核種の出力（xs-nuclides/<核種>.msgpack）の 560 K のポイントワイズの捕獲の極大（両隣の点より大きい点）ごとに、上限がその
エネルギーになる温度 E / (倍数 k) が 560〜1350 K に入るものを、エネルギー、極大の捕獲、温度の順に出す。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

BOLTZMANN = 8.617333262e-5  # eV/K
LOW, HIGH = 560.0, 1350.0


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    d = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
    (nuclide,) = d["nuclides"]
    (table,) = [t for t in nuclide["temperatures"] if t["temperature"] == LOW]
    e = np.array(table["pointwise"]["energy"])
    c = np.array(table["pointwise"]["capture"])
    peaks = np.flatnonzero((c[1:-1] > c[:-2]) & (c[1:-1] > c[2:])) + 1
    for multiple in map(float, sys.argv[2:]):
        print(f"{nuclide['name']}、上限 {multiple:g} kT:")
        for i in peaks:
            t = e[i] / (multiple * BOLTZMANN)
            if LOW <= t <= HIGH:
                print(f"  {e[i]:9.3f} eV  {c[i]:8.1f} b  {t:7.1f} K")


if __name__ == "__main__":
    main()
