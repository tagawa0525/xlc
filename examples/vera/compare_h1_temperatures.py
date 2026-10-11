"""565 K で処理した H-1 の自由原子の断面積と、550 K と 573.6 K の断面積の温度についての線形の内挿を比べる。

使い方: python3 compare_h1_temperatures.py <OpenMC の H1.h5>

H1.h5 は openmc_library.py で、565 K と S(α,β) の温度（550 K と 573.6 K）を含めて作ったもの。参照解の OpenMC は
565 K の H-1 を使い、xlc は 550 K と 573.6 K の H-1 の表全体を温度で内挿する（ADR-0080）。その違いの上限として、
565 K のエネルギーの格子の点ごとに、全断面積、弾性散乱、捕獲の相対差（内挿 / 565 K − 1）の最大をエネルギーの
範囲ごとに表示する。
"""

import sys

import numpy as np
import openmc.data

# 比べる反応（MT）
REACTIONS = {"total": 1, "elastic": 2, "capture": 102}
# エネルギーの範囲 (eV)。4 eV は H in H2O の S(α,β) を使う上限の目安
RANGES = [(1e-5, 4.0), (4.0, 1e5), (1e5, 2e7)]


def main() -> None:
    data = openmc.data.IncidentNeutron.from_hdf5(sys.argv[1])
    temperatures = {
        name: kt / openmc.data.K_BOLTZMANN
        for name, kt in zip(data.temperatures, data.kTs, strict=True)
    }
    for name in ("550K", "565K", "574K"):
        if name not in temperatures:
            raise SystemExit(
                f"{sys.argv[1]} has no {name}; it has {list(temperatures)}"
            )
    low, high = temperatures["550K"], temperatures["574K"]
    f = (temperatures["565K"] - low) / (high - low)
    energy = data.energy["565K"]
    print(f"temperatures {temperatures}; weight of 574K {f:.6f}")
    for label, mt in REACTIONS.items():
        xs = data[mt].xs
        exact = xs["565K"](energy)
        interpolated = (1 - f) * xs["550K"](energy) + f * xs["574K"](energy)
        relative = interpolated / exact - 1
        for lo, hi in RANGES:
            inside = (energy >= lo) & (energy < hi)
            worst = np.max(np.abs(relative[inside]))
            print(
                f"{label:8s} {lo:8.0e}–{hi:8.0e} eV: max |relative difference| {worst:.2e}"
            )


if __name__ == "__main__":
    main()
