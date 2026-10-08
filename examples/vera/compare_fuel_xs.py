"""xlc と OpenMC の、燃料の核種ごとの実効断面積を群ごとに比べる（自己遮蔽の検証）。

使い方: python3 compare_fuel_xs.py <xlc の fuel-xs.msgpack> <OpenMC の fuel-xs.msgpack> [核種 ...]

xlc 側は crates/lattice/examples/fuel_xs.rs、OpenMC 側は openmc_model.py で作る。核種と反応ごとに、
相対差 (xlc / OpenMC − 1) が最も大きい群と、共鳴の群（WIMS 69 群の 15〜27 群）の差を表示する。
核種は OpenMC の名前（U238 など）で、省略すると U238 と U235。
"""

import sys
from pathlib import Path

import msgpack

# WIMS 69 群の分離共鳴の群（1 始まり）
RESONANCE = range(15, 28)


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())["nuclides"]


def main() -> None:
    ours, reference = read(sys.argv[1]), read(sys.argv[2])
    for side in (ours, reference):
        for rates in side.values():
            # 非弾性（全 − 弾性）。共鳴の群では捕獲 + 核分裂（吸収）に等しい
            rates["nonelastic"] = [
                t - e for t, e in zip(rates["total"], rates["elastic"])
            ]
    nuclides = sys.argv[3:] or ["U238", "U235"]
    for nuclide in nuclides:
        for score in ["total", "elastic", "nonelastic", "fission"]:
            a, b = ours[nuclide][score], reference[nuclide][score]
            diffs = [(g + 1, x / y - 1) for g, (x, y) in enumerate(zip(a, b)) if y > 0]
            if not diffs:
                continue
            worst = max(diffs, key=lambda d: abs(d[1]))
            resonance = " ".join(f"{d * 100:+.1f}" for g, d in diffs if g in RESONANCE)
            print(f"{nuclide} {score}: worst group {worst[0]} {worst[1] * 100:+.2f}%")
            print(f"  groups 15-27 (%): {resonance}")


if __name__ == "__main__":
    main()
