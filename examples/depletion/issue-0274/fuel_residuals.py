"""直積の分岐のテーブルの燃料の温度の軸で、2 群の定数が √T の滑らかな曲線から外れる大きさを見る（#0274）。

使い方: python3 fuel_residuals.py <直積のテーブル>

燃焼度、ホウ素、減速材の点の組ごとに、燃料の温度のすべての点の定数（高速群と熱群の吸収と νΣf）を √T の 3 次式に最小 2 乗で
当てはめ、残差を点の値の平均で割った値（1e-5 単位）の、組の中の最小と最大を燃料の温度の点ごとに出す。
"""

import itertools
import sys
from pathlib import Path

import msgpack
import numpy as np


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    t = msgpack.unpackb(Path(sys.argv[1]).read_bytes())
    names = [a["name"] for a in t["axes"]]
    if names != ["burnup", "boron", "moderator", "fuel"]:
        raise SystemExit(
            f"the table needs the axes burnup, boron, moderator, fuel, not {names}"
        )
    points = [a["points"] for a in t["axes"]]
    indices = list(itertools.product(*[range(len(p)) for p in points]))
    if len(indices) != len(t["states"]):
        raise SystemExit("the table must hold the product of its axes")
    flat = {x: n for n, x in enumerate(indices)}
    g = np.sqrt(points[3])
    print("fuel (K): " + " ".join(f"{f:>13g}" for f in points[3]))
    for key, group in (
        ("absorption", 0),
        ("absorption", 1),
        ("nu_fission", 0),
        ("nu_fission", 1),
    ):
        residuals = []
        for rest in itertools.product(*[range(len(p)) for p in points[:3]]):
            y = np.array(
                [t["states"][flat[(*rest, f)]][key][group] for f in range(len(g))]
            )
            fit = np.polyval(np.polyfit(g, y, 3), g)
            residuals.append((y - fit) / y.mean() * 1e5)
        r = np.array(residuals)
        cells = " ".join(
            f"{lo:+6.1f}〜{hi:+6.1f}" for lo, hi in zip(r.min(0), r.max(0), strict=True)
        )
        print(f"{key} {group + 1} 群（{len(residuals)} 組、1e-5）: {cells}")


if __name__ == "__main__":
    main()
