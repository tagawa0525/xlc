"""xlc と OpenMC の、集合体（VERA 2B〜2D）のピンごとの核分裂率を比べる。

使い方: python3 compare_pins.py <xlc の pins.msgpack> <OpenMC の pins.msgpack>

どちらもセルごとの核分裂率（17 × 17、上の行から）。燃料ピン（核分裂率が正のセル）の平均で規格化し直し、
相対差 (xlc / OpenMC − 1) の二乗平均と最大を表示する。
"""

import math
import sys
from pathlib import Path

import msgpack


def read(path: str) -> list[list[float]]:
    rows = msgpack.unpackb(Path(path).read_bytes())["rows"]
    values = [x for row in rows for x in row if x > 0]
    mean = sum(values) / len(values)
    return [[x / mean for x in row] for row in rows]


def main() -> None:
    ours, reference = read(sys.argv[1]), read(sys.argv[2])
    diffs = [
        (a / b - 1, i, j)
        for i, (row_a, row_b) in enumerate(zip(ours, reference))
        for j, (a, b) in enumerate(zip(row_a, row_b))
        if a > 0 and b > 0
    ]
    rms = math.sqrt(sum(d * d for d, _, _ in diffs) / len(diffs))
    worst = max(diffs, key=lambda x: abs(x[0]))
    print(
        f"{len(diffs)} pins: rms {100 * rms:.2f}%, max {100 * worst[0]:+.2f}% at {worst[1:]}"
    )


if __name__ == "__main__":
    main()
