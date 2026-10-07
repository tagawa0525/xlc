"""xlc と OpenMC の、集合体（VERA 2B〜2D）のピンごとの核分裂率を比べる。

使い方: python3 compare_pins.py <xlc の pins.msgpack> <OpenMC の pins.msgpack>

どちらもセルごとの核分裂率（17 × 17、上の行から）。燃料ピン（核分裂率が正のセル）の平均で規格化し直し、
相対差 (xlc / OpenMC − 1) の二乗平均と最大、水の穴（燃料のないセル）の上下左右の隣、斜めの隣、それ以外の
ピンの平均を表示する。
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
    shape = [len(row) for row in ours]
    if shape != [len(row) for row in reference]:
        raise SystemExit("the two maps have different shapes")
    fuel = [[x > 0 for x in row] for row in ours]
    if fuel != [[x > 0 for x in row] for row in reference]:
        raise SystemExit("the two maps have fission in different cells")
    diffs = [
        (a / b - 1, i, j)
        for i, (row_a, row_b) in enumerate(zip(ours, reference))
        for j, (a, b) in enumerate(zip(row_a, row_b))
        if a > 0 and b > 0
    ]
    n = len(reference)
    holes = {(i, j) for i in range(n) for j in range(n) if not fuel[i][j]}

    def touches(i: int, j: int, steps: tuple[tuple[int, int], ...]) -> bool:
        return any((i + di, j + dj) in holes for di, dj in steps)

    side = ((1, 0), (-1, 0), (0, 1), (0, -1))
    corner = ((1, 1), (1, -1), (-1, 1), (-1, -1))
    groups = {"beside holes": [], "diagonal to holes": [], "others": []}
    for d, i, j in diffs:
        if touches(i, j, side):
            groups["beside holes"].append(d)
        elif touches(i, j, corner):
            groups["diagonal to holes"].append(d)
        else:
            groups["others"].append(d)
    rms = math.sqrt(sum(d * d for d, _, _ in diffs) / len(diffs))
    worst = max(diffs, key=lambda x: abs(x[0]))
    print(
        f"{len(diffs)} pins: rms {100 * rms:.2f}%, max {100 * worst[0]:+.2f}% at {worst[1:]}"
    )
    for name, values in groups.items():
        if values:
            mean = 100 * sum(values) / len(values)
            print(f"{name} ({len(values)} pins): mean {mean:+.2f}%")


if __name__ == "__main__":
    main()
