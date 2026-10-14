"""OpenMC のピンの核分裂率を、VERA の仕様書の KENO-VI の分布（1/8 の表）と比べる。

使い方: python3 compare_keno.py <spec.txt> <問題（2B など）> <OpenMC の pins.msgpack>

仕様書の表は、集合体の中央から k 番目（1〜8）の行に、中央の列から k 番目までの値を並べたもので、案内管と計装管の位置は
空欄になる。OpenMC の分布は燃料ピンの平均で 1 にし、対称な 8 本を平均して表の位置と比べる（KENO-VI も 1/4 の計算を
対称な 8 本で平均している）。
"""

import math
import re
import sys
from pathlib import Path

import msgpack

CAPTIONS = {
    "2A": "Figure P2-3: Problem 2A (565K) CE KENO-VI Power Distribution Results",
    "2B": "Figure P2-4: Problem 2B (600K) CE KENO-VI Power Distribution Results",
    "2C": "Figure P2-5: Problem 2C (900K) CE KENO-VI Power Distribution Results",
    "2D": "Figure P2-6: Problem 2D (1200K) CE KENO-VI Power Distribution Results",
    "2E": "Figure P2-7: Problem 2E (12 Pyrex) CE KENO-VI Power Distribution Results",
    "2F": "Figure P2-8: Problem 2F (24 Pyrex) CE KENO-VI Power Distribution Results",
    "2G": "Figure P2-9: Problem 2G (AIC) CE KENO-VI Power Distribution Results",
    "2H": "Figure P2-10: Problem 2H (B4C) CE KENO-VI Power Distribution Results",
    "2I": "Figure P2-11: Problem 2I (Instrument Thimble) CE KENO-VI Power Distribution Results",
    "2J": "Figure P2-12: Problem 2J (Instrument and Pyrex) CE KENO-VI Power Distribution Results",
    "2K": "Figure P2-13: Problem 2K (Zoned and 24 Pyrex) CE KENO-VI Power Distribution Results",
    "2Q": "Figure P2-19: Problem 2Q (Zircaloy Spacer Grid) CE KENO-VI Power Distribution Results",
}
# 1/8 の中の案内管の位置（中央からの行と列の数）。中央 (0, 0) は計装管
TUBES = {(0, 0), (3, 0), (3, 3), (5, 5), (6, 0), (6, 3)}
VALUE = re.compile(r"\d\.\d{4}")


def keno_octant(spec: list[str], caption: str) -> dict[tuple[int, int], float]:
    """図の見出しの上にある 8 行の表（左半分が核分裂率、右半分が統計誤差）を、(行, 列) ごとの値にする。"""
    lines = [i for i, line in enumerate(spec) if caption in line]
    if not lines:
        raise SystemExit(f"caption not found: {caption}")
    rows: list[str] = []
    i = lines[-1]
    while len(rows) < 8:
        i -= 1
        if i < 0:
            raise SystemExit(f"table above {caption!r} is incomplete")
        line = spec[i]
        if line.strip() and "Max:" not in line:
            rows.append(line)
    octant = {}
    for k, line in enumerate(reversed(rows), start=1):
        values = [float(v) for v in VALUE.findall(line.split("%")[0])]
        columns = [j for j in range(k + 1) if (k, j) not in TUBES]
        if len(values) != len(columns):
            raise SystemExit(f"row {k} has {len(values)} values for columns {columns}")
        octant.update({(k, j): v for j, v in zip(columns, values)})
    return octant


def compare(
    spec: list[str], problem: str, rows: list[list[float]]
) -> tuple[int, float, float, tuple[int, int]]:
    """対称な 8 本の平均を KENO-VI の表と比べる。位置の数、相対差の二乗平均、絶対値が最大の相対差とその位置を返す。"""
    if problem not in CAPTIONS:
        raise SystemExit(f"unknown problem {problem}; expected one of {list(CAPTIONS)}")
    if len(rows) != 17 or any(len(row) != 17 for row in rows):
        raise SystemExit("expected a 17 x 17 map")
    fuel = [x for row in rows for x in row if x > 0]
    mean = sum(fuel) / len(fuel)
    differences = []
    for (i, j), keno in keno_octant(spec, CAPTIONS[problem]).items():
        images = {
            (8 + a, 8 + b)
            for a, b in [
                (i, j),
                (j, i),
                (-i, j),
                (j, -i),
                (i, -j),
                (-j, i),
                (-i, -j),
                (-j, -i),
            ]
        }
        ours = sum(rows[r][c] for r, c in images) / len(images) / mean
        differences.append((ours / keno - 1, (i, j)))
    rms = math.sqrt(sum(d * d for d, _ in differences) / len(differences))
    worst, at = max(differences, key=lambda d: abs(d[0]))
    return len(differences), rms, worst, at


def main() -> None:
    spec = Path(sys.argv[1]).read_text().split("\n")
    problem = sys.argv[2]
    rows = msgpack.unpackb(Path(sys.argv[3]).read_bytes())["rows"]
    n, rms, worst, at = compare(spec, problem, rows)
    print(
        f"{problem}: {n} positions, rms {100 * rms:.3f}%, max {100 * worst:+.3f}% at {at}"
    )


if __name__ == "__main__":
    main()
