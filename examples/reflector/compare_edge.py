"""2 つのピン出力（`xlc-pin-powers`）の相対の差 (a / b − 1) を、燃料のピン全体と、右の端の列、上の端の行、
右上の角のピンに分けてまとめる（I 型の格子では右の列、角型の格子では右の列と上の行が反射体に接する）。どちらも燃料のピン（出力が正のセル）の平均で 1 に規格化し直す（../colorset/compare.py と同じ）。

使い方: python3 compare_edge.py <a.msgpack> <b（参照）.msgpack>
"""

import sys
from pathlib import Path

import msgpack
import numpy as np


def read(path: str) -> np.ndarray:
    rows = np.array(msgpack.unpackb(Path(path).read_bytes(), raw=False)["rows"])
    rows = rows[rows.any(axis=1)][:, rows.any(axis=0)]
    return rows / rows[rows > 0].mean()


def summary(name: str, d: np.ndarray) -> None:
    print(
        f"{name}: {d.size} pins, rms {100 * np.sqrt((d * d).mean()):.3f}%, max {100 * d[np.abs(d).argmax()]:+.3f}%"
    )


def main() -> None:
    a, b = read(sys.argv[1]), read(sys.argv[2])
    fuel = b > 0
    if not np.array_equal(fuel, a > 0):
        raise SystemExit("the two maps have different fuel cells")
    d = np.zeros_like(b)
    d[fuel] = a[fuel] / b[fuel] - 1
    summary("all", d[fuel])
    edge = fuel[:, -1]
    summary(
        f"right column (power {b[edge, -1].min():.3f}–{b[edge, -1].max():.3f})",
        d[edge, -1],
    )
    top = fuel[0, :]
    summary(
        f"top row (power {b[0, top].min():.3f}–{b[0, top].max():.3f})",
        d[0, top],
    )
    print(f"top right pin: power {b[0, -1]:.3f}, {100 * d[0, -1]:+.3f}%")


if __name__ == "__main__":
    main()
