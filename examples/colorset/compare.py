"""2 つのピン出力（`xlc-pin-powers`、行は上から）を比べ、相対の差 (a / b − 1) の RMS、絶対値の最大、b の出力が
0.8 以上のピンの絶対値の最大を出す。どちらも燃料のピン（出力が正のセル）の平均で 1 に規格化し直す。

`--blocks n` を付けると、行と列を n × n の集合体に分け、差を集合体の平均の出力の差と、集合体の中の形（集合体の平均で
割った値）の差に分けても出す。形の差の「出力 0.8 以上」の判定には、形ではなく b の出力そのものを使う（ADR-0118 の判定と
同じピンを選ぶ）。地図が正方形で n で割り切れなければエラーにする。

使い方: python3 compare.py <a.msgpack> <b（参照）.msgpack> [--blocks n]
"""

import sys
from pathlib import Path

import msgpack
import numpy as np


def read(path: str) -> np.ndarray:
    rows = np.array(msgpack.unpackb(Path(path).read_bytes(), raw=False)["rows"])
    rows = rows[rows.any(axis=1)][:, rows.any(axis=0)]
    return rows / rows[rows > 0].mean()


def summary(name: str, a: np.ndarray, b: np.ndarray, power: np.ndarray) -> None:
    """a / b − 1 をまとめる。「出力 0.8 以上」は `power`（参照の出力）で選ぶ。"""
    fuel = b > 0
    if not np.array_equal(fuel, a > 0):
        raise SystemExit("the two maps have different fuel cells")
    d = a[fuel] / b[fuel] - 1
    hot = power[fuel] >= 0.8
    if not hot.any():
        raise SystemExit("no pin of the reference has a power of 0.8 or more")
    print(
        f"{name}: {fuel.sum()} pins, rms {100 * np.sqrt((d * d).mean()):.3f}%, "
        f"max {100 * d[np.abs(d).argmax()]:+.3f}%, max at power >= 0.8 {100 * np.abs(d[hot]).max():.3f}%"
    )


def main() -> None:
    a, b = read(sys.argv[1]), read(sys.argv[2])
    if a.shape != b.shape:
        raise SystemExit(f"shapes differ: {a.shape} and {b.shape}")
    summary("pins", a, b, b)
    if "--blocks" in sys.argv:
        n = int(sys.argv[sys.argv.index("--blocks") + 1])
        if a.shape[0] != a.shape[1] or a.shape[0] % n != 0:
            raise SystemExit(
                f"a {a.shape} map does not split into {n} x {n} assemblies"
            )
        size = a.shape[0] // n
        mean_a, mean_b = np.ones_like(a), np.ones_like(b)
        blocks = []
        for r in range(n):
            for c in range(n):
                s = (slice(r * size, (r + 1) * size), slice(c * size, (c + 1) * size))
                fa, fb = a[s][a[s] > 0].mean(), b[s][b[s] > 0].mean()
                mean_a[s], mean_b[s] = fa, fb
                blocks.append(f"{100 * (fa / fb - 1):+.3f}%")
        print("assembly means:", " ".join(blocks))
        summary("shape in assemblies", a / mean_a, b / mean_b, b)


if __name__ == "__main__":
    main()
