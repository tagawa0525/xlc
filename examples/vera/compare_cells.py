"""燃料ピンを格子の入力のセルの名前（2K の 3.1% の `pin` と 3.6% の `pin36` など）で分け、セルの種類ごとに、xlc と
OpenMC の核分裂率の差の平均と、基準の計算からの変化の平均と、xlc と OpenMC の
変化の比の平均を、OpenMC の統計誤差つきで表示する。

使い方: python3 compare_cells.py <格子の入力> <xlc の pins.msgpack> <OpenMC の pins.msgpack> <OpenMC の作業ディレクトリ>
        [<基準の xlc の pins.msgpack> <基準の OpenMC の pins.msgpack> <基準の OpenMC の作業ディレクトリ>]

核分裂率は、どの分布も燃料ピンの平均で 1 にする（変化は集合体の平均に対する相対の値）。セルの位置は格子の入力の
`rows` から読む。統計誤差の扱いは compare_neighbours.py と同じで、ピンと世代を独立とした下限の見積もり。
"""

import sys
from pathlib import Path

import numpy as np
import tomllib
from compare_neighbours import read_pins, relative_sigma
from openmc_model import N


def cell_names(path: str) -> list[list[str]]:
    rows = [row.split() for row in tomllib.loads(Path(path).read_text())["rows"]]
    if len(rows) != N or any(len(row) != N for row in rows):
        raise SystemExit(f"{path}: expected {N} x {N} rows")
    return rows


def summary(name: str, ratios: list[float], sigmas: list[float]) -> None:
    """比の平均、最小、最大と、ピンと世代を独立とした平均の統計誤差（下限）。"""
    mean = np.mean(ratios) - 1
    sigma = np.sqrt(np.sum(np.square(sigmas))) / len(sigmas)
    print(
        f"{name}: mean {100 * mean:+.3f} % (min {100 * (min(ratios) - 1):+.3f}, "
        f"max {100 * (max(ratios) - 1):+.3f}), "
        f"sigma assuming independent pins and generations {100 * sigma:.3f} %"
    )


def main() -> None:
    if len(sys.argv) not in (5, 8):
        raise SystemExit(__doc__)
    names = cell_names(sys.argv[1])
    xlc = read_pins(sys.argv[2])
    pins, sigma = read_pins(sys.argv[3]), relative_sigma(sys.argv[4])
    if not np.array_equal(xlc > 0, pins > 0):
        raise SystemExit("xlc and OpenMC have fission in different cells")
    base = None
    if len(sys.argv) == 8:
        base_xlc = read_pins(sys.argv[5])
        base, base_sigma = read_pins(sys.argv[6]), relative_sigma(sys.argv[7])
        if not np.array_equal(base_xlc > 0, pins > 0) or not np.array_equal(
            base > 0, pins > 0
        ):
            raise SystemExit("the base maps have fission in different cells")
    by_name: dict[str, list[tuple[int, int]]] = {}
    for i in range(N):
        for j in range(N):
            if pins[i, j] > 0:
                by_name.setdefault(names[i][j], []).append((i, j))
    for name, cells in sorted(by_name.items()):
        print(f"{name}: {len(cells)} pins")
        summary(
            "  xlc / OpenMC - 1",
            [xlc[p] / pins[p] for p in cells],
            [sigma[p] for p in cells],
        )
        if base is not None:
            summary(
                "  OpenMC / base - 1",
                [pins[p] / base[p] for p in cells],
                [np.hypot(sigma[p], base_sigma[p]) for p in cells],
            )
            change = np.mean([xlc[p] / base_xlc[p] for p in cells]) - 1
            print(f"  xlc / base - 1: mean {100 * change:+.3f} %")
            # xlc の分布は誤差を持たないとするので、変化の比の統計誤差は OpenMC の 2 つの計算のものだけ
            summary(
                "  (xlc / base) / (OpenMC / base) - 1",
                [(xlc[p] / base_xlc[p]) / (pins[p] / base[p]) for p in cells],
                [np.hypot(sigma[p], base_sigma[p]) for p in cells],
            )


if __name__ == "__main__":
    main()
