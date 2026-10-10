"""セルの上下左右の隣の 4 本のピンについて、xlc と OpenMC の核分裂率の差と、OpenMC の計算どうしの変化を、
OpenMC の統計誤差つきで表示する（計装管のシンブルのような、1 つのセルに入れたものの影響を見るため）。

使い方: python3 compare_neighbours.py <行> <列> <xlc の pins.msgpack> <OpenMC の pins.msgpack> <OpenMC の作業ディレクトリ>
        [<基準の OpenMC の pins.msgpack> <基準の OpenMC の作業ディレクトリ>]

行と列は 0 始まりで上の行から（lattice のピン出力と同じ向き）。核分裂率は pins.msgpack（openmc_model.py の書くもの、
p2/openmc-*.msgpack）から読み、燃料ピン（核分裂率が正のセル）の平均で規格化する。作業ディレクトリ（statepoint の番号が
最も大きいものを読む）からは、名前が fission の tally の相対統計誤差だけを取る。規格化の分の誤差は入れない。
4 本の平均の統計誤差は、4 本を独立とし、世代間の相関も無視して合わせたもので、実際の統計誤差の下限の見積もりになる。
xlc の分布は誤差を持たないとする。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc
from openmc_model import N, pin_rows

SIDES = ((0, -1), (-1, 0), (0, 1), (1, 0))


def read_pins(path: str) -> np.ndarray:
    rows = np.asarray(msgpack.unpackb(Path(path).read_bytes())["rows"])
    if rows.shape != (N, N):
        raise SystemExit(f"{path}: expected a {N} x {N} map")
    return rows / rows[rows > 0].mean()


def relative_sigma(work: str) -> np.ndarray:
    statepoints = sorted(
        Path(work).glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1])
    )
    if not statepoints:
        raise SystemExit(f"{work}: no statepoint")
    sp = openmc.StatePoint(str(statepoints[-1]))
    mean, std = pin_rows(sp), pin_rows(sp, "std_dev")
    return np.divide(std, mean, out=np.zeros_like(mean), where=mean > 0)


def summary(name: str, ratios: list[float], sigmas: list[float]) -> None:
    each = ", ".join(f"{100 * (r - 1):+.3f}" for r in ratios)
    mean = np.mean(ratios) - 1
    sigma = np.sqrt(np.sum(np.square(sigmas))) / len(sigmas)
    print(
        f"{name}: {each} %; mean {100 * mean:+.3f} %, "
        f"sigma assuming independent pins and generations {100 * sigma:.3f} %"
    )


def main() -> None:
    if len(sys.argv) not in (6, 8):
        raise SystemExit(__doc__)
    row, column = int(sys.argv[1]), int(sys.argv[2])
    xlc = read_pins(sys.argv[3])
    pins, sigma = read_pins(sys.argv[4]), relative_sigma(sys.argv[5])
    neighbours = [(row + di, column + dj) for di, dj in SIDES]
    for i, j in neighbours:
        if not (0 <= i < N and 0 <= j < N) or pins[i, j] <= 0 or xlc[i, j] <= 0:
            raise SystemExit(f"neighbour {(i, j)} is not a fuel pin in both maps")
    print("neighbours " + " ".join(str(p) for p in neighbours))
    summary(
        "xlc / OpenMC - 1",
        [xlc[p] / pins[p] for p in neighbours],
        [sigma[p] for p in neighbours],
    )
    if len(sys.argv) == 8:
        base, base_sigma = read_pins(sys.argv[6]), relative_sigma(sys.argv[7])
        summary(
            "OpenMC / base - 1",
            [pins[p] / base[p] for p in neighbours],
            [np.hypot(sigma[p], base_sigma[p]) for p in neighbours],
        )


if __name__ == "__main__":
    main()
