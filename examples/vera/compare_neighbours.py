"""セルの上下左右の隣の 4 本のピンについて、xlc と OpenMC の核分裂率の差と、OpenMC の計算どうしの変化を、
OpenMC の統計誤差つきで表示する（計装管のシンブルのような、1 つのセルに入れたものの影響を見るため）。

使い方: python3 compare_neighbours.py <行> <列> <xlc の pins.msgpack> <OpenMC の作業ディレクトリ>
        [<基準の OpenMC の作業ディレクトリ>]

行と列は 0 始まりで上の行から（lattice のピン出力と同じ向き）。作業ディレクトリは openmc_model.py のもの（statepoint の
番号が最も大きいものを読む）。どの分布も燃料ピン（核分裂率が正のセル）の平均で規格化する。統計誤差は名前が fission の
tally の標準偏差の相対値で、規格化の分の誤差は入れない。4 本の平均の統計誤差は、ピンごとの統計誤差を独立として合わせる。
xlc の分布はそのまま使い、誤差を持たないとする。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

N = 17
SIDES = ((0, -1), (-1, 0), (0, 1), (1, 0))


def openmc_pins(work: Path) -> tuple[np.ndarray, np.ndarray]:
    """燃料ピンの平均で 1 にした核分裂率と、その相対統計誤差（17 × 17、上の行から）。"""
    statepoints = sorted(
        work.glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1])
    )
    if not statepoints:
        raise SystemExit(f"{work}: no statepoint")
    tally = openmc.StatePoint(str(statepoints[-1])).get_tally(name="fission")
    # メッシュは x が速く変わり、y は下から上。上の行から並べる
    mean = tally.mean.ravel().reshape(N, N)[::-1, :]
    std = tally.std_dev.ravel().reshape(N, N)[::-1, :]
    relative = np.divide(std, mean, out=np.zeros_like(mean), where=mean > 0)
    return mean / mean[mean > 0].mean(), relative


def summary(name: str, ratios: list[float], sigmas: list[float]) -> None:
    each = ", ".join(f"{100 * (r - 1):+.3f}" for r in ratios)
    mean = np.mean(ratios) - 1
    sigma = np.sqrt(np.sum(np.square(sigmas))) / len(sigmas)
    print(f"{name}: {each} %; mean {100 * mean:+.3f} +- {100 * sigma:.3f} %")


def main() -> None:
    if len(sys.argv) not in (5, 6):
        raise SystemExit(__doc__)
    row, column = int(sys.argv[1]), int(sys.argv[2])
    xlc = np.asarray(msgpack.unpackb(Path(sys.argv[3]).read_bytes())["rows"])
    if xlc.shape != (N, N):
        raise SystemExit("expected a 17 x 17 map")
    xlc = xlc / xlc[xlc > 0].mean()
    pins, sigma = openmc_pins(Path(sys.argv[4]))
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
    if len(sys.argv) == 6:
        base, base_sigma = openmc_pins(Path(sys.argv[5]))
        summary(
            "OpenMC / base - 1",
            [pins[p] / base[p] for p in neighbours],
            [np.hypot(sigma[p], base_sigma[p]) for p in neighbours],
        )


if __name__ == "__main__":
    main()
