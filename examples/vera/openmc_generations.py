"""OpenMC の固有値計算の、世代ごとの k から統計誤差と計算どうしの相関を見積もり、ピン出力の統計誤差を表示する。

使い方: python3 openmc_generations.py <基準の作業ディレクトリ> <作業ディレクトリ> ...

作業ディレクトリは openmc_model.py のもの（statepoint の番号が最も大きいものを読む）。計算ごとに、表示の k（3 つの推定量の
組み合わせ）と統計誤差、世代ごとの k（衝突の推定量）の平均の統計誤差（世代ごと、5〜50 世代でまとめたもの）、隣の世代との
相関係数、ピンの核分裂率（名前が fission の tally）の相対統計誤差の中央値と最大を表示する。基準でない計算ごとに、基準との
世代ごとの k の相関係数（世代ごと、10 と 20 世代の平均）と、差の平均の統計誤差（世代ごと、10〜30 世代でまとめたもの）を
表示する。まとめる世代の数を変えて統計誤差が大きくなれば、世代間の相関で表示の統計誤差が小さめに出ている。
"""

import sys
from pathlib import Path

import numpy as np
import openmc

BLOCKS = (1, 5, 10, 20, 30, 50)
CORRELATION_BLOCKS = (1, 10, 20)
DIFFERENCE_BLOCKS = (1, 10, 20, 30)


def block_means(values: np.ndarray, size: int) -> np.ndarray:
    n = len(values) // size * size
    return values[:n].reshape(-1, size).mean(axis=1)


def sigma_of_mean(values: np.ndarray, size: int) -> float:
    means = block_means(values, size)
    return float(means.std(ddof=1) / np.sqrt(len(means)))


def read(work: Path) -> tuple[openmc.StatePoint, np.ndarray]:
    statepoints = sorted(
        work.glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1])
    )
    if not statepoints:
        raise SystemExit(f"{work}: no statepoint")
    sp = openmc.StatePoint(str(statepoints[-1]))
    return sp, np.asarray(sp.k_generation)[sp.n_inactive :]


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    works = [Path(a) for a in sys.argv[1:]]
    runs = {}
    for work in works:
        sp, k = read(work)
        runs[work] = k
        tally = sp.get_tally(name="fission")
        mean, std = tally.mean.ravel(), tally.std_dev.ravel()
        relative = std[mean > 0] / mean[mean > 0]
        sigmas = " ".join(f"{b}:{1e5 * sigma_of_mean(k, b):.1f}" for b in BLOCKS)
        lag = np.corrcoef(k[:-1], k[1:])[0, 1]
        print(
            f"{work.name}: k {sp.keff.nominal_value:.6f} +- {1e5 * sp.keff.std_dev:.1f} pcm; "
            f"generations {len(k)}, sigma (pcm, by block size) {sigmas}; lag-1 {lag:+.3f}; "
            f"pin sigma median {100 * np.median(relative):.3f}% max {100 * relative.max():.3f}%"
        )
    base = runs[works[0]]
    for work in works[1:]:
        k = runs[work]
        if len(k) != len(base):
            raise SystemExit(f"{work}: {len(k)} generations against {len(base)}")
        correlations = " ".join(
            f"{b}:{np.corrcoef(block_means(base, b), block_means(k, b))[0, 1]:+.3f}"
            for b in CORRELATION_BLOCKS
        )
        sigmas = " ".join(
            f"{b}:{1e5 * sigma_of_mean(k - base, b):.1f}" for b in DIFFERENCE_BLOCKS
        )
        print(
            f"{work.name} - {works[0].name}: correlation (by block size) {correlations}; "
            f"sigma of the difference (pcm, by block size) {sigmas}"
        )


if __name__ == "__main__":
    main()
