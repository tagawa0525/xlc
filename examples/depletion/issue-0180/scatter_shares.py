"""燃料の核種ごとの、自由ガスの範囲の散乱源への寄与を、エネルギーの範囲ごとに出す（#0180）。

使い方: python3 scatter_shares.py <履歴> <燃焼度 (MWd/kg)> <核種の出力のディレクトリ> <温度 (K)> <核種>...

履歴（lattice の燃焼の msgpack）の燃焼度の点の、燃料のすべての燃焼の領域の個数密度を体積で平均し、核種の出力
（xs-nuclides/<核種>.msgpack）の温度のポイントワイズの弾性散乱断面積と全断面積から、核種ごとに

- 平均: 範囲の中の局所の割合 N_j σs,j / Σs（Σs は与えた核種の和）のレサジーの平均
- NR: 狭い共鳴の中性子束 1/(Σt E) で重みを付けた散乱の率の割合 ∫ N_j σs,j / Σt du / ∫ Σs / Σt du
- 最大: 範囲の中の局所の割合の最大

を出す。格子は与えた核種のポイントワイズの格子の和で、断面積は lin-lin に補間する（lattice の UFG と同じ）。割合は 1e−6 で切らず
そのまま出す。Σs と Σt は与えた核種だけの和なので、核種を落とすと分母が変わる（lattice の UFG のリングの核種をすべて与える）。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

# エネルギーの範囲 (eV)。31 eV は分離共鳴を持つ核種の自由ガスの上限（ADR-0243）、1000 eV は持たない核種の上限（ADR-0246）、
# 4 eV は UFG の範囲の下限
RANGES = [(4.0, 31.0), (31.0, 100.0), (100.0, 1000.0)]


def load(
    path: Path, temperature: float
) -> tuple[str, float, bool, np.ndarray, np.ndarray, np.ndarray]:
    """核種の出力の（名前、質量比、分離共鳴を持つか、エネルギー、弾性散乱、全断面積）"""
    d = msgpack.unpackb(path.read_bytes(), strict_map_key=False)
    (nuclide,) = d["nuclides"]
    tables = [t for t in nuclide["temperatures"] if t["temperature"] == temperature]
    if len(tables) != 1:
        raise SystemExit(f"{path}: no table at {temperature} K")
    p = tables[0]["pointwise"]
    return (
        nuclide["name"],
        nuclide["awr"],
        nuclide["resolved"] is not None,
        np.array(p["energy"]),
        np.array(p["elastic"]),
        np.array(p["total"]),
    )


def fuel_densities(history: dict, burnup: float) -> dict[str, float]:
    """燃焼度の点の、燃料（物質の名前が fuel で始まる）の燃焼の領域の個数密度の体積の平均"""
    points = [p for p in history["points"] if p["burnup"] == burnup]
    if len(points) != 1:
        raise SystemExit(f"no point at {burnup} MWd/kg")
    regions = [
        (i, r["volume"])
        for i, r in enumerate(history["regions"])
        if r["material"].startswith("fuel")
    ]
    volume = sum(v for _, v in regions)
    densities = points[0]["densities"]
    return {
        name: sum(densities[i][k] * v for i, v in regions) / volume
        for k, name in enumerate(history["nuclides"])
    }


def shares(
    grid: np.ndarray,
    scatter: dict[str, np.ndarray],
    total: np.ndarray,
    low: float,
    high: float,
) -> dict[str, tuple[float, float, float]]:
    """範囲 [low, high] の核種ごとの（平均、NR、最大）。`grid` は昇順のエネルギー、`scatter` は核種の N σs、`total` は Σt"""
    inside = (grid >= low) & (grid <= high)
    e = grid[inside]
    u = -np.log(e)
    sigma_s = sum(s[inside] for s in scatter.values())
    sigma_t = total[inside]

    def integral(y: np.ndarray) -> float:
        # レサジーは e の昇順で降順なので符号を戻す
        return -float(np.trapezoid(y, u))

    width = integral(np.ones_like(e))
    rate = integral(sigma_s / sigma_t)
    out = {}
    for name, s in scatter.items():
        local = s[inside] / sigma_s
        out[name] = (
            integral(local) / width,
            integral(s[inside] / sigma_t) / rate,
            float(local.max()),
        )
    return out


def main() -> None:
    if len(sys.argv) < 6:
        raise SystemExit(__doc__)
    history = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
    burnup = float(sys.argv[2])
    directory = Path(sys.argv[3])
    temperature = float(sys.argv[4])
    densities = fuel_densities(history, burnup)
    data = [load(directory / f"{name}.msgpack", temperature) for name in sys.argv[5:]]
    grid = np.unique(
        np.concatenate(
            [d[3][(d[3] >= RANGES[0][0]) & (d[3] <= RANGES[-1][1])] for d in data]
        )
    )
    grid = np.unique(np.concatenate([grid, [b for r in RANGES for b in r]]))
    scatter = {}
    total = np.zeros_like(grid)
    for name, _, _, energy, elastic, tot in data:
        n = densities[name]
        scatter[name] = n * np.interp(grid, energy, elastic)
        total += n * np.interp(grid, energy, tot)
    tables = [shares(grid, scatter, total, low, high) for low, high in RANGES]
    header = "核種 質量比 分離共鳴 個数密度 " + " ".join(
        f"{low:g}-{high:g}eV(平均,NR,最大)" for low, high in RANGES
    )
    print(header)
    for name, awr, resolved, *_ in sorted(data, key=lambda d: -tables[0][d[0]][1]):
        cells = " ".join(" ".join(f"{x:.3e}" for x in t[name]) for t in tables)
        print(
            f"{name} {awr:.4f} {'あり' if resolved else 'なし'} {densities[name]:.4e} {cells}"
        )


if __name__ == "__main__":
    main()
