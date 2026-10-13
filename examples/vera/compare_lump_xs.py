"""xlc の構造材の塊の、核種ごとの実効断面積を OpenMC と比べ、正味の吸収の差を k の差に換算する（吸収体の自己遮蔽の検証）。

使い方: python3 compare_lump_xs.py <xlc の lump-xs.msgpack> <OpenMC の balance.msgpack> <格子の入力 (TOML)>
        <物質> <xlc のライブラリ> [核種 ...]

xlc 側は crates/lattice/examples/lump_xs.rs、OpenMC 側は openmc_balance.py の `nuclides`（物質の核種ごとの群ごとの率）と
`materials`（物質の中性子束）で作る。OpenMC の実効断面積は、率 / (個数密度 × 物質の中性子束)。xlc は塊ごとに Dancoff
係数が違うので、塊の領域の数で重みづけした平均を使う（塊の体積はどれも同じなので、体積の平均。中性子束の重みでは
ない）。個数密度は格子の入力の物質のもの。

核種ごとに、正味の吸収の実効断面積の比 (xlc / OpenMC − 1) を群ごとに出し、OpenMC の中性子束に固定した 1 次の近似で、
その差による k の差 Δk ≈ −k² N (σa,xlc − σa,OpenMC) φ_OpenMC（× 1e5。k は OpenMC のもの）を群ごとと範囲ごとに出す。
群は、OpenMC の吸収の率が物質の吸収の 0.5% 以上のものだけを出す。核種を省略すると、物質のすべての核種。σ は OpenMC の
実効断面積の相対標準偏差の上限 √(率² + 中性子束²)（中性子束の σ は世代ごとの値から）。
"""

import math
import sys
from pathlib import Path

import msgpack
import tomllib

# 範囲の境 (eV)
EDGES = (15.03e3, 4.0, 0.625)


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())


def main() -> None:
    if len(sys.argv) < 6:
        raise SystemExit(__doc__)
    ours, reference = read(sys.argv[1]), read(sys.argv[2])
    lattice = tomllib.loads(Path(sys.argv[3]).read_text())
    name = sys.argv[4]
    bounds = read(sys.argv[5])["group_bounds"]
    densities = {
        n.replace("-", ""): d for n, d in lattice["material"][name]["nuclides"].items()
    }
    nuclides = sys.argv[6:] or sorted(densities)
    groups = len(bounds) - 1
    k = reference["k"]
    flux = reference["materials"][name]["flux"]
    batches = [b["materials"][name]["flux"] for b in reference["batches"]]
    flux_std = [
        math.sqrt(
            sum((b[g] - flux[g]) ** 2 for b in batches)
            / (len(batches) - 1)
            / len(batches)
        )
        for g in range(groups)
    ]
    absorbed = reference["materials"][name]["absorption"]
    total_absorbed = sum(absorbed)
    lumps = ours["lumps"]
    weight = sum(lump["regions"] for lump in lumps)

    def ours_xs(nuclide: str, score: str) -> list[float]:
        return [
            sum(lump["regions"] * lump["nuclides"][nuclide][score][g] for lump in lumps)
            / weight
            for g in range(groups)
        ]

    def range_of(g: int) -> int:
        return sum(1 for e in EDGES if bounds[g] <= e)

    labels = (
        [f"> {EDGES[0]:g} eV"]
        + [f"{EDGES[i + 1]:g}〜{EDGES[i]:g} eV" for i in range(len(EDGES) - 1)]
        + [f"< {EDGES[-1]:g} eV"]
    )
    print(
        f"{name}: net absorption {total_absorbed:.5f} per production, k {k:.6f}; "
        f"{len(lumps)} lumps"
    )
    total_by_range = [0.0] * len(labels)
    for nuclide in nuclides:
        n = densities[nuclide]
        rates = reference["nuclides"][name][nuclide]
        theirs = [r / (n * f) for r, f in zip(rates["absorption"], flux)]
        sigma = [
            math.sqrt((s / r) ** 2 + (fs / f) ** 2) if r > 0 else 0.0
            for r, s, f, fs in zip(
                rates["absorption"], rates["absorption_std"], flux, flux_std
            )
        ]
        mine = ours_xs(nuclide, "absorption")
        dk = [-(k**2) * n * (a - b) * f * 1e5 for a, b, f in zip(mine, theirs, flux)]
        by_range = [0.0] * len(labels)
        for g in range(groups):
            by_range[range_of(g)] += dk[g]
        for i, v in enumerate(by_range):
            total_by_range[i] += v
        share = sum(rates["absorption"]) / total_absorbed
        print(
            f"\n{nuclide}: {share * 100:.1f}% of the absorption, Δk "
            + ", ".join(f"{label} {v:+.1f}" for label, v in zip(labels, by_range))
            + f", sum {sum(by_range):+.1f}"
        )
        for g in range(groups):
            if rates["absorption"][g] < 0.005 * total_absorbed:
                continue
            print(
                f"  group {g + 1:2d} ({bounds[g + 1]:.4g}〜{bounds[g]:.4g} eV): "
                f"σa xlc {mine[g]:.5g} OpenMC {theirs[g]:.5g} b, "
                f"{(mine[g] / theirs[g] - 1) * 100:+.2f}% ± {sigma[g] * 100:.2f}%, "
                f"Δk {dk[g]:+.1f}"
            )
    print(
        "\nall listed nuclides: "
        + ", ".join(f"{label} {v:+.1f}" for label, v in zip(labels, total_by_range))
        + f", sum {sum(total_by_range):+.1f}"
    )


if __name__ == "__main__":
    main()
