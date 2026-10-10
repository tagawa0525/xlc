"""xlc と OpenMC の balance.msgpack を、物質（または核種）ごと・群ごとに比べる。

使い方: python3 compare_groups.py <xlc の balance.msgpack> <OpenMC の balance.msgpack> <xlc のライブラリ>

群ごとに、正味の吸収の差の k の差への寄与（−k_xlc k_OpenMC ΔA_g × 1e5。`compare_balance.py` の分け方を群ごとに
したもの）と、xlc / OpenMC − 1 の中性子束、正味の吸収・生成・ほかの群へ出る移行・同じ群に残る移行の、反応率 /
中性子束の比を出す。OpenMC の値は世代ごとの値の平均で、σ は世代ごとの値から求める（比は 1 次の近似で、中性子束と
率の σ と相関を含む）。OpenMC の balance は世代ごとの `out_scatter` と `within_scatter` を持つこと
（`openmc_homogeneous.py`）。
"""

import math
import sys
from pathlib import Path

import msgpack

KEYS = ("absorption", "production", "out_scatter", "within_scatter")


def mean(batches: list[dict], key: str, g: int) -> float:
    return sum(t[key][g] for t in batches) / len(batches)


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    ours, reference = read(sys.argv[1]), read(sys.argv[2])
    bounds = read(sys.argv[3])["group_bounds"]
    groups = len(bounds) - 1
    names = sorted(ours["materials"])
    if names != sorted(reference["materials"]):
        raise SystemExit("the materials differ")
    batches = reference["batches"]
    n = len(batches)
    for b in batches:
        for name in names:
            if any(not f > 0 for f in b["materials"][name]["flux"]):
                raise SystemExit(f"{name}: a batch has no flux in some group")

    def spread(values: list[float]) -> float:
        m = sum(values) / n
        return math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)

    def k_of(materials: dict, below: float, above: float) -> float:
        return 1 / (
            sum(sum(r["absorption"]) for r in materials.values()) + below + above
        )

    k_ours = k_of(ours["materials"], 0.0, 0.0)
    k_refs = [k_of(b["materials"], b["below"], b["above"]) for b in batches]
    k_ref = sum(k_refs) / n
    print(
        f"k xlc {ours['k']:.6f} (1/A {k_ours:.6f}); OpenMC {reference['k']:.6f} ± "
        f"{reference['k_std']:.6f} (1/A {k_ref:.6f} ± {spread(k_refs):.6f})"
    )

    for name in names:
        mine = ours["materials"][name]
        theirs = [b["materials"][name] for b in batches]
        moved = mine["transfer"]
        xs = {
            "absorption": [
                mine["absorption"][g] / mine["flux"][g] for g in range(groups)
            ],
            "production": [
                mine["production"][g] / mine["flux"][g] for g in range(groups)
            ],
            "out_scatter": [
                (sum(moved[g]) - moved[g][g]) / mine["flux"][g] for g in range(groups)
            ],
            "within_scatter": [moved[g][g] / mine["flux"][g] for g in range(groups)],
        }
        print(
            f"{name}: group, energy (eV), Δk (pcm), xlc / OpenMC − 1 of flux and rate per flux "
            "(absorption, production, out-scatter, within-group)"
        )

        for g in range(groups):
            phi = mean(theirs, "flux", g)
            to_k = [
                -k_ours * k * 1e5 * (mine["absorption"][g] - t["absorption"][g])
                for k, t in zip(k_refs, theirs)
            ]
            flux = mine["flux"][g] / phi - 1
            flux_s = spread([t["flux"][g] / phi for t in theirs]) * (flux + 1)
            parts = [
                f"{g + 1:3d} {bounds[g + 1]:9.4g}〜{bounds[g]:<9.4g}",
                f"{sum(to_k) / n:+7.2f} ± {spread(to_k):4.2f}",
                f"φ {flux:+7.3%} ± {flux_s:.3%}",
            ]
            for key in KEYS:
                rate = mean(theirs, key, g)
                if rate == 0.0:
                    parts.append(f"{key[:3]} {'—':>17}")
                    continue
                r = xs[key][g] * phi / rate
                linear = [t["flux"][g] / phi - t[key][g] / rate for t in theirs]
                parts.append(f"{key[:3]} {r - 1:+7.3%} ± {r * spread(linear):.3%}")
            print("  " + "  ".join(parts))


if __name__ == "__main__":
    main()
