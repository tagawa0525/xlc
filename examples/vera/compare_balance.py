"""xlc と OpenMC の balance.msgpack を比べ、k の差を物質とエネルギーの範囲ごとの正味の吸収の差に分ける。

使い方: python3 compare_balance.py <xlc の balance.msgpack> <OpenMC の balance.msgpack> <xlc のライブラリ>
        [範囲の境 (eV) ...]

どちらも生成の和を 1 に規格化してあり（OpenMC はすべてのエネルギーの和）、外周が鏡面反射なので正味の吸収の和が
1/k になる。そこで k_xlc − k_OpenMC = −k_xlc k_OpenMC Σ ΔA が厳密に成り立ち、k の差（× 1e5、検証の記録の pcm と同じ
絶対の差）を物質と範囲に分けられる。OpenMC の群の範囲の外（10 MeV より上）の吸収は xlc にないので、別の行にする。

中性子束は、xlc も物理的な値（SPH 因子を掛けた領域と群では MOC の値に因子を掛けたもの）なので、反応率 / 中性子束は
SPH 因子を掛ける前の物理的な実効断面積になる。1 次の近似では、群ごとの実効断面積の差を OpenMC の中性子束に掛けて
断面積の差による分を求め（Δk ≈ k (ΔP − k ΔA)、中性子束を OpenMC に固定）、全体の差との残りを中性子束の形の差による分と
する。σ は OpenMC の統計誤差だけ（群の間の相関は無視する）。範囲の境の既定は 15.03 keV（UFG の範囲の上）、4 eV（下）、
0.625 eV。
"""

import math
import sys
from pathlib import Path

import msgpack


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    ours, reference = read(sys.argv[1]), read(sys.argv[2])
    bounds = msgpack.unpackb(Path(sys.argv[3]).read_bytes())["group_bounds"]
    edges = sorted((float(e) for e in sys.argv[4:]), reverse=True) or [
        15.03e3,
        4.0,
        0.625,
    ]
    names = sorted(ours["materials"])
    if names != sorted(reference["materials"]):
        raise SystemExit(
            f"materials differ: {names} and {sorted(reference['materials'])}"
        )
    groups = len(bounds) - 1
    for side in (ours, reference):
        for rates in side["materials"].values():
            if any(
                len(rates[key]) != groups
                for key in ("flux", "absorption", "production")
            ):
                raise SystemExit("the group counts differ from the library")
    if "outside" not in reference:
        raise SystemExit("the OpenMC balance has no `outside` part")

    # 群 g は bounds[g]〜bounds[g+1]（高いエネルギーから）。範囲は群の上の境で決める
    labels = (
        [f"> {edges[0]:g} eV"]
        + [f"{edges[i + 1]:g}〜{edges[i]:g} eV" for i in range(len(edges) - 1)]
        + [f"< {edges[-1]:g} eV"]
    )
    ranges = [[] for _ in labels]
    for g in range(groups):
        ranges[sum(1 for e in edges if bounds[g] <= e * (1 + 1e-9))].append(g)

    a_ours = sum(sum(r["absorption"]) for r in ours["materials"].values())
    outside = reference["outside"]["absorption"]
    a_ref = sum(sum(r["absorption"]) for r in reference["materials"].values()) + outside
    k_ours, k_ref = 1 / a_ours, 1 / a_ref
    to_k = -k_ours * k_ref * 1e5
    print(
        f"k xlc {ours['k']:.6f} (1/A {k_ours:.6f}); OpenMC {reference['k']:.6f} "
        f"(1/A {k_ref:.6f}, outside the groups: production {reference['outside']['production']:.2e}, "
        f"absorption {outside:.2e}); k difference (1/A) {(k_ours - k_ref) * 1e5:+.1f}"
    )

    def std(values) -> float:
        return math.sqrt(sum(v * v for v in values))

    print(f"{'material':9s} {'range':>18s} {'Δk':>14s} {'ΔA/A':>8s} {'Δφ/φ':>8s}")
    total = 0.0
    total_var = 0.0
    for name in names:
        a, b = ours["materials"][name], reference["materials"][name]
        for label, sel in zip(labels, ranges):
            da = sum(a["absorption"][g] - b["absorption"][g] for g in sel)
            sa = std(b["absorption_std"][g] for g in sel)
            ab = sum(b["absorption"][g] for g in sel)
            fa = sum(a["flux"][g] for g in sel)
            fb = sum(b["flux"][g] for g in sel)
            total += to_k * da
            total_var += (to_k * sa) ** 2
            print(
                f"{name:9s} {label:>18s} {to_k * da:+7.1f} ± {abs(to_k) * sa:4.1f} "
                f"{(da / ab if ab else 0) * 100:+7.2f}% {(fa / fb - 1) * 100 if fb else 0:+7.2f}%"
            )
    total += to_k * (0.0 - outside)
    print(f"{'outside':9s} {'> 10 MeV (OpenMC)':>18s} {to_k * -outside:+7.1f}")
    print(f"{'sum':9s} {'':>18s} {total:+7.1f} ± {math.sqrt(total_var):4.1f}")

    # 1 次の近似: 群ごとの物理的な実効断面積の差を、OpenMC の中性子束に掛ける
    print("first order by effective cross sections on the OpenMC flux:")
    first = 0.0
    for name in names:
        a, b = ours["materials"][name], reference["materials"][name]
        parts = []
        for label, sel in zip(labels, ranges):
            dp = da = var = 0.0
            for g in sel:
                if a["flux"][g] <= 0 or b["flux"][g] <= 0:
                    continue
                phi = b["flux"][g]
                sx_a = a["absorption"][g] / a["flux"][g]
                sx_p = a["production"][g] / a["flux"][g]
                da += sx_a * phi - b["absorption"][g]
                dp += sx_p * phi - b["production"][g]
                # OpenMC の反応率と、xlc の断面積 × OpenMC の中性子束の誤差
                var += (k_ref * b["absorption_std"][g]) ** 2 + b["production_std"][
                    g
                ] ** 2
                var += ((sx_p - k_ref * sx_a) * b["flux_std"][g]) ** 2
            dk = k_ref * (dp - k_ref * da) * 1e5
            first += dk
            parts.append(f"{label}: {dk:+7.1f} ± {k_ref * math.sqrt(var) * 1e5:4.1f}")
        print(f"  {name:9s} " + " | ".join(parts))
    print(
        f"  first order sum {first:+.1f}, the rest (flux shape and outside) {total - first:+.1f}"
    )

    # 群ごとの寄与の大きいもの（σ の 3 倍を超えるものだけ）
    rows = []
    for name in names:
        a, b = ours["materials"][name], reference["materials"][name]
        for g in range(groups):
            da = a["absorption"][g] - b["absorption"][g]
            sa = b["absorption_std"][g]
            if abs(da) <= 3 * sa:
                continue
            ab, fb = b["absorption"][g], b["flux"][g]
            rows.append(
                (
                    to_k * da,
                    abs(to_k) * sa,
                    name,
                    g + 1,
                    bounds[g + 1],
                    bounds[g],
                    da / ab if ab else 0.0,
                    a["flux"][g] / fb - 1 if fb else 0.0,
                    b["flux_std"][g] / fb if fb else 0.0,
                )
            )
    print("largest group contributions beyond 3σ (group from 1, energies in eV):")
    largest = sorted(rows, key=lambda r: -abs(r[0]))[:20]
    for pcm, sigma, name, g, low, high, rel, flux, flux_sigma in largest:
        print(
            f"{name:9s} {g:3d} {low:10.4g}〜{high:<10.4g} {pcm:+7.1f} ± {sigma:4.1f}  "
            f"ΔA/A {rel * 100:+6.2f}%  Δφ/φ {flux * 100:+6.2f}% ± {flux_sigma * 100:.2f}%"
        )


if __name__ == "__main__":
    main()
