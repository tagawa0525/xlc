"""xlc と OpenMC の balance.msgpack を比べ、k の差を物質とエネルギーの範囲ごとの正味の吸収の差に分ける。

使い方: python3 compare_balance.py <xlc の balance.msgpack> <OpenMC の balance.msgpack> <xlc のライブラリ>
        [範囲の境 (eV) ...]

どちらも生成の和を 1 に規格化してあり（OpenMC はすべてのエネルギーの和）、外周が鏡面反射なので正味の吸収の和が
1/k になる。そこで k_xlc − k_OpenMC = −k_xlc k_OpenMC Σ ΔA が厳密に成り立ち、k の差（× 1e5、検証の記録の pcm と同じ
絶対の差）を物質と範囲に分けられる。OpenMC の群の範囲の外（いちばん下の境より下と、いちばん上の境より上）の吸収は
xlc にないので、別の行にする。

中性子束は、xlc も物理的な値（SPH 因子を戻したもの）なので、反応率 / 中性子束は物理的な実効断面積になる。1 次の近似で
は、群ごとの吸収と生成の実効断面積の差を OpenMC の中性子束に掛けて、それらの断面積の差による分を求める（Δk ≈
k (ΔP − k ΔA)、中性子束を OpenMC に固定）。全体の差との残りは、散乱による移行の差による分のほか、範囲の外の項、2 次の
項、断面積の差が下流の中性子束に及ぼす間接の影響を含むので、原因ごとの分け方ではない。

σ は OpenMC の世代ごとの値の標準誤差（同じ粒子の追跡による量の間の相関、規格化の生成、範囲の外の項を含む。世代の
間の相関は含まないので小さめに出る）。範囲の境の既定は 15.03 keV（UFG の範囲の上）、4 eV（下）、0.625 eV。

最後に、燃料の群の間の散乱による移行（反応率 / 中性子束）を、共鳴の範囲の群で比べる（群から出る割合と、群に残る割合）。
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
    keys = ("flux", "absorption", "production")
    for side in (ours, reference):
        for rates in side["materials"].values():
            if any(len(rates[key]) != groups for key in keys):
                raise SystemExit("the group counts differ from the library")
    if not reference.get("batches") or "outside" not in reference:
        raise SystemExit("the OpenMC balance needs `batches` and `outside`")

    # 群 g は bounds[g]〜bounds[g+1]（高いエネルギーから）。範囲は群の上の境で決める
    labels = (
        [f"> {edges[0]:g} eV"]
        + [f"{edges[i + 1]:g}〜{edges[i]:g} eV" for i in range(len(edges) - 1)]
        + [f"< {edges[-1]:g} eV"]
    )
    ranges = [[] for _ in labels]
    for g in range(groups):
        ranges[sum(1 for e in edges if bounds[g] <= e * (1 + 1e-9))].append(g)

    k_ours = 1 / sum(sum(r["absorption"]) for r in ours["materials"].values())

    def split(ref: dict, below: float, above: float) -> dict:
        """OpenMC の 1 組の値（平均か 1 世代）に対する、物質と範囲ごとの寄与と 1 次の近似の寄与。"""
        a_ref = sum(sum(r["absorption"]) for r in ref.values()) + below + above
        k_ref = 1 / a_ref
        to_k = -k_ours * k_ref * 1e5
        parts = {"k_ref": k_ref, "below": to_k * -below, "above": to_k * -above}
        for name in names:
            a, b = ours["materials"][name], ref[name]
            for r, sel in enumerate(ranges):
                parts[(name, r)] = to_k * sum(
                    a["absorption"][g] - b["absorption"][g] for g in sel
                )
                dp = da = 0.0
                for g in sel:
                    if a["flux"][g] > 0 and b["flux"][g] > 0:
                        phi = b["flux"][g]
                        da += (
                            a["absorption"][g] / a["flux"][g] * phi - b["absorption"][g]
                        )
                        dp += (
                            a["production"][g] / a["flux"][g] * phi - b["production"][g]
                        )
                parts[("first", name, r)] = k_ref * (dp - k_ref * da) * 1e5
        return parts

    outside = reference["outside"]
    mean = split(
        reference["materials"],
        outside["below"]["absorption"],
        outside["above"]["absorption"],
    )
    samples = [
        split(b["materials"], b["below"], b["above"]) for b in reference["batches"]
    ]
    n = len(samples)

    def sigma(key) -> float:
        """`key` は寄与の名前か、1 組の寄与から値を求める関数"""
        values = [key(s) if callable(key) else s[key] for s in samples]
        m = sum(values) / n
        return math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)

    k_ref = mean["k_ref"]
    print(
        f"k xlc {ours['k']:.6f} (1/A {k_ours:.6f}); OpenMC {reference['k']:.6f} ± {reference['k_std']:.6f} "
        f"(1/A {k_ref:.6f} ± {sigma('k_ref'):.6f}); "
        f"k difference by 1/A {(k_ours - k_ref) * 1e5:+.1f} ± {sigma('k_ref') * 1e5:.1f}, "
        f"by the k estimate {(ours['k'] - reference['k']) * 1e5:+.1f} ± {reference['k_std'] * 1e5:.1f}"
    )
    print(
        f"outside the groups: production below {outside['below']['production']:.2e}, "
        f"above {outside['above']['production']:.2e}"
    )

    def total(s: dict) -> float:
        return (
            s["below"]
            + s["above"]
            + sum(s[(name, r)] for name in names for r in range(len(labels)))
        )

    def first(s: dict) -> float:
        return sum(s[("first", name, r)] for name in names for r in range(len(labels)))

    print(f"{'material':9s} {'range':>18s} {'Δk':>14s} {'Δφ/φ':>8s}")
    for name in names:
        a, b = ours["materials"][name], reference["materials"][name]
        for r, (label, sel) in enumerate(zip(labels, ranges)):
            fa = sum(a["flux"][g] for g in sel)
            fb = sum(b["flux"][g] for g in sel)
            print(
                f"{name:9s} {label:>18s} {mean[(name, r)]:+7.1f} ± {sigma((name, r)):4.1f} "
                f"{(fa / fb - 1) * 100 if fb else 0:+7.2f}%"
            )
    for side in ("below", "above"):
        print(f"{'outside':9s} {side:>18s} {mean[side]:+7.1f} ± {sigma(side):4.1f}")
    print(f"{'sum':9s} {'':>18s} {total(mean):+7.1f} ± {sigma(total):4.1f}")

    print("first order by absorption and production cross sections on the OpenMC flux:")
    for name in names:
        parts = [
            f"{label}: {mean[('first', name, r)]:+7.1f} ± {sigma(('first', name, r)):4.1f}"
            for r, label in enumerate(labels)
        ]
        print(f"  {name:9s} " + " | ".join(parts))
    rest = lambda s: total(s) - first(s)
    print(
        f"  first order sum {first(mean):+.1f} ± {sigma(first):.1f}, "
        f"the rest {rest(mean):+.1f} ± {sigma(rest):.1f}"
    )

    # 燃料の散乱による移行（反応率 / 中性子束）を共鳴の範囲の群で比べる
    if "transfer" in ours["materials"].get("fuel", {}) and "transfer" in reference[
        "materials"
    ].get("fuel", {}):
        a, b = ours["materials"]["fuel"], reference["materials"]["fuel"]
        print(
            "fuel scattering transfer in the resonance range (xlc / OpenMC − 1, OpenMC σ):"
        )
        for g in ranges[1]:
            out_a = (
                sum(a["transfer"][g][h] for h in range(groups) if h != g) / a["flux"][g]
            )
            out_b = (
                sum(b["transfer"][g][h] for h in range(groups) if h != g) / b["flux"][g]
            )
            out_s = (
                math.sqrt(
                    sum(b["transfer_std"][g][h] ** 2 for h in range(groups) if h != g)
                )
                / b["flux"][g]
            )
            in_a = a["transfer"][g][g] / a["flux"][g]
            in_b = b["transfer"][g][g] / b["flux"][g]
            in_s = b["transfer_std"][g][g] / b["flux"][g]
            print(
                f"  {g + 1:3d} {bounds[g + 1]:9.4g}〜{bounds[g]:<9.4g} "
                f"out {out_a / out_b - 1:+7.2%} ± {out_s / out_b:.2%}  "
                f"in {in_a / in_b - 1:+7.2%} ± {in_s / in_b:.2%}"
            )


if __name__ == "__main__":
    main()
