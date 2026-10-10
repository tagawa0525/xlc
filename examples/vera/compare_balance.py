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

最後に、燃料の群の間の散乱による移行を、4 eV〜15.03 keV に入る群（範囲の境の引数によらない）で比べる。群 g から
ほかの群へ出る率を O_g、同じ群に残る率を W_g、中性子束を φ_g（どれも生成 1 あたり）として、群ごとに
(O_g/φ_g)_xlc / (O_g/φ_g)_OpenMC − 1（W も同じ）と、出る率の不足の計 D = Σ_g [(O_g/φ_g)_xlc φ_g,OpenMC − O_g,OpenMC]
（xlc の移行の断面積を OpenMC の中性子束に掛けた率と OpenMC の率の差、生成 1 あたり）と D / Σ_g O_g,OpenMC を出す。
OpenMC の値は世代ごとの値の平均で、σ は世代ごとの値から求める（比は 1 次の近似で、中性子束と率の σ と相関を含む）。
"""

import math
import sys
from pathlib import Path

import msgpack

# 燃料の移行を比べる範囲 (eV)
RESONANCE = (4.0, 15.03e3)


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
    # 反応率 / 中性子束を作るので、中性子束が 0 の群があれば止める
    for label, side in [
        ("xlc", ours["materials"]),
        ("OpenMC", reference["materials"]),
    ] + [
        (f"OpenMC batch {i + 1}", b["materials"])
        for i, b in enumerate(reference["batches"])
    ]:
        for name, rates in side.items():
            zero = [g + 1 for g, f in enumerate(rates["flux"]) if not f > 0]
            if zero:
                raise SystemExit(f"{label}: {name} has no flux in groups {zero}")

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
                    phi = b["flux"][g]
                    da += a["absorption"][g] / a["flux"][g] * phi - b["absorption"][g]
                    dp += a["production"][g] / a["flux"][g] * phi - b["production"][g]
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
                f"{(fa / fb - 1) * 100:+7.2f}%"
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

    # 燃料の散乱による移行を 4 eV〜15.03 keV の群で比べる
    low, high = RESONANCE
    resonance = [
        g
        for g in range(groups)
        if bounds[g + 1] >= low * (1 - 1e-9) and bounds[g] <= high * (1 + 1e-9)
    ]
    fuel = ours["materials"]["fuel"]
    batches = [b["materials"]["fuel"] for b in reference["batches"]]
    if "transfer" not in fuel or any("out_scatter" not in b for b in batches):
        raise SystemExit("the fuel transfer needs `transfer` and batch `out_scatter`")

    def mean_of(key: str, g: int) -> float:
        return sum(b[key][g] for b in batches) / n

    def spread(values: list[float]) -> float:
        m = sum(values) / n
        return math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)

    def ratio(xs: float, key: str, g: int) -> tuple[float, float]:
        """xs / (key / φ)_OpenMC − 1 と、その σ（1 次の近似）"""
        rate, phi = mean_of(key, g), mean_of("flux", g)
        r = xs * phi / rate
        linear = [b["flux"][g] / phi - b[key][g] / rate for b in batches]
        return r - 1, r * spread(linear)

    print(
        f"fuel scattering transfer in {low:g}〜{high:g} eV "
        "(xlc / OpenMC − 1 of rate per flux; σ from batches):"
    )
    deficit = [0.0] * n
    out_total = [0.0] * n
    for g in resonance:
        out_xs = (
            sum(fuel["transfer"][g][h] for h in range(groups) if h != g)
            / fuel["flux"][g]
        )
        within_xs = fuel["transfer"][g][g] / fuel["flux"][g]
        out, out_s = ratio(out_xs, "out_scatter", g)
        within, within_s = ratio(within_xs, "within_scatter", g)
        print(
            f"  {g + 1:3d} {bounds[g + 1]:9.4g}〜{bounds[g]:<9.4g} "
            f"out {out:+7.2%} ± {out_s:.2%}  within {within:+7.2%} ± {within_s:.2%}"
        )
        for i, b in enumerate(batches):
            deficit[i] += out_xs * b["flux"][g] - b["out_scatter"][g]
            out_total[i] += b["out_scatter"][g]
    d, o = sum(deficit) / n, sum(out_total) / n
    share = [x / o - d * y / o**2 for x, y in zip(deficit, out_total)]
    print(
        f"  groups {resonance[0] + 1}〜{resonance[-1] + 1}: out-scatter deficit "
        f"D {d:+.3e} ± {spread(deficit):.1e} per unit production, "
        f"D / OpenMC out {d / o:+.2%} ± {spread(share):.2%}"
    )


if __name__ == "__main__":
    main()
