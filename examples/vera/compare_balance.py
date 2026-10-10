"""xlc と OpenMC の balance.msgpack を比べ、k の差を物質とエネルギーの範囲ごとの正味の吸収の差に分ける。

使い方: python3 compare_balance.py <xlc の balance.msgpack> <OpenMC の balance.msgpack> <xlc のライブラリ>
        [範囲の境 (eV) ...]

どちらも生成の和を 1 に規格化してあり、外周が鏡面反射なので正味の吸収の和が 1/k になる。そこで
1/k_xlc − 1/k_OpenMC = Σ（正味の吸収の差）が厳密に成り立ち、k の差 (pcm) は −k_OpenMC × k_xlc × ΔA × 1e5 で物質と範囲に
分けられる。生成の差（どこで中性子が生まれるか）と中性子束の差も同じ範囲で表示する。範囲の境の既定は 15.03 keV
（UFG の範囲の上）、4 eV（下）、0.625 eV。
"""

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
    edges = [float(e) for e in sys.argv[4:]] or [15.03e3, 4.0, 0.625]
    edges = sorted(edges, reverse=True)
    if set(ours["materials"]) != set(reference["materials"]):
        raise SystemExit(
            f"materials differ: {sorted(ours['materials'])} and {sorted(reference['materials'])}"
        )
    groups = len(bounds) - 1
    for side in (ours, reference):
        for rates in side["materials"].values():
            if any(len(v) != groups for v in rates.values()):
                raise SystemExit("the group counts differ from the library")

    # 群 g は bounds[g]〜bounds[g+1]（高いエネルギーから）。範囲は群の上の境で決める
    def range_of(g: int) -> int:
        return sum(1 for e in edges if bounds[g] <= e * (1 + 1e-9))

    labels = (
        [f"> {edges[0]:g} eV"]
        + [f"{edges[i + 1]:g}〜{edges[i]:g} eV" for i in range(len(edges) - 1)]
        + [f"< {edges[-1]:g} eV"]
    )
    absorbed = lambda side: sum(
        sum(r["absorption"]) for r in side["materials"].values()
    )
    # 和が厳密に閉じるよう、どちらも 1/（正味の吸収の和）を k とする（OpenMC の k の推定値とは統計誤差の程度で違う）
    k_ours, k_ref = 1 / absorbed(ours), 1 / absorbed(reference)
    print(
        f"k xlc {ours['k']:.6f} (1/A {k_ours:.6f}), OpenMC {reference['k']:.6f} (1/A {k_ref:.6f}), "
        f"diff of 1/A {(k_ours - k_ref) / k_ref * 1e5:+.1f} pcm"
    )
    scale = -k_ours * 1e5
    header = f"{'material':10s} {'range':>18s} {'Δk (pcm)':>9s} {'ΔA/A':>8s} {'ΔP':>9s} {'Δφ/φ':>8s}"
    print(header)
    total = 0.0
    for name in sorted(ours["materials"]):
        a, b = ours["materials"][name], reference["materials"][name]
        for r, label in enumerate(labels):
            sel = [g for g in range(groups) if range_of(g) == r]
            da = sum(a["absorption"][g] - b["absorption"][g] for g in sel)
            ab = sum(b["absorption"][g] for g in sel)
            dp = sum(a["production"][g] - b["production"][g] for g in sel)
            fa = sum(a["flux"][g] for g in sel)
            fb = sum(b["flux"][g] for g in sel)
            pcm = scale * da
            total += pcm
            print(
                f"{name:10s} {label:>18s} {pcm:+9.1f} "
                f"{(da / ab if ab else 0) * 100:+7.2f}% {dp:+9.2e} {(fa / fb - 1) * 100 if fb else 0:+7.2f}%"
            )
    print(f"{'sum':10s} {'':>18s} {total:+9.1f}")
    # 1 次の近似: 群ごとの実効断面積（反応率 / 中性子束）の差を OpenMC の中性子束に掛け、断面積の差による分を求める。
    # Δk/k ≈ Σ ΔνΣf φ / P − Σ ΔΣa φ / A（OpenMC の P = 1、A = 1/k）。全体の差との残りは、中性子束の形の差
    # （散乱による移行、核分裂スペクトル、輸送）による
    print("first order by effective cross sections on the OpenMC flux (pcm):")
    a_ref = 1 / k_ref
    first_total = 0.0
    for name in sorted(ours["materials"]):
        a, b = ours["materials"][name], reference["materials"][name]
        parts = []
        for r, label in enumerate(labels):
            sel = [g for g in range(groups) if range_of(g) == r]
            dabs = dprod = 0.0
            for g in sel:
                if a["flux"][g] > 0 and b["flux"][g] > 0:
                    phi = b["flux"][g]
                    dabs += (
                        a["absorption"][g] / a["flux"][g] - b["absorption"][g] / phi
                    ) * phi
                    dprod += (
                        a["production"][g] / a["flux"][g] - b["production"][g] / phi
                    ) * phi
            pcm = (dprod - dabs / a_ref) * 1e5
            first_total += pcm
            parts.append(
                f"{label}: {pcm:+7.1f} (P {dprod * 1e5:+7.1f}, A {-dabs / a_ref * 1e5:+7.1f})"
            )
        print(f"  {name:10s} " + " | ".join(parts))
    print(
        f"  first order sum {first_total:+.1f}, the rest (flux shape) {total - first_total:+.1f}"
    )
    # 燃料の核種ごと（xlc は群ごとの燃料の値を核種の断面積の比で分けた近似）
    shared = sorted(set(ours.get("nuclides", {})) & set(reference.get("nuclides", {})))
    if shared:
        print("fuel nuclides (xlc split by shielded nuclide shares):")
        for key in shared:
            a, b = ours["nuclides"][key], reference["nuclides"][key]
            parts = []
            for r, label in enumerate(labels):
                sel = [g for g in range(groups) if range_of(g) == r]
                da = sum(a["absorption"][g] - b["absorption"][g] for g in sel)
                dp = sum(a["production"][g] - b["production"][g] for g in sel)
                parts.append(f"{label}: A {scale * da:+7.1f} P {dp:+.2e}")
            print(f"  {key:14s} " + " | ".join(parts))
    # 群ごとの寄与の大きいもの
    rows = []
    for name in sorted(ours["materials"]):
        a, b = ours["materials"][name], reference["materials"][name]
        for g in range(groups):
            da = a["absorption"][g] - b["absorption"][g]
            ab = b["absorption"][g]
            fb = b["flux"][g]
            rows.append(
                (
                    scale * da,
                    name,
                    g + 1,
                    bounds[g + 1],
                    bounds[g],
                    da / ab if ab else 0.0,
                    a["flux"][g] / fb - 1 if fb else 0.0,
                )
            )
    print("largest group contributions (group from 1, energies in eV):")
    for pcm, name, g, low, high, rel, flux in sorted(rows, key=lambda r: -abs(r[0]))[
        :20
    ]:
        print(
            f"{name:10s} {g:3d} {low:10.4g}〜{high:<10.4g} {pcm:+8.1f} pcm  ΔA/A {rel * 100:+6.2f}%  Δφ/φ {flux * 100:+6.2f}%"
        )


if __name__ == "__main__":
    main()
