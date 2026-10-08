"""xlc と OpenMC の燃料の核種ごとの実効捕獲・核分裂断面積を群ごとに比べ、捕獲率の差を k の差に換算する。

使い方: python3 compare_capture.py <xlc のライブラリ> <xlc の fuel-xs.msgpack> <OpenMC の fuel-xs.msgpack>
        [核種 ...]

xlc 側は crates/lattice/examples/fuel_xs.rs、OpenMC 側は openmc_model.py で作る。捕獲は全 − 弾性 − 核分裂
（UFG の範囲より下では非弾性の反応は捕獲と核分裂だけ）。群ごとの相対差 (xlc / OpenMC − 1) を、
UFG の範囲（WIMS 69 群の 15〜27 群）と、その下の 4〜0.1 eV（28〜56 群）について表示する。
捕獲率の差は OpenMC の燃料の中性子束で Σ_g φ_g N (σ_xlc − σ_OpenMC) と足し、セル全体の吸収率に対する比
（中性子束を固定した 1 次の Δk/k ≈ −ΔC/A）で表す。核種は OpenMC の名前で、省略すると
Pu240, Pu239, Pu241, Pu242, U238。
"""

import sys
from pathlib import Path

import msgpack

UFG = range(15, 28)
BELOW = range(28, 57)


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())


def reactions(rates: dict) -> dict[str, list[float]]:
    capture = [
        t - e - f for t, e, f in zip(rates["total"], rates["elastic"], rates["fission"])
    ]
    return {"capture": capture, "fission": rates["fission"]}


def main() -> None:
    bounds = read(sys.argv[1])["group_bounds"]
    ours = read(sys.argv[2])["nuclides"]
    reference = read(sys.argv[3])
    nuclides = sys.argv[4:] or ["Pu240", "Pu239", "Pu241", "Pu242", "U238"]
    phi = reference["flux"]
    absorption = reference["absorption"]
    k = reference["keff"]
    print(f"OpenMC k {k[0]:.6f} +- {k[1]:.6f}")

    for score in ["capture", "fission"]:
        print(f"\n{score}: xlc / OpenMC - 1 (%), OpenMC rate rsd (%) in brackets")
        print("group  E_low (eV)  " + "  ".join(f"{n:>15}" for n in nuclides))
        for g in [*UFG, *BELOW]:
            cells = []
            for n in nuclides:
                a = reactions(ours[n])[score][g - 1]
                b = reactions(reference["nuclides"][n])[score][g - 1]
                rsd = reference["nuclides"][n][f"{score}-rsd"][g - 1]
                cells.append(
                    f"{(a / b - 1) * 100:+7.2f} [{rsd * 100:4.2f}]"
                    if b > 0
                    else f"{'-':>15}"
                )
            print(f"{g:5d}  {bounds[g]:10.4g}  " + "  ".join(cells))

    print("\ncapture/fission rate difference (sum over groups, OpenMC flux)")
    print("nuclide  reaction  group range  rel. diff (%)  dR / A (pcm)")
    density = reference["densities"]
    for n in nuclides:
        for score in ["capture", "fission"]:
            a = reactions(ours[n])[score]
            b = reactions(reference["nuclides"][n])[score]
            for label, groups in [
                ("all", range(1, len(phi) + 1)),
                ("15-27", UFG),
                ("28-56", BELOW),
            ]:
                ref = sum(phi[g - 1] * density[n] * b[g - 1] for g in groups)
                diff = sum(
                    phi[g - 1] * density[n] * (a[g - 1] - b[g - 1]) for g in groups
                )
                if ref <= 0:
                    continue
                print(
                    f"{n:7s}  {score:8s}  {label:11s}  {diff / ref * 100:+13.2f}"
                    f"  {diff / absorption * 1e5:+12.1f}"
                )


if __name__ == "__main__":
    main()
