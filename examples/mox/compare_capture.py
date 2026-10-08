"""xlc と OpenMC の燃料の核種ごとの実効捕獲・核分裂断面積を群ごとに比べ、捕獲率の差を k の差に換算する。

使い方: python3 compare_capture.py <xlc のライブラリ> <xlc の fuel-xs.msgpack> <OpenMC の fuel-xs.msgpack>
        [核種 ...]

xlc 側は crates/lattice/examples/fuel_xs.rs、OpenMC 側は openmc_model.py で作る。捕獲は全 − 弾性 − 核分裂
（UFG の範囲より下では非弾性の反応は捕獲と核分裂だけ）。群ごとの相対差 (xlc / OpenMC − 1) を、
UFG の範囲（WIMS 69 群の 15〜27 群）と、その下の 4〜0.1 eV（28〜56 群）について表示する。
反応率の差は OpenMC の燃料の中性子束で ΔR = Σ_g φ_g N (σ_xlc − σ_OpenMC) と足す。k = νF / A（F は核分裂率、
A はセル全体の吸収率）なので、中性子束を固定した 1 次の近似で Δk/k ≈ ΔνF / νF − (ΔC + ΔF) / A に換算する
（ν はライブラリの無限希釈の群の値）。全 − 弾性 − 核分裂が捕獲に等しい 15 群（9.1 keV）より下だけを足す。核種は OpenMC の名前で、省略すると Pu240, Pu239, Pu241, Pu242, U238。
"""

import sys
from pathlib import Path

import msgpack

UFG = range(15, 28)
BELOW = range(28, 57)
THERMAL = range(57, 70)


def read(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes())


def reactions(rates: dict) -> dict[str, list[float]]:
    capture = [
        t - e - f for t, e, f in zip(rates["total"], rates["elastic"], rates["fission"])
    ]
    return {"capture": capture, "fission": rates["fission"]}


def nu(library: dict, nuclide: str) -> list[float]:
    """無限希釈の群ごとの ν（核分裂のない群は 0）。"""
    (data,) = [n for n in library["nuclides"] if n["name"].replace("-", "") == nuclide]
    dilute = data["temperatures"][0]["shielded"][0]
    return [
        p / f if f > 0 else 0.0 for p, f in zip(dilute["nu_fission"], dilute["fission"])
    ]


def main() -> None:
    library = read(sys.argv[1])
    bounds = library["group_bounds"]
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

    print("\nrate differences (sum over groups, OpenMC flux) and first-order k effect")
    print("nuclide  groups  capture (%)  fission (%)  dk/k (pcm)")
    total = {}
    density = reference["densities"]
    production = reference["nu-fission"]
    for n in nuclides:
        a, b = reactions(ours[n]), reactions(reference["nuclides"][n])
        v = nu(library, n)
        for label, groups in [
            ("15-27", UFG),
            ("28-56", BELOW),
            ("57-69", THERMAL),
            ("15-69", range(UFG.start, THERMAL.stop)),
        ]:
            cells, dk = [], 0.0
            for score in ["capture", "fission"]:
                ref = sum(phi[g - 1] * density[n] * b[score][g - 1] for g in groups)
                diff = [
                    phi[g - 1] * density[n] * (a[score][g - 1] - b[score][g - 1])
                    for g in groups
                ]
                cells.append(
                    f"{sum(diff) / ref * 100:+11.2f}" if ref > 0 else f"{'-':>11}"
                )
                dk -= sum(diff) / absorption
                if score == "fission":
                    dk += sum(v[g - 1] * d for g, d in zip(groups, diff)) / production
            print(f"{n:7s}  {label:6s}  {cells[0]}  {cells[1]}  {dk * 1e5:+10.1f}")
            total[label] = total.get(label, 0.0) + dk
    for label, dk in total.items():
        print(f"{'sum':7s}  {label:6s}  {'':11s}  {'':11s}  {dk * 1e5:+10.1f}")


if __name__ == "__main__":
    main()
