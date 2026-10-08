"""xlc と OpenMC の燃料の核種ごとの実効捕獲・核分裂断面積を群ごとに比べ、反応率の差を k の差に換算する。

使い方: python3 compare_capture.py <xs.toml> <xlc の fuel-xs.msgpack> <OpenMC の fuel-xs.msgpack> [核種 ...]

xlc 側は crates/lattice/examples/fuel_xs.rs、OpenMC 側は openmc_model.py で作る。`xs.toml` からライブラリ
（群の境界と ν）と核種の ENDF のファイルを読む。核種は OpenMC の名前で、省略すると Pu240, Pu239, Pu241, Pu242, U238。

xlc の出力には捕獲がないので、全 − 弾性 − 核分裂を捕獲として OpenMC の (n,γ) の tally と比べる。両者が等しいのは、
弾性散乱・核分裂・捕獲のほかの反応（非弾性散乱、(n,α) など）がない群だけなので、核種ごとに ENDF の MF3 から
ほかの反応の断面積が 0 でなくなる最も低いエネルギー（しきい値）を求め、上の境界がしきい値以下の群だけを比べる
（Pu-239 は第1準位の 7.9 keV、U-235 は 77 eV の準位があり、UFG の範囲の上の群が外れる）。核分裂はどの群でも比べる。

群ごとの相対差 (xlc / OpenMC − 1) を、UFG の範囲（4 eV〜15.03 keV、WIMS 69 群の 14〜27 群）と、その下の
4〜0.1 eV（28〜56 群）について表示する。括弧の中は OpenMC の実効断面積の相対標準偏差の上限
√(反応率² + 中性子束²)（反応率と中性子束は同じ飛跡から数えるので正に相関し、実際の誤差はこれより小さい）。

反応率の差は OpenMC の燃料の中性子束で ΔR = Σ_g φ_g N (σ_xlc − σ_OpenMC) と足す。k = νF / A（F は核分裂率、
A はセル全体の吸収率）なので、中性子束を固定した 1 次の近似で Δk/k ≈ ΔνF / νF − (ΔC + ΔF) / A に換算する
（ν はライブラリの無限希釈の群の値。吸収は捕獲と核分裂だけで、捕獲は比べた群だけを足す）。14 群（15.03 keV）
より下だけを足す。
"""

import math
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib

# WIMS 69 群（1 始まり）の UFG の範囲（4 eV〜15.03 keV）、4〜0.1 eV、熱群
UFG = range(14, 28)
BELOW = range(28, 57)
THERMAL = range(57, 70)

# しきい値を求めるときに除く MF3 の MT: 和（1 全, 3 非弾性, 4 非弾性散乱, 27 吸収, 101 消滅）、
# 弾性散乱（2）、核分裂（18〜21, 38）、捕獲（102）
EXCLUDED = {1, 2, 3, 4, 18, 19, 20, 21, 27, 38, 101, 102}


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes())


def endf_float(field: str) -> float:
    """ENDF の 11 文字の数（"1.234567+5" のように e を省くことがある）。"""
    text = field.strip()
    if not text:
        return 0.0
    for i in range(len(text) - 1, 0, -1):
        if text[i] in "+-" and text[i - 1] not in "eE":
            return float(f"{text[:i]}e{text[i:]}")
    return float(text)


def fields(line: str) -> list[str]:
    return [line[11 * i : 11 * (i + 1)] for i in range(6)]


Table = tuple[list[float], list[float]]


def competing(path: Path) -> tuple[float, list[Table]]:
    """弾性散乱・核分裂・捕獲のほかの反応（MF3）の、しきい値と、しきい値のない反応の表。

    しきい値は、しきい値のある反応の断面積が 0 でなくなる最も低いエネルギー (eV)（TAB1 で最初に断面積が正になる
    点の一つ前の点。線形に補間すると、そこから上で正になる）。最も低いエネルギーから断面積が正の反応（U-235 と
    U-238 の MT5 と (n,α)、O-16 の (n,p) など、中性子を出さない発熱反応）は、表 (エネルギー, 断面積) を返す。
    """
    lines = path.read_text().splitlines()
    threshold, tables = math.inf, []
    i = 0
    while i < len(lines):
        line = lines[i]
        mf, mt = int(line[70:72].strip() or 0), int(line[72:75].strip() or 0)
        if mf != 3 or mt == 0:
            i += 1
            continue
        # HEAD の次の TAB1 の制御の行: QM, QI, 0, LR, NR, NP
        control = fields(lines[i + 1])
        regions, points = int(control[4]), int(control[5])
        start = i + 2 + math.ceil(2 * regions / 6)
        values = [
            endf_float(f)
            for row in lines[start : start + math.ceil(2 * points / 6)]
            for f in fields(row)
        ][: 2 * points]
        energies, sigmas = values[0::2], values[1::2]
        if (mt < 200 or 600 <= mt < 850) and mt not in EXCLUDED:
            first = next((j for j, s in enumerate(sigmas) if s > 0), None)
            if first == 0:
                tables.append((energies, sigmas))
            elif first is not None:
                threshold = min(threshold, energies[first - 1])
        # 次の MT（SEND の行の後）まで飛ばす
        i = start + math.ceil(2 * points / 6)
    return threshold, tables


def maximum(table: Table, low: float, high: float) -> float:
    """線形に補間した断面積の [low, high] での最大（最大は内側の点か両端にある）。"""
    energies, sigmas = table
    inside = [s for e, s in zip(energies, sigmas) if low <= e <= high]
    ends = [float(np.interp(e, energies, sigmas)) for e in (low, high)]
    return max(inside + ends)


def nu(library: dict, nuclide: str) -> list[float]:
    """無限希釈の群ごとの ν（核分裂のない群は 0）。"""
    (data,) = [n for n in library["nuclides"] if n["name"].replace("-", "") == nuclide]
    dilute = data["temperatures"][0]["shielded"][0]
    return [
        p / f if f > 0 else 0.0 for p, f in zip(dilute["nu_fission"], dilute["fission"])
    ]


def main() -> None:
    toml = Path(sys.argv[1])
    config = tomllib.loads(toml.read_text())
    library = read(toml.parent / config["output"])
    endf = toml.parent / config["endf-directory"]
    bounds = library["group_bounds"]
    ours = read(Path(sys.argv[2]))["nuclides"]
    reference = read(Path(sys.argv[3]))
    nuclides = sys.argv[4:] or ["Pu240", "Pu239", "Pu241", "Pu242", "U238"]
    files = {
        name.replace("-", ""): endf / entry["endf"]
        for name, entry in config["nuclide"].items()
    }
    phi = reference["flux"]
    absorption = reference["absorption"]
    k = reference["keff"]
    print(f"OpenMC k {k[0]:.6f} +- {k[1]:.6f}")

    xlc, openmc = {}, {}
    for n in nuclides:
        rates = ours[n]
        xlc[n] = {
            "capture": [
                t - e - f
                for t, e, f in zip(rates["total"], rates["elastic"], rates["fission"])
            ],
            "fission": rates["fission"],
        }
        openmc[n] = {s: reference["nuclides"][n][s] for s in ["capture", "fission"]}

    # 捕獲を比べる群: 上の境界（bounds[g - 1]、高い順）がしきい値以下の群
    # しきい値のない反応による差 (n,γ) と全 − 弾性 − 核分裂 の違いの上限は、群の中の断面積の最大 / OpenMC の
    # 実効捕獲断面積（比べた群の中の最大を表示する）
    compared = {}
    print(
        "\ncompeting-reaction threshold (eV), first compared capture group, and"
        " bound of the relative effect of reactions without threshold"
    )
    for n in nuclides:
        threshold, tables = competing(files[n])
        compared[n] = {g for g in range(1, len(bounds)) if bounds[g - 1] <= threshold}
        bound = max(
            (
                sum(maximum(t, bounds[g], bounds[g - 1]) for t in tables)
                / openmc[n]["capture"][g - 1]
                for g in compared[n]
                if openmc[n]["capture"][g - 1] > 0
            ),
            default=0.0,
        )
        first = min(compared[n], default="-")
        print(f"{n:7s}  {threshold:12.5g}  {first:>5}  {bound:9.2g}")

    def rsd(n: str, score: str, g: int) -> float:
        rate = reference["nuclides"][n][f"{score}-rsd"][g - 1]
        return math.hypot(rate, reference["flux-rsd"][g - 1])

    def used(n: str, score: str, g: int) -> bool:
        return score == "fission" or g in compared[n]

    for score in ["capture", "fission"]:
        print(
            f"\n{score}: xlc / OpenMC - 1 (%), bound of OpenMC rsd (%) in brackets"
            + (
                " (- where inelastic or other reactions compete)"
                if score == "capture"
                else ""
            )
        )
        print("group  E_low (eV)  " + "  ".join(f"{n:>15}" for n in nuclides))
        for g in [*UFG, *BELOW]:
            cells = []
            for n in nuclides:
                a, b = xlc[n][score][g - 1], openmc[n][score][g - 1]
                cells.append(
                    f"{(a / b - 1) * 100:+7.2f} [{rsd(n, score, g) * 100:4.2f}]"
                    if b > 0 and used(n, score, g)
                    else f"{'-':>15}"
                )
            print(f"{g:5d}  {bounds[g]:10.4g}  " + "  ".join(cells))

    print("\nrate differences (sum over groups, OpenMC flux) and first-order k effect")
    print("nuclide  groups  capture (%)  fission (%)  dk/k (pcm)")
    total = {}
    density = reference["densities"]
    production = reference["nu-fission"]
    for n in nuclides:
        a, b = xlc[n], openmc[n]
        v = nu(library, n)
        for label, groups in [
            ("14-27", UFG),
            ("28-56", BELOW),
            ("57-69", THERMAL),
            ("14-69", range(UFG.start, THERMAL.stop)),
        ]:
            cells, dk = [], 0.0
            for score in ["capture", "fission"]:
                selected = [g for g in groups if used(n, score, g)]
                ref = sum(phi[g - 1] * density[n] * b[score][g - 1] for g in selected)
                diff = [
                    phi[g - 1] * density[n] * (a[score][g - 1] - b[score][g - 1])
                    for g in selected
                ]
                cells.append(
                    f"{sum(diff) / ref * 100:+11.2f}" if ref > 0 else f"{'-':>11}"
                )
                dk -= sum(diff) / absorption
                if score == "fission":
                    dk += sum(v[g - 1] * d for g, d in zip(selected, diff)) / production
            print(f"{n:7s}  {label:6s}  {cells[0]}  {cells[1]}  {dk * 1e5:+10.1f}")
            total[label] = total.get(label, 0.0) + dk
    for label, dk in total.items():
        print(f"{'sum':7s}  {label:6s}  {'':11s}  {'':11s}  {dk * 1e5:+10.1f}")


if __name__ == "__main__":
    main()
