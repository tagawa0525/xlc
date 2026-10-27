"""#0133 の試行: xlc の格子の入力（lattice.toml）の燃料の組成で 1C のピンセルを OpenMC で解き、物質（燃料、被覆管、減速材）ごとの
群ごとの中性子束と、群の間の散乱の移行の行列（nu-scatter、出来事ごとの推定）を、`swap_transfer` の例が読む形（balance.msgpack の
`k` と `materials` の `flux`、`transfer`）で書く。散乱の移行の差の k への寄与を、物質と入射の群の範囲ごとに測るため。同じファイルに、
燃料の群ごとの中性子束と核種ごとの反応率（`(n,gamma)`、`fission`、`nu-fission` だけ）を `swap_absorption` の例が読むキー（`flux`、
`bounds`、`densities`、`rates` の `mean`）で書く（吸収の置き換えも同じ実行の値で測れるように）。標準偏差のキーは `std_dev` で、
openmc_pin_rates.py の `std` と違い、物質ごとの吸収も無いので、openmc_pin_rates.py の出力を読むほかのツールには渡せない。

使い方（references の devShell で）: python3 openmc_pin_transfer.py <openmc-1C-rings1.toml> <lattice.toml> <xlc のライブラリ>
    <作業ディレクトリ> <粒子数> <世代数> <非アクティブ> <スレッド数> <シード>（ここまですべて要る）
    [<細かい入射の下端 eV> <上端 eV> <1 桁あたりの区間の数>]

細かい入射の範囲（3 つとも書くか、どれも書かない。下端と上端は群の境界）を書くと、その範囲の群をそれぞれレサジーの等間隔の
区間（群の幅の桁数 × 1 桁あたりの数を丸めた数、1 以上）に分け、物質ごとの区間の中性子束と、燃料と減速材の区間から群への
移行（nu-scatter）を `fine` に書く（#0152。群の行を、群の中の入射のエネルギーごとの行と、その重みの中性子束に分けて、
群の中のスペクトルの違いが行をどれだけ変えるかを測るため）。区間の行の、区間の中性子束での和が群の行になる。

模型は `openmc_depletion.build`（燃料 1 リング）で、燃料の組成だけを lattice.toml の `[material.fuel]` の `nuclides` に置き換える。
群構造は xlc のライブラリのもので、群 1 の上端は 10 MeV より上に広げる（xlc は 10 MeV より上に生まれた中性子を群 1 に入れる。
openmc_pin_rates.py の ABOVE と同じ）。値は世代の生成で規格化せず、tally の平均のまま書く（`swap_transfer` は移行 / 中性子束の
比だけを使う）。
"""

import itertools
import math
import os
import sys
from pathlib import Path

import msgpack
import openmc
import tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
sys.path.insert(0, str(ROOT / "examples" / "vera"))
import openmc_depletion
from compare_depletion import openmc_name
from openmc_common import group_filter

ABOVE = 1.0e9
SCORES = ["(n,gamma)", "fission", "nu-fission"]
# OpenMC の物質の名前から、xlc の格子の入力の物質の名前へ（間隙は真空で物質を持たない）
NAMES = {"fuel-0": "fuel", "zircaloy-4": "zircaloy", "moderator": "water"}
# 細かい入射の移行を数える物質（xlc の名前）
FINE_MATERIALS = ("fuel", "water")


def fine_edges(bounds, low, high, per_decade):
    """群の境界 `bounds`（昇順）のうち [low, high] の群を、レサジーの等間隔の区間に分けた境界（昇順）。"""
    inside = [b for b in bounds if low * (1 - 1e-9) <= b <= high * (1 + 1e-9)]
    if len(inside) < 2 or not (
        math.isclose(inside[0], low, rel_tol=1e-9)
        and math.isclose(inside[-1], high, rel_tol=1e-9)
    ):
        raise SystemExit(f"{low}-{high} eV must start and end at group bounds")
    edges = [inside[0]]
    for a, b in itertools.pairwise(inside):
        count = max(1, round(math.log10(b / a) * per_decade))
        edges += [a * (b / a) ** ((j + 1) / count) for j in range(count - 1)] + [b]
    return edges


def main():
    if len(sys.argv) not in (10, 13):
        raise SystemExit(__doc__)
    source, lattice, library, work = sys.argv[1:5]
    particles, batches, inactive, threads, seed = map(int, sys.argv[5:10])
    fine_range = (
        (float(sys.argv[10]), float(sys.argv[11]), int(sys.argv[12]))
        if len(sys.argv) == 13
        else None
    )
    run = openmc_depletion.read(Path(source).resolve())
    if run["rings"] != 1:
        raise SystemExit(f"the model must have one fuel ring, got {run['rings']}")
    run.update(particles=particles, batches=batches, inactive=inactive, seed=seed)
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    composition = tomllib.loads(Path(lattice).read_text())["material"]["fuel"][
        "nuclides"
    ]
    model, fuels = openmc_depletion.build(run)
    (fuel,) = fuels
    for name in list(fuel.get_nuclides()):
        fuel.remove_nuclide(name)
    for name, n in composition.items():
        fuel.add_nuclide(openmc_name(name), n)
    fuel.set_density("sum")
    order = list(model.materials)
    unknown = [m.name for m in order if m.name not in NAMES]
    if unknown:
        raise SystemExit(f"materials without an xlc name: {unknown}")
    bounds = list(group_filter(Path(library).resolve()).values)
    groups = openmc.EnergyFilter(bounds[:-1] + [ABOVE])
    out_groups = openmc.EnergyoutFilter(bounds[:-1] + [ABOVE])
    materials = openmc.MaterialFilter(order)
    flux = openmc.Tally(name="flux")
    flux.filters = [materials, groups]
    flux.scores = ["flux"]
    # 出来事ごとの推定が要る score は別の tally にする（同じ tally だと全体が出来事ごとの推定になる）
    transfer = openmc.Tally(name="transfer")
    transfer.filters = [materials, groups, out_groups]
    transfer.scores = ["nu-scatter"]
    nuclides = [openmc_name(name) for name in composition]
    rates = openmc.Tally(name="rates")
    rates.filters = [openmc.MaterialFilter([fuel]), groups]
    rates.nuclides = nuclides
    rates.scores = SCORES
    tallies = [flux, transfer, rates]
    if fine_range is not None:
        fine_bounds = fine_edges(bounds, *fine_range)
        fine_groups = openmc.EnergyFilter(fine_bounds)
        fine_flux = openmc.Tally(name="fine-flux")
        fine_flux.filters = [materials, fine_groups]
        fine_flux.scores = ["flux"]
        fine_transfer = openmc.Tally(name="fine-transfer")
        scattering = [m for m in order if NAMES[m.name] in FINE_MATERIALS]
        fine_transfer.filters = [
            openmc.MaterialFilter(scattering),
            fine_groups,
            out_groups,
        ]
        fine_transfer.scores = ["nu-scatter"]
        tallies += [fine_flux, fine_transfer]
    model.tallies = openmc.Tallies(tallies)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    openmc.config["cross_sections"] = run["cross-sections"]
    sp_path = model.run(cwd=work, output=False)
    n = len(bounds) - 1
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        phi = sp.get_tally(name="flux")
        moved = sp.get_tally(name="transfer")
        out = {
            "format": "xlc-balance",
            "k": k.nominal_value,
            "k_std": k.std_dev,
            "bounds": sorted(bounds, reverse=True),
            "materials": {},
            "densities": dict(composition),
            "rates": {},
        }
        t = sp.get_tally(name="rates")
        for name in composition:
            out["rates"][name] = {
                score: {
                    value: t.get_values(
                        scores=[score], nuclides=[openmc_name(name)], value=value
                    )
                    .ravel()[::-1]
                    .tolist()
                    for value in ("mean", "std_dev")
                }
                for score in SCORES
            }
        for i, m in enumerate(order):
            values = {
                value: (
                    phi.get_values(value=value).reshape(len(order), n)[i][::-1],
                    moved.get_values(value=value).reshape(len(order), n, n)[i][
                        ::-1, ::-1
                    ],
                )
                for value in ("mean", "std_dev")
            }
            # 群は高いエネルギーから。行列は [入射][出射]
            out["materials"][NAMES[m.name]] = {
                "flux": values["mean"][0].tolist(),
                "flux_std": values["std_dev"][0].tolist(),
                "transfer": values["mean"][1].tolist(),
                "transfer_std": values["std_dev"][1].tolist(),
            }
        if fine_range is not None:
            nf = len(fine_bounds) - 1
            f = sp.get_tally(name="fine-flux")
            ft = sp.get_tally(name="fine-transfer")
            # 区間も群と同じく高いエネルギーから。行列は [区間][出射の群]
            out["fine"] = {"bounds": sorted(fine_bounds, reverse=True), "materials": {}}
            for i, m in enumerate(order):
                out["fine"]["materials"][NAMES[m.name]] = {
                    value: f.get_values(value=value)
                    .reshape(len(order), nf)[i][::-1]
                    .tolist()
                    for value in ("mean", "std_dev")
                }
            for i, m in enumerate(scattering):
                entry = out["fine"]["materials"][NAMES[m.name]]
                for value in ("mean", "std_dev"):
                    entry[f"transfer_{value}"] = (
                        ft.get_values(value=value)
                        .reshape(len(scattering), nf, n)[i][::-1, ::-1]
                        .tolist()
                    )
    out["flux"] = out["materials"]["fuel"]["flux"]
    (work / "transfer.msgpack").write_bytes(msgpack.packb(out))
    print(f"k {k}")


main()
