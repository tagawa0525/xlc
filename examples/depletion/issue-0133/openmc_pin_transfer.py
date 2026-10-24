"""#0133 の試行: xlc の格子の入力（lattice.toml）の燃料の組成で 1C のピンセルを OpenMC で解き、物質（燃料、被覆管、減速材）ごとの
群ごとの中性子束と、群の間の散乱の移行の行列（nu-scatter、出来事ごとの推定）を、`swap_transfer` の例が読む形（balance.msgpack の
`k` と `materials` の `flux`、`transfer`）で書く。散乱の移行の差の k への寄与を、物質と入射の群の範囲ごとに測るため。同じファイルに、
燃料の群ごとの中性子束と核種ごとの反応率（`(n,gamma)`、`fission`、`nu-fission`）を `swap_absorption` の例が読む形（openmc_pin_rates.py
の `flux`、`bounds`、`densities`、`rates`）で書く（吸収の置き換えも同じ実行の値で測れるように）。

使い方（references の devShell で）: python3 openmc_pin_transfer.py <openmc-1C-rings1.toml> <lattice.toml> <xlc のライブラリ>
    <作業ディレクトリ> <粒子数> <世代数> <非アクティブ> <スレッド数> <シード>（すべて要る）

模型は `openmc_depletion.build`（燃料 1 リング）で、燃料の組成だけを lattice.toml の `[material.fuel]` の `nuclides` に置き換える。
群構造は xlc のライブラリのもので、群 1 の上端は 10 MeV より上に広げる（xlc は 10 MeV より上に生まれた中性子を群 1 に入れる。
openmc_pin_rates.py の ABOVE と同じ）。値は世代の生成で規格化せず、tally の平均のまま書く（`swap_transfer` は移行 / 中性子束の
比だけを使う）。
"""

import os
import re
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
from openmc_common import group_filter

ABOVE = 1.0e9
SCORES = ["(n,gamma)", "fission", "nu-fission"]
# OpenMC の物質の名前から、xlc の格子の入力の物質の名前へ（間隙は真空で物質を持たない）
NAMES = {"fuel-0": "fuel", "zircaloy-4": "zircaloy", "moderator": "water"}


def openmc_name(name: str) -> str:
    """xlc の核種の名前（U-235、Pm-148m）を OpenMC の名前（U235、Pm148_m1）に"""
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m)?", name).groups()
    return f"{symbol}{mass}" + ("_m1" if isomer else "")


def main():
    if len(sys.argv) != 10:
        raise SystemExit(__doc__)
    source, lattice, library, work = sys.argv[1:5]
    particles, batches, inactive, threads, seed = map(int, sys.argv[5:10])
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
    model.tallies = openmc.Tallies([flux, transfer, rates])
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
    out["flux"] = out["materials"]["fuel"]["flux"]
    (work / "transfer.msgpack").write_bytes(msgpack.packb(out))
    print(f"k {k}")


main()
