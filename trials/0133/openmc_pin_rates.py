"""#0133 の試行: OpenMC の 1C の燃焼の結果のある燃焼点の燃料の組成（輸送する核種だけ）で、1C のピンセルを OpenMC で解き、燃料の
群ごとの中性子束と、核種ごとの群ごとの反応率を書く。xlc 側は同じ組成の格子の入力（lattice.toml、このスクリプトが書く）を
`balance` と `fuel_xs` の例で解く。

使い方（references の devShell で）: python3 openmc_pin_rates.py <openmc-1C-rings1.toml> <燃焼の作業ディレクトリ> <燃焼点の番号>
    <xlc のライブラリ> <xlc の格子の入力のテンプレート> <作業ディレクトリ> <粒子数> <世代数> <非アクティブ> <スレッド数>
"""

import json
import os
import re
import sys
from pathlib import Path

import msgpack
import openmc
import openmc.deplete

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
sys.path.insert(0, str(ROOT / "examples" / "vera"))
import openmc_depletion
from openmc_common import group_filter

SCORES = ["(n,gamma)", "fission", "total", "(n,elastic)", "absorption", "nu-fission"]


def xlc_name(name):
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)(\d+)(_m1)?", name).groups()
    return f"{symbol}-{mass}" + ("m" if isomer else "")


def main():
    source, work_dep, index, library, template, work = sys.argv[1:7]
    particles, batches, inactive, threads = map(int, sys.argv[7:11])
    seed = int(sys.argv[11]) if len(sys.argv) > 11 else 1
    index = int(index)
    run = openmc_depletion.read(Path(source).resolve())
    run.update(particles=particles, batches=batches, inactive=inactive, seed=seed)
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    transportable = json.loads((Path(work_dep) / "transport-seconds.json").read_text())[
        "transportable"
    ]
    results = openmc.deplete.Results(Path(work_dep) / "depletion_results.h5")
    model, fuels = openmc_depletion.build(run)
    fuel = fuels[0]
    mat = str(fuel.id)
    densities = {}
    for name in transportable:
        n = results.get_atoms(mat, name, nuc_units="atom/b-cm")[1][index]
        if n > 0.0:
            densities[name] = float(n)
    for name in list(fuel.get_nuclides()):
        fuel.remove_nuclide(name)
    for name, n in densities.items():
        fuel.add_nuclide(name, n)
    fuel.set_density("sum")
    groups = group_filter(Path(library).resolve())
    mfilter = openmc.MaterialFilter([fuel])
    flux = openmc.Tally(name="flux")
    flux.filters = [mfilter, groups]
    flux.scores = ["flux"]
    rates = openmc.Tally(name="rates")
    rates.filters = [mfilter, groups]
    NUCLIDES = list(densities)
    rates.nuclides = NUCLIDES
    rates.scores = SCORES
    total = openmc.Tally(name="total")
    total.filters = [mfilter]
    total.nuclides = NUCLIDES
    total.scores = ["(n,gamma)", "fission"]
    # 物質ごとの群ごとの吸収と生成（燃料、被覆管、減速材）
    materials = openmc.Tally(name="materials")
    materials.filters = [openmc.MaterialFilter(list(model.materials)), groups]
    materials.scores = ["flux", "absorption", "nu-fission"]
    model.tallies = openmc.Tallies([flux, rates, total, materials])
    os.environ["OMP_NUM_THREADS"] = str(threads)
    openmc.config["cross_sections"] = run["cross-sections"]
    sp_path = model.run(cwd=work, output=False)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        phi = sp.get_tally(name="flux").get_values(scores=["flux"]).ravel()
        phi_std = (
            sp.get_tally(name="flux")
            .get_values(scores=["flux"], value="std_dev")
            .ravel()
        )
        t = sp.get_tally(name="rates")
        out = {
            "k": k.nominal_value,
            "k_std": k.std_dev,
            "flux": phi[::-1].tolist(),
            "flux_std": phi_std[::-1].tolist(),
            "bounds": sorted(groups.values, reverse=True),
            "densities": {xlc_name(n): v for n, v in densities.items()},
            "rates": {},
            "one_group": {},
        }
        tt = sp.get_tally(name="total")
        for name in NUCLIDES:
            out["rates"][xlc_name(name)] = {
                s: {
                    "mean": t.get_values(scores=[s], nuclides=[name])
                    .ravel()[::-1]
                    .tolist(),
                    "std": t.get_values(scores=[s], nuclides=[name], value="std_dev")
                    .ravel()[::-1]
                    .tolist(),
                }
                for s in SCORES
            }
            out["one_group"][xlc_name(name)] = {
                s: [
                    float(tt.get_values(scores=[s], nuclides=[name]).ravel()[0]),
                    float(
                        tt.get_values(
                            scores=[s], nuclides=[name], value="std_dev"
                        ).ravel()[0]
                    ),
                ]
                for s in ["(n,gamma)", "fission"]
            }
        m = sp.get_tally(name="materials")
        out["materials"] = {}
        for material in model.materials:
            out["materials"][material.name] = {
                s: m.get_values(
                    scores=[s],
                    filters=[openmc.MaterialFilter],
                    filter_bins=[(material.id,)],
                )
                .ravel()[::-1]
                .tolist()
                for s in ["flux", "absorption", "nu-fission"]
            }
    (work / "openmc-rates.msgpack").write_bytes(msgpack.packb(out))
    # xlc の格子の入力: テンプレートの燃料の組成を、この組成（xlc の名前）に置き換える
    text = Path(template).read_text()
    composition = ", ".join(f'"{xlc_name(n)}" = {v:.8e}' for n, v in densities.items())
    text = re.sub(
        r"(\[material\.fuel\]\ntemperature = \d+\n)nuclides = \{[^}]*\}",
        lambda m: m.group(1) + f"nuclides = {{ {composition} }}",
        text,
    )
    (work / "lattice.toml").write_text(text)
    print(f"k {k}")


main()
