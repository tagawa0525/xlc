"""#0133 の試行（issue-0133/openmc_pin_rates.py の組成の出どころだけを変えた版）: OpenMC の 2A の燃焼の結果の燃焼点の燃料の組成
（39 の組の体積の重みの平均、1C の燃焼で輸送する核種だけ）で、1C のピンセルを OpenMC で解き、燃料の
群ごとの中性子束と、核種ごとの群ごとの反応率を書く。群の中の重みを見るため、燃料の細かいエネルギーの中性子束と
`FINE_NUCLIDES` の捕獲も書く（`fine`）。xlc 側は同じ組成の格子の入力（lattice.toml、このスクリプトが書く）を
`balance` と `fuel_xs` の例で解く。

使い方（references の devShell で）: python3 openmc_pin_rates_2a.py <openmc-1C-rings1.toml> <燃焼の作業ディレクトリ> <燃焼点の番号>
    <xlc のライブラリ> <xlc の格子の入力のテンプレート> <作業ディレクトリ> <粒子数> <世代数> <非アクティブ> <スレッド数> <シード> <openmc-2A.msgpack>（すべて要る）
"""

import json
import math
import os
import re
import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc
import openmc.deplete

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
sys.path.insert(0, str(ROOT / "examples" / "vera"))
import openmc_depletion
from openmc_common import group_filter

SCORES = ["(n,gamma)", "fission", "total", "(n,elastic)", "absorption", "nu-fission"]
# 群構造の上端（10 MeV）より上の中性子束と反応は、群 1 に数える（openmc_pin_rates.py の ABOVE と同じ）
ABOVE = 1.0e9
# 群の中の重みを見る細かいエネルギーの tally（1e-3 eV〜20 keV をレサジーで等分、1 桁に 500 区間）の、捕獲を数える核種
FINE_NUCLIDES = [
    "U236",
    "U238",
    "Pu239",
    "Pu240",
    "Pu241",
    "Pu242",
    "Am241",
    "Am243",
    "Tc99",
    "Ag109",
    "Nd145",
    "Sm149",
    "Sm152",
    "Xe135",
]
FINE_BOUNDS = (1e-3, 2e4, 500)


def xlc_name(name):
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)(\d+)(_m1)?", name).groups()
    return f"{symbol}-{mass}" + ("m" if isomer else "")


def main():
    source, work_dep, index, library, template, work = sys.argv[1:7]
    particles, batches, inactive, threads = map(int, sys.argv[7:11])
    if len(sys.argv) != 13:
        raise SystemExit(__doc__)
    seed = int(sys.argv[11])
    index = int(index)
    run = openmc_depletion.read(Path(source).resolve())
    run.update(particles=particles, batches=batches, inactive=inactive, seed=seed)
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    transportable = json.loads((Path(work_dep) / "transport-seconds.json").read_text())[
        "transportable"
    ]
    # 組成は 2A の結果（引数 12）から取る。work_dep は輸送する核種の一覧にだけ使う
    assembly = msgpack.unpackb(Path(sys.argv[12]).read_bytes())
    if assembly["format"] != "openmc-depletion-assembly":
        raise SystemExit(f"{sys.argv[12]}: not an OpenMC assembly depletion result")
    if not 0 <= index < len(assembly["burnup"]):
        raise SystemExit(f"point {index} is not in {len(assembly['burnup'])} points")
    print(f"burnup {assembly['burnup'][index]} MWd/kg")
    weights = [r["volume"] for r in assembly["regions"]]
    total_volume = sum(weights)
    model, fuels = openmc_depletion.build(run)
    fuel = fuels[0]
    densities = {}
    # 2A の結果に無い核種（1C だけで輸送する核種）は、両方の組成から外す。外した核種を出す
    missing = [n for n in transportable if n not in assembly["regions"][0]["densities"]]
    print(f"not in the 2A result: {missing}")
    for name in transportable:
        if name in missing:
            continue
        n = (
            sum(
                w * r["densities"][name][index]
                for w, r in zip(weights, assembly["regions"])
            )
            / total_volume
        )
        if n > 0.0:
            densities[name] = float(n)
    for name in list(fuel.get_nuclides()):
        fuel.remove_nuclide(name)
    for name, n in densities.items():
        fuel.add_nuclide(name, n)
    fuel.set_density("sum")
    # 群 1 の上端を広げ、10 MeV より上の中性子束と反応を群 1 に数える（openmc_pin_rates.py と同じ）
    bounds = list(group_filter(Path(library).resolve()).values)
    groups = openmc.EnergyFilter(bounds[:-1] + [ABOVE])
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
    low, high, per_decade = FINE_BOUNDS
    decades = math.log10(high / low)
    fine_filter = openmc.EnergyFilter(
        low * 10 ** (np.arange(round(decades * per_decade) + 1) / per_decade)
    )
    fine = openmc.Tally(name="fine")
    fine.filters = [mfilter, fine_filter]
    fine.scores = ["flux"]
    fine_capture = openmc.Tally(name="fine-capture")
    fine_capture.filters = [mfilter, fine_filter]
    fine_capture.nuclides = [n for n in FINE_NUCLIDES if n in densities]
    fine_capture.scores = ["(n,gamma)"]
    model.tallies = openmc.Tallies([flux, rates, total, materials, fine, fine_capture])
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
            "bounds": sorted(bounds, reverse=True),
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
        f = sp.get_tally(name="fine")
        fc = sp.get_tally(name="fine-capture")
        # 細かい tally はエネルギーの昇順のまま書く
        out["fine"] = {
            "bounds": fine_filter.values.tolist(),
            "flux": f.get_values(scores=["flux"]).ravel().tolist(),
            "capture": {
                xlc_name(n): fc.get_values(scores=["(n,gamma)"], nuclides=[n])
                .ravel()
                .tolist()
                for n in fc.nuclides
            },
        }
    (work / "openmc-rates.msgpack").write_bytes(msgpack.packb(out))
    # xlc の格子の入力: テンプレートの燃料の組成を、この組成（xlc の名前）に置き換える
    text = Path(template).read_text()
    composition = ", ".join(f'"{xlc_name(n)}" = {v:.8e}' for n, v in densities.items())
    text, count = re.subn(
        r"(\[material\.fuel\]\ntemperature = \d+\n)nuclides = \{[^}]*\}",
        lambda m: m.group(1) + f"nuclides = {{ {composition} }}",
        text,
    )
    if count != 1:
        raise SystemExit(f"{template}: no fuel composition to replace")
    (work / "lattice.toml").write_text(text)
    print(f"k {k}")


main()
