"""VERA Problem 1B〜1D, 2B〜2D を、xlc と同じ多群定数で OpenMC の多群モードで解く（MOC の輸送の検証）。

使い方: python3 openmc_mg.py <多群定数 (msgpack)> <問題> <作業ディレクトリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

多群定数は crates/lattice/examples/mg_materials.rs で、同じ問題の格子の入力から作る（無限希釈、散乱は P0〜P_L の
ルジャンドルのモーメント）。形状は openmc_model.py と同じ。2B〜2D は lattice と同じ形式のピンごとの核分裂率を
pins.msgpack に書く。

OpenMC は散乱角の分布の負の部分を 0 に切るので、xlc と同じ散乱核を解くには、ライブラリを
positive_library.py で負にならない形にしてから多群定数を作る。散乱角は表に直さず、ルジャンドル展開から
直接サンプリングする（表の点の数による近似を入れない）。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc
from openmc_model import geometry_of, pin_mesh, write_pins

# 物質の名前（xlc の格子の入力と同じ）
FUEL, CLAD, WATER = "fuel", "zircaloy", "water"


def library(path: Path, work: Path) -> Path:
    """多群定数を OpenMC の MGXS ライブラリ（HDF5）にして、そのパスを返す。"""
    data = msgpack.unpackb(path.read_bytes())
    if data.get("format") != "xlc-mg-materials":
        raise SystemExit(f"{path}: not an xlc-mg-materials file")
    # OpenMC の境界は昇順。群の番号は xlc と同じく高いエネルギーから
    structure = openmc.mgxs.EnergyGroups(sorted(data["group_bounds"]))
    mgxs = openmc.MGXSLibrary(structure)
    for name, xs in data["materials"].items():
        total = np.array(xs["total"])
        # 散乱で生まれる中性子の行列（(n,2n) などの多重度を含む）[l][入射][出射]
        production = np.array(xs["scatter"])
        # OpenMC は散乱の反応と多重度を分けて持つ。生まれる中性子の和が全断面積を超える群（(n,2n) の効く
        # 高速群）では、散乱を全断面積（吸収 0）とし、残りを行ごとに一様な多重度にする。衝突あたりの
        # 出射の重みの期待値は νΣs(g→g')/Σt のままなので、分け方は輸送の結果に効かない
        made = production[0].sum(axis=1)
        multiplicity = np.where(made > total, made / total, 1.0)
        scatter = production / multiplicity[np.newaxis, :, np.newaxis]
        x = openmc.XSdata(name, structure)
        x.order = production.shape[0] - 1
        x.set_total(total)
        x.set_absorption(total - scatter[0].sum(axis=1))
        # OpenMC の並びは [入射][出射][l]
        x.set_scatter_matrix(np.moveaxis(scatter, 0, -1))
        x.set_multiplicity_matrix(
            np.repeat(multiplicity[:, np.newaxis], len(total), axis=1)
        )
        nu_fission = np.array(xs["nu_fission"])
        if nu_fission.any():
            x.set_fission(np.array(xs["fission"]))
            x.set_nu_fission(nu_fission)
        mgxs.add_xsdata(x)
    file = work / "mgxs.h5"
    mgxs.export_to_hdf5(str(file))
    return file


def main() -> None:
    constants = Path(sys.argv[1]).resolve()
    problem = sys.argv[2]
    work = Path(sys.argv[3]).resolve()
    particles = int(sys.argv[4]) if len(sys.argv) > 4 else 100_000
    batches = int(sys.argv[5]) if len(sys.argv) > 5 else 300
    inactive = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    work.mkdir(parents=True, exist_ok=True)

    def macroscopic(name: str) -> openmc.Material:
        m = openmc.Material(name=name)
        m.set_density("macro", 1.0)
        m.add_macroscopic(name)
        return m

    fuel, clad, water = (macroscopic(n) for n in (FUEL, CLAD, WATER))
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(library(constants, work))
    geometry, half = geometry_of(problem, fuel, clad, water)

    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.tabular_legendre = {"enable": False}
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
        constraints={"fissionable": True},
    )
    tallies = openmc.Tallies()
    if problem.startswith("2"):
        tally = openmc.Tally(name="fission")
        tally.filters = [openmc.MeshFilter(pin_mesh())]
        tally.scores = ["fission"]
        tallies.append(tally)

    model = openmc.Model(geometry, materials, settings, tallies)
    statepoint = model.run(cwd=work, output=False)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        if problem.startswith("2"):
            write_pins(sp, work / "pins.msgpack")
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
