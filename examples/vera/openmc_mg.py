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
from openmc_model import (
    COLD,
    FUEL_TEMPERATURE,
    MODERATOR_OF,
    MODERATOR_TEMPERATURE,
    geometry_of,
    pin_mesh,
    write_pins,
)
from positive_library import minimum

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import eigenvalue_settings

# 物質の名前（xlc の格子の入力と同じ）
FUEL, CLAD, WATER = "fuel", "zircaloy", "water"
# 断面積が全て 0 の物質の名前
VOID = "void"


class Run:
    """コマンドラインの 2〜7 番目の引数: <問題> <作業ディレクトリ> [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]。"""

    def __init__(self, argv: list[str]) -> None:
        self.problem = argv[2]
        if any(
            self.problem not in table
            for table in (FUEL_TEMPERATURE, MODERATOR_TEMPERATURE, MODERATOR_OF)
        ):
            raise SystemExit(
                f"unknown problem {self.problem}; expected one of {list(FUEL_TEMPERATURE)}"
            )
        self.work = Path(argv[3]).resolve()
        self.particles = int(argv[4]) if len(argv) > 4 else 100_000
        self.batches = int(argv[5]) if len(argv) > 5 else 300
        self.inactive = int(argv[6]) if len(argv) > 6 else 50
        # 統計誤差の見積もりを確かめるには、種を変えて解き直す
        self.seed = int(argv[7]) if len(argv) > 7 else 1

    def settings(self, half: float) -> openmc.Settings:
        """世代の数と乱数の種、燃料に一様な初期の中性子源（`half` は外周の箱の半分の幅）。"""
        return eigenvalue_settings(
            self.particles, self.batches, self.inactive, self.seed, half
        )


def check_positive(name: str, scatter: np.ndarray) -> None:
    """散乱角の分布が負になる移行があれば止める（OpenMC は負の部分を 0 に切るので、xlc と違う散乱核になる）。"""
    weights = (2 * np.arange(scatter.shape[0]) + 1) / 2
    if (scatter[0] < 0).any() or scatter[1:, scatter[0] == 0].any():
        raise SystemExit(f"{name}: negative or moment-only scattering transfers")
    for g, h in zip(*np.nonzero(scatter[0] > 0)):
        if minimum(weights * scatter[:, g, h] / scatter[0, g, h]) < -1e-12:
            raise SystemExit(
                f"{name}: the scattering kernel from group {g} to {h} is negative "
                "(make the library with positive_library.py at this order)"
            )


def library(path: Path, problem: str, work: Path) -> Path:
    """多群定数を OpenMC の MGXS ライブラリ（HDF5）にして、そのパスを返す。物質の温度が問題と合うこと。"""
    data = msgpack.unpackb(path.read_bytes())
    if data.get("format") != "xlc-mg-materials":
        raise SystemExit(f"{path}: not an xlc-mg-materials file")
    expected = {FUEL: FUEL_TEMPERATURE[problem], CLAD: COLD, WATER: COLD}
    if data["temperatures"] != expected:
        raise SystemExit(
            f"{path}: temperatures {data['temperatures']} are not those of {problem} {expected}"
        )
    # OpenMC の境界は昇順。群の番号は xlc と同じく高いエネルギーから
    structure = openmc.mgxs.EnergyGroups(sorted(data["group_bounds"]))
    ng = structure.num_groups
    for name, xs in data["materials"].items():
        shapes = [
            np.shape(xs[k]) for k in ("total", "fission", "nu_fission", "scatter")
        ]
        if shapes[:3] != [(ng,), (ng,), (ng, ng)] or shapes[3][1:] != (ng, ng):
            raise SystemExit(f"{name}: arrays {shapes} do not match {ng} groups")
    mgxs = openmc.MGXSLibrary(structure)
    for name, xs in data["materials"].items():
        total = np.array(xs["total"])
        # 散乱で生まれる中性子の行列（(n,2n) などの多重度を含む）[l][入射][出射]
        production = np.array(xs["scatter"])
        check_positive(name, production)
        # OpenMC は散乱の反応と多重度を分けて持つ。生まれる中性子の和が全断面積を超える群（(n,2n) の効く
        # 高速群）では、散乱を全断面積（吸収 0）とし、残りを行ごとに一様な多重度にする。衝突あたりの
        # 出射の重みの期待値は νΣs(g→g')/Σt のままなので、分け方は輸送の結果に効かない
        made = production[0].sum(axis=1)
        multiplicity = np.where(made > total, made / total, 1.0)
        scatter = production / multiplicity[np.newaxis, :, np.newaxis]
        x = openmc.XSdata(name, structure)
        x.order = production.shape[0] - 1
        x.set_total(total)
        # 多重度で割った行は丸めで吸収がわずかに負になりうるので、ちょうど 0 にする
        x.set_absorption(np.where(made > total, 0.0, total - scatter[0].sum(axis=1)))
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
    # 真空のギャップの代わりの、断面積が全て 0 の物質（`main` を参照）
    void = openmc.XSdata(VOID, structure)
    void.order = 0
    void.set_total(np.zeros(ng))
    void.set_absorption(np.zeros(ng))
    void.set_scatter_matrix(np.zeros((ng, ng, 1)))
    mgxs.add_xsdata(void)
    file = work / "mgxs.h5"
    mgxs.export_to_hdf5(str(file))
    return file


def main() -> None:
    constants = Path(sys.argv[1]).resolve()
    run = Run(sys.argv)
    problem, work = run.problem, run.work
    work.mkdir(parents=True, exist_ok=True)

    def macroscopic(name: str) -> openmc.Material:
        m = openmc.Material(name=name)
        m.set_density("macro", 1.0)
        m.add_macroscopic(name)
        return m

    fuel, clad, water, void = (macroscopic(n) for n in (FUEL, CLAD, WATER, VOID))
    materials = openmc.Materials([fuel, clad, water, void])
    materials.cross_sections = str(library(constants, problem, work))
    # ギャップを物質のない真空にすると、OpenMC 0.15.2 の多群モードは、核分裂の tally を数えるときに
    # 物質の番号 -1 で断面積の配列を引いて落ちる（src/tallies/tally_scoring.cpp の score_general_mg）。
    # 断面積が全て 0 の物質で埋める（衝突しないので真空と同じ）
    geometry, half = geometry_of(problem, fuel, clad, water, void)

    settings = run.settings(half)
    settings.energy_mode = "multi-group"
    settings.tabular_legendre = {"enable": False}
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
