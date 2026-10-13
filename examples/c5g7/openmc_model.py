"""C5G7 の 2D の問題を OpenMC の多群モードで解き、lattice と同じ形式のピンごとの核分裂率を書く。

使い方: python3 openmc_model.py run <入力> <作業ディレクトリ>
        python3 openmc_model.py extract <入力> <作業ディレクトリ> <結果のディレクトリ>

<入力> は openmc-c5g7.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 1 世代の粒子数
（`particles`）、世代数（`batches`）、非アクティブバッチ数（`inactive`）、乱数のシード（`seed`）、スレッド数（`threads`）。
run は OpenMC を解いて作業ディレクトリに statepoint を残し、extract は statepoint から k を k.json に、
ピンごとの核分裂率を pins.msgpack に結果のディレクトリへ書く（ADR-0123）。

断面積は c5g7-xs.toml（lattice と同じファイル）を読んで OpenMC の多群データに変換する。
形状は c5g7.toml と同じく、左上と右下が UO2、右上と左下が MOX の集合体、右と下が水の反射体。
"""

import sys
from pathlib import Path

import numpy as np
import openmc
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import read_input
from openmc_run import RUN_KEYS, fuel_mean, run_or_extract, write_k, write_pin_powers

HERE = Path(__file__).resolve().parent
PITCH = 1.26
RADIUS = 0.54
N = 17


def read_xs(path: Path) -> tuple[int, dict[str, dict[str, np.ndarray]]]:
    """c5g7-xs.toml を物質名 → {total, fission, nu-fission, chi, scatter} の辞書にする。"""
    data = tomllib.loads(path.read_text())
    materials = {
        name: {key: np.array(values) for key, values in entry.items()}
        for name, entry in data["material"].items()
    }
    return data["groups"], materials


# 入力のキー。この問題は断面積が c5g7-xs.toml にあるので、問題と cross_sections.xml は要らない
KEYS = RUN_KEYS - {"problem", "cross-sections"}


def solve(source: Path, work: Path) -> None:
    run = read_input(source, KEYS, set())
    particles, batches, inactive, seed = (
        run[key] for key in ("particles", "batches", "inactive", "seed")
    )
    work.mkdir(parents=True, exist_ok=True)

    groups, data = read_xs(HERE / "c5g7-xs.toml")
    # 群境界は使わないが、OpenMC は昇順の境界を求める（群 1 が最も高いエネルギー）
    structure = openmc.mgxs.EnergyGroups(np.logspace(-5, 7, groups + 1))
    library = openmc.MGXSLibrary(structure)
    materials = {}
    for name, xs in data.items():
        x = openmc.XSdata(name, structure)
        x.order = 0
        x.set_total(xs["total"])
        x.set_absorption(xs["total"] - np.sum(xs["scatter"], axis=1))
        # OpenMC の散乱行列は [入射][出射][ルジャンドル次数]
        x.set_scatter_matrix(np.array(xs["scatter"])[:, :, np.newaxis])
        if "nu-fission" in xs:
            x.set_fission(xs["fission"])
            x.set_nu_fission(xs["nu-fission"])
            x.set_chi(xs["chi"] / xs["chi"].sum())
        library.add_xsdata(x)
        m = openmc.Material(name=name)
        m.set_density("macro", 1.0)
        m.add_macroscopic(name)
        materials[name] = m
    library.export_to_hdf5(str(work / "mgxs.h5"))
    material_list = openmc.Materials(materials.values())
    material_list.cross_sections = str(work / "mgxs.h5")

    circle = openmc.ZCylinder(r=RADIUS)

    def pin(fuel: str) -> openmc.Universe:
        inside = openmc.Cell(fill=materials[fuel], region=-circle)
        outside = openmc.Cell(fill=materials["water"], region=+circle)
        return openmc.Universe(cells=[inside, outside])

    water = openmc.Universe(cells=[openmc.Cell(fill=materials["water"])])
    pins = {name: pin(name) for name in ["UO2", "MOX43", "MOX70", "MOX87", "GT", "FC"]}

    # c5g7.toml と同じ集合体の配置（1: 4.3%, 2: 7.0%, 3: 8.7%, 4: 案内管, 5: 核分裂計数管）
    mox_map = np.array(
        [
            [1] * 17,
            [1] + [2] * 15 + [1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 5, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1],
            [1, 2, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 2, 1],
            [1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1],
            [1, 2, 2, 4, 2, 3, 3, 3, 3, 3, 3, 3, 2, 4, 2, 2, 1],
            [1, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 1],
            [1] + [2] * 15 + [1],
            [1] * 17,
        ]
    )
    mox_names = {1: "MOX43", 2: "MOX70", 3: "MOX87", 4: "GT", 5: "FC"}
    uo2_names = {1: "UO2", 2: "UO2", 3: "UO2", 4: "GT", 5: "FC"}

    def block(names: dict[int, str] | None) -> list[list[openmc.Universe]]:
        if names is None:
            return [[water] * N for _ in range(N)]
        return [[pins[names[v]] for v in row] for row in mox_map]

    uo2, mox, refl = block(uo2_names), block(mox_names), block(None)
    universes = []
    for row_blocks in [(uo2, mox, refl), (mox, uo2, refl), (refl, refl, refl)]:
        for r in range(N):
            universes.append([u for b in row_blocks for u in b[r]])
    width = 3 * N * PITCH
    lattice = openmc.RectLattice()
    lattice.lower_left = (-width / 2, -width / 2)
    lattice.pitch = (PITCH, PITCH)
    # OpenMC の RectLattice の行は上から並べる（c5g7.toml の rows と同じ順）
    lattice.universes = universes

    left = openmc.XPlane(-width / 2, boundary_type="reflective")
    right = openmc.XPlane(width / 2, boundary_type="vacuum")
    bottom = openmc.YPlane(-width / 2, boundary_type="vacuum")
    top = openmc.YPlane(width / 2, boundary_type="reflective")
    root = openmc.Cell(fill=lattice, region=+left & -right & +bottom & -top)
    geometry = openmc.Geometry([root])

    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    settings.source = openmc.IndependentSource(
        # 燃料のある左上の 2 × 2 集合体（x は左から、y は上から 2/3）
        space=openmc.stats.Box((-width / 2, -width / 6, -1), (width / 6, width / 2, 1)),
        constraints={"fissionable": True},
    )

    mesh = openmc.RegularMesh()
    mesh.dimension = (3 * N, 3 * N)
    mesh.lower_left = (-width / 2, -width / 2)
    mesh.upper_right = (width / 2, width / 2)
    tally = openmc.Tally(name="fission")
    tally.filters = [openmc.MeshFilter(mesh)]
    tally.scores = ["fission"]

    model = openmc.Model(geometry, material_list, settings, openmc.Tallies([tally]))
    model.run(cwd=work, threads=run["threads"], output=False)


def extract(source: Path, work: Path, results: Path) -> None:
    run = read_input(source, KEYS, set())
    results.mkdir(parents=True, exist_ok=True)
    with openmc.StatePoint(work / f"statepoint.{run['batches']}.h5") as sp:
        k = sp.keff
        values = sp.get_tally(name="fission").get_values(scores=["fission"]).ravel()
    # メッシュは x が速く変わり、y は下から上。lattice と同じく上の行から並べる
    rates = values.reshape(3 * N, 3 * N)[::-1, :]
    write_pin_powers(results / "pins.msgpack", rates / fuel_mean(rates))
    write_k(results / "k.json", "C5G7", k)


if __name__ == "__main__":
    run_or_extract(solve, extract, __doc__)
