"""VERA Problem 3A と 3B（3D の 17×17 集合体、零出力）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model_3d.py run [--restart] [--threads N] <入力> <作業ディレクトリ>
        python3 openmc_model_3d.py extract <入力> <作業ディレクトリ> <結果のディレクトリ>

<入力> は p3/openmc-3A.toml の形の TOML で、openmc_model.py の入力と同じキー（`problem`、`cross-sections`、`particles`、
`batches`、`inactive`、`seed`、`threads`）をすべて書く（既定値はない。ADR-0130）。途中の statepoint を書く間隔は、任意のキー `statepoint-interval`（世代数）。
`--restart` と `--threads` は openmc_resume.py。run と extract の分け方も openmc_model.py
と同じ（ADR-0123）。

径方向の形（ピン、案内管、計装管、Pyrex の棒の寸法と配置、支持格子の枠の置き方）と、問題 2 と同じ組成は
openmc_model.py から import する。軸方向の層の並びと、層ごとにピンの中身を選ぶ部分はここに書く（openmc_model.py の
形状は z に一様な 2D なので）。

軸方向（z は下部炉心板の上面から、cm。仕様書の 1.1〜1.7 節、表 1、3〜6、9、図 P3-1）:
- 下から水 50 cm、下部炉心板（−5〜0）、下部ノズル（0〜6.053）、上部ノズル（397.51〜406.337）、上部炉心板（〜413.937）、
  水 50 cm。炉心板とノズルは表 P3-3 の均質な組成で、集合体の外周の水のギャップを含む箱全体を占める。上下の端は真空
- 燃料棒は 10.281〜395.381: 下の端栓（Zircaloy-4 の中実の棒、1.67 cm）、燃料（11.951〜377.711）、プレナム（16.0 cm、
  被覆管の中は真空。ばねは仕様書のとおり除く）、上の端栓（1.67 cm）
- 案内管と計装管は下部ノズルの上から上部ノズルの下まで（6.053〜397.51）
- 支持格子は、表 3 の質量を 289 のセルに等しく分け、セルの外周の正方形の帯の枠にする（2Q と同じ。端の 2 つは
  Inconel-718、中間の 6 つは Zircaloy-4）
- 3B の Pyrex の棒（16 本、図 5）は、下の端栓（SS304 の中実の棒、半径は被覆管の外側、2.54 cm）、毒物（15.761〜376.441 cm、
  高さ 360.68 cm）、プレナム（22.2 cm、被覆管の中は真空）で、上部ノズルの下端で切る（ノズルの中は均質な組成）。
  端栓の高さ 2.54 cm は下の端栓とする（毒物の上のプレナムの上端 398.641 cm が上部ノズルの中にあるので）。
  プレナムの内管は寸法がないので置かない
- 3B のシンブルプラグ（Pyrex の棒のない 8 本の案内管）は SS304 の中実の棒で、表 9 の位置 383.31 cm を下端とする
  （383.31〜394.31。表 9 の位置を下端とするのは、燃料と Pyrex の毒物の位置がどれも下端なので）

径方向の外周は鏡面反射（集合体 1 体、外周の水のギャップ 0.04 cm を含む）。Hf とギャップの He は除いて真空とし
（ADR-0054）、Pyrex の棒の中の He も同じとする（ADR-0120）。3A は 600 K、3B は 565 K の等温（表 P3-1）。

extract は結果のディレクトリに次を書く:
- k.json: k と統計誤差
- pins.msgpack: 燃料の全長で積分したピンの核分裂率（`xlc-pin-powers`。行は上から、燃料のピンの平均で 1）
- pins-std.msgpack: その統計誤差（同じ形式と規格化。規格化の定数の誤差は含めない）
- layers.msgpack: 仕様書の表 P3-5 の 49 層ごとの出力（`format` が "xlc-axial-pin-powers"）。`bounds` は層の境界の z
  （50 個、下から）、`axial` は層ごとの径方向に積分した出力の線密度（燃料の全長での平均で 1）、`layers` は層ごとの
  ピンの出力の線密度（[層（下から）][行（上から）][列]。燃料のピンと全長での平均で 1）、`axial-std` と `layers-std` は
  その統計誤差
"""

import math
import sys
from collections.abc import Callable
from functools import cache
from itertools import pairwise
from pathlib import Path
from typing import Any

import msgpack
import numpy as np
import openmc
import openmc_tallies

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import material
from openmc_model import (
    CLAD_INNER,
    CLAD_OUTER,
    FUEL,
    GUIDE_INNER,
    GUIDE_OUTER,
    GUIDE_TUBES,
    HALF_GAP,
    INSTRUMENT_INNER,
    INSTRUMENT_OUTER,
    INSTRUMENT_TUBE,
    MODERATOR_DENSE,
    PITCH,
    PYREX,
    PYREX_OCTANTS,
    PYREX_RADII,
    SS304,
    ZIRCALOY,
    Frame,
    Layers,
    N,
    concentric,
    fuel_layers,
    octant_images,
    pin_lattice,
    pyrex_rod_layers,
    temperature_mode,
    with_carbon,
)
from openmc_resume import (
    OPTIONAL,
    Options,
    run_model,
    run_or_extract,
    set_statepoints,
)
from openmc_run import (
    box_settings,
    fuel_mean,
    mesh_values,
    write_k,
    write_pin_powers,
)
from openmc_run import read_run as common_read_run

# 表 P3-3 の個数密度 (/(barn cm))
# 3B の濃縮度 2.619% の燃料
FUEL_2619 = {
    "U234": 5.09503e-06,
    "U235": 6.06733e-04,
    "U236": 2.76809e-06,
    "U238": 2.22663e-02,
    "O16": 4.57617e-02,
}
# 3B の 1066 ppm の減速材（0.743 g/cc）
MODERATOR_1066 = {
    "H1": 4.96340e-02,
    "B10": 8.77976e-06,
    "B11": 3.53397e-05,
    "O16": 2.48170e-02,
}
INCONEL = {
    "Si28": 4.04885e-03,
    "Si29": 2.05685e-04,
    "Si30": 1.35748e-04,
    "Ti46": 2.12518e-04,
    "Ti47": 1.91652e-04,
    "Ti48": 1.89901e-03,
    "Ti49": 1.39360e-04,
    "Ti50": 1.33435e-04,
    "Cr50": 6.18222e-04,
    "Cr52": 1.19218e-02,
    "Cr53": 1.35184e-03,
    "Cr54": 3.36501e-04,
    "Fe54": 3.61353e-04,
    "Fe56": 5.67247e-03,
    "Fe57": 1.31002e-04,
    "Fe58": 1.74340e-05,
    "Ni58": 4.17608e-02,
    "Ni60": 1.60862e-02,
    "Ni61": 6.99255e-04,
    "Ni62": 2.22953e-03,
    "Ni64": 5.67796e-04,
}
# ノズルと炉心板の SS304 の部分（3A と 3B で同じ）。上部ノズル、下部ノズル、炉心板の順
STEEL_PARTS = {
    "Si28": (3.02920e-04, 4.41720e-04, 7.90985e-04),
    "Si29": (1.53886e-05, 2.24397e-05, 4.01826e-05),
    "Si30": (1.01561e-05, 1.48097e-05, 2.65197e-05),
    "P31": (1.34026e-05, 1.95438e-05, 3.49969e-05),
    "Cr50": (1.46468e-04, 2.13581e-04, 3.82458e-04),
    "Cr52": (2.82449e-03, 4.11869e-03, 7.37532e-03),
    "Cr53": (3.20275e-04, 4.67027e-04, 8.36302e-04),
    "Cr54": (7.97232e-05, 1.16253e-04, 2.08173e-04),
    "Mn55": (3.35836e-04, 4.89719e-04, 8.76936e-04),
    "Fe54": (6.60188e-04, 9.62690e-04, 1.72388e-03),
    "Fe56": (1.03635e-02, 1.51122e-02, 2.70613e-02),
    "Fe57": (2.39339e-04, 3.49006e-04, 6.24963e-04),
    "Fe58": (3.18517e-05, 4.64463e-05, 8.31710e-05),
    "Ni58": (1.01650e-03, 1.48226e-03, 2.65427e-03),
    "Ni60": (3.91552e-04, 5.70964e-04, 1.02242e-03),
    "Ni61": (1.70205e-05, 2.48194e-05, 4.44439e-05),
    "Ni62": (5.42688e-05, 7.91351e-05, 1.41707e-04),
    "Ni64": (1.38207e-05, 2.01534e-05, 3.60885e-05),
}
STEEL_CARBON = (6.14459e-05, 8.96008e-05, 1.60447e-04)
# ノズルと炉心板の水の部分（H-1, B-10, B-11, O-16）。問題ごとに、上部ノズル、下部ノズル、炉心板の順
WATER_PARTS = {
    "3A": (
        (4.01211e-02, 8.65222e-06, 3.48263e-05, 2.00606e-02),
        (3.57661e-02, 7.70514e-06, 3.10142e-05, 1.78830e-02),
        (2.48112e-02, 5.33040e-06, 2.14555e-05, 1.24056e-02),
    ),
    "3B": (
        (4.01305e-02, 7.09198e-06, 2.85461e-05, 2.00653e-02),
        (3.57744e-02, 6.32374e-06, 2.54539e-05, 1.78872e-02),
        (2.48171e-02, 4.40970e-06, 1.77496e-05, 1.24085e-02),
    ),
}
STRUCTURES = ("top nozzle", "bottom nozzle", "core plate")


def mixture(problem: str, structure: str) -> dict[str, float]:
    """ノズルと炉心板（`STRUCTURES` の 1 つ）の、SS304 と水を均質にした組成（表 P3-3）。"""
    k = STRUCTURES.index(structure)
    water = dict(zip(("H1", "B10", "B11", "O16"), WATER_PARTS[problem][k]))
    steel = {nuclide: parts[k] for nuclide, parts in STEEL_PARTS.items()}
    return with_carbon({**water, **steel}, STEEL_CARBON[k])


# 問題ごとの等温の温度 (K)、燃料、減速材（表 P3-1、P3-3）
TEMPERATURE = {"3A": 600.0, "3B": 565.0}
FUEL_OF = {"3A": FUEL, "3B": FUEL_2619}
MODERATOR_OF = {"3A": MODERATOR_DENSE, "3B": MODERATOR_1066}

# 軸方向の位置 (cm)。表 1、4、5
PLATE_BOTTOM = -5.0
NOZZLE_BOTTOM_TOP = 6.053
ROD_BOTTOM = 10.281
FUEL_BOTTOM = 11.951
FUEL_TOP = FUEL_BOTTOM + 365.76
PLENUM_TOP = FUEL_TOP + 16.0
ROD_TOP = PLENUM_TOP + 1.67
NOZZLE_TOP_BOTTOM = 397.51
NOZZLE_TOP_TOP = 406.337
PLATE_TOP = NOZZLE_TOP_TOP + 7.6
# 炉心板の外の水の厚さ（仕様書の問題 3 の「Geometry」）
WATER = 50.0

# 支持格子（表 3）: 物質の名前、質量 (g)、密度 (g/cc。2 節)、高さ (cm)、中心の位置 (cm)
GRIDS = (
    ("inconel", 1017.0, 8.19, 3.866, (13.884, 388.2)),
    ("zircaloy-4", 875.0, 6.56, 3.810, (75.2, 127.4, 179.6, 231.8, 284.0, 336.2)),
)


def frame_thickness(mass: float, density: float, height: float) -> float:
    """支持格子 1 つの質量を 289 のセルに等しく分けた面積 A を、ピッチ p の正方形の外周の帯 4pt − 4t² で表す厚さ t
    （openmc_model.py の GRID_FRAME と同じ置き方）。"""
    area = mass / density / height / (N * N)
    return (PITCH - math.sqrt(PITCH * PITCH - area)) / 2


# 3B の Pyrex の棒（表 6）: 毒物の下端と高さ、その上のプレナムの高さ、下の端栓の高さ
POISON_BOTTOM = 15.761
POISON_HEIGHT = 360.68
POISON_PLENUM = 22.2
POISON_PLUG = 2.54
# 3B のシンブルプラグ（表 9）: 半径、高さ、下端
PLUG_RADIUS = 0.538
PLUG_HEIGHT = 11.0
PLUG_BOTTOM = 383.31

# 図 5 の 16 本の Pyrex の配置。案内管のほかの 1/8 の位置は (3, 3) と (6, 0)
PYREX_16 = octant_images(PYREX_OCTANTS[16])

# 表 P3-5 の 49 層の上端の z (cm)。下端は燃料の下端
LEVEL_TOPS = (
    15.817,
    24.0281,
    32.2393,
    40.4504,
    48.6616,
    56.8727,
    65.0839,
    73.295,
    77.105,
    85.17,
    93.235,
    101.3,
    109.365,
    117.43,
    125.495,
    129.305,
    137.37,
    145.435,
    153.5,
    161.565,
    169.63,
    177.695,
    181.505,
    189.57,
    197.635,
    205.7,
    213.765,
    221.83,
    229.895,
    233.705,
    241.77,
    249.835,
    257.9,
    265.965,
    274.03,
    282.095,
    285.905,
    293.97,
    302.035,
    310.1,
    318.165,
    326.23,
    334.295,
    338.105,
    346.0262,
    353.9474,
    361.8686,
    369.7898,
    377.711,
)
LEVELS = (FUEL_BOTTOM, *LEVEL_TOPS)
if not math.isclose(LEVELS[-1], FUEL_TOP):
    raise SystemExit("the axial levels do not end at the top of the fuel")


class Materials:
    """問題の物質。"""

    def __init__(self, problem: str):
        t = TEMPERATURE[problem]
        self.fuel = material("fuel", FUEL_OF[problem], t)
        self.clad = material("zircaloy-4", ZIRCALOY, t)
        self.water = material("moderator", MODERATOR_OF[problem], t)
        self.inconel = material("inconel", INCONEL, t)
        self.steel = material("ss304", SS304, t)
        self.pyrex = material("pyrex", PYREX, t)
        self.structures = {s: material(s, mixture(problem, s), t) for s in STRUCTURES}
        for m in (self.water, *self.structures.values()):
            m.add_s_alpha_beta("c_H_in_H2O")


def frames(m: Materials) -> list[tuple[float, float, Frame]]:
    """支持格子の (下端, 上端, 枠) の並び。"""
    fills = {"inconel": m.inconel, "zircaloy-4": m.clad}
    out = []
    for name, mass, density, height, centres in GRIDS:
        frame = (fills[name], frame_thickness(mass, density, height))
        out.extend((c - height / 2, c + height / 2, frame) for c in centres)
    return out


def fuel_rod(m: Materials, z: float) -> Layers:
    """高さ `z` の燃料棒の層（燃料棒のない高さでは空）。"""
    if ROD_BOTTOM <= z < FUEL_BOTTOM or PLENUM_TOP <= z < ROD_TOP:
        return [(m.clad, CLAD_OUTER)]
    if FUEL_BOTTOM <= z < FUEL_TOP:
        return fuel_layers(m.fuel, m.clad)
    if FUEL_TOP <= z < PLENUM_TOP:
        return [(None, CLAD_INNER), (m.clad, CLAD_OUTER)]
    return []


def pyrex_rod(m: Materials, z: float) -> Layers:
    """高さ `z` の Pyrex の棒の層（内管、Pyrex、被覆管。表 6）。上部ノズルの下端で切る。"""
    clad_in, clad_out = PYREX_RADII[4:]
    poison_top = POISON_BOTTOM + POISON_HEIGHT
    if POISON_BOTTOM - POISON_PLUG <= z < POISON_BOTTOM:
        return [(m.steel, clad_out)]
    if POISON_BOTTOM <= z < poison_top:
        return pyrex_rod_layers(m.pyrex, m.steel)
    if poison_top <= z < min(poison_top + POISON_PLENUM, NOZZLE_TOP_BOTTOM):
        return [(None, clad_in), (m.steel, clad_out)]
    return []


def thimble_plug(m: Materials, z: float) -> Layers:
    """高さ `z` のシンブルプラグの層。"""
    if PLUG_BOTTOM <= z < PLUG_BOTTOM + PLUG_HEIGHT:
        return [(m.steel, PLUG_RADIUS)]
    return []


def boundaries(problem: str, m: Materials) -> list[float]:
    """集合体の部分（ノズルの間）で、中身が変わる z の並び（昇順）。"""
    zs = {NOZZLE_BOTTOM_TOP, ROD_BOTTOM, FUEL_BOTTOM, FUEL_TOP, PLENUM_TOP}
    zs |= {ROD_TOP, NOZZLE_TOP_BOTTOM}
    for bottom, top, _ in frames(m):
        zs |= {bottom, top}
    if problem == "3B":
        poison_top = POISON_BOTTOM + POISON_HEIGHT
        zs |= {POISON_BOTTOM - POISON_PLUG, POISON_BOTTOM, poison_top}
        zs |= {poison_top + POISON_PLENUM, PLUG_BOTTOM, PLUG_BOTTOM + PLUG_HEIGHT}
    zs = sorted(z for z in zs if NOZZLE_BOTTOM_TOP <= z <= NOZZLE_TOP_BOTTOM)
    return zs


def lattice_at(
    problem: str, m: Materials, z: float, universes_of: Callable
) -> openmc.RectLattice:
    """高さ `z`（ノズルの間）の 17×17 の格子。外周の水のギャップは格子の外（`outer`）。`universes_of(layers, frame)` は
    セルのユニバースを作る（同じ中身のセルを層の間で再利用する）。"""
    frame = None
    for bottom, top, f in frames(m):
        if bottom <= z < top:
            frame = f

    def universe(layers: tuple) -> openmc.Universe:
        return universes_of(layers, frame)

    def tube(inside: Layers, inner: float, outer: float) -> openmc.Universe:
        return universe((*inside, (m.water, inner), (m.clad, outer)))

    rod = universe(tuple(fuel_rod(m, z)))
    universes = [[rod] * N for _ in range(N)]
    for i, j in GUIDE_TUBES:
        inside: Layers = []
        if problem == "3B":
            inside = pyrex_rod(m, z) if (i, j) in PYREX_16 else thimble_plug(m, z)
        universes[i][j] = tube(inside, GUIDE_INNER, GUIDE_OUTER)
    i, j = INSTRUMENT_TUBE
    universes[i][j] = tube([], INSTRUMENT_INNER, INSTRUMENT_OUTER)
    return pin_lattice(universes, m.water)


def geometry_of(problem: str, m: Materials) -> tuple[openmc.Geometry, float]:
    """問題の形状と、外周の箱の半分の幅 (cm)。"""
    if PYREX_16 - set(GUIDE_TUBES):
        raise SystemExit("Pyrex rods outside the guide tubes")
    half = N * PITCH / 2 + HALF_GAP
    box = -openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    bottom = PLATE_BOTTOM - WATER
    top = PLATE_TOP + WATER
    # (上端, 中身) を下から
    zones: list[tuple[float, Any]] = [
        (PLATE_BOTTOM, m.water),
        (0.0, m.structures["core plate"]),
        (NOZZLE_BOTTOM_TOP, m.structures["bottom nozzle"]),
    ]
    zs = boundaries(problem, m)
    universes_of = cache(lambda layers, frame: concentric(list(layers), m.water, frame))
    for lower, upper in pairwise(zs):
        z = (lower + upper) / 2
        zones.append((upper, lattice_at(problem, m, z, universes_of)))
    zones += [
        (NOZZLE_TOP_TOP, m.structures["top nozzle"]),
        (PLATE_TOP, m.structures["core plate"]),
        (top, m.water),
    ]
    cells = []
    lower_plane = openmc.ZPlane(z0=bottom, boundary_type="vacuum")
    for upper, fill in zones:
        kind = "vacuum" if upper == top else "transmission"
        upper_plane = openmc.ZPlane(z0=upper, boundary_type=kind)
        cells.append(openmc.Cell(fill=fill, region=box & +lower_plane & -upper_plane))
        lower_plane = upper_plane
    return openmc.Geometry(cells), half


def pin_mesh(z_grid: list[float]) -> openmc.RectilinearMesh:
    """ピンのセルに一致し、z を `z_grid` で分けるメッシュ。"""
    edges = [(k - N / 2) * PITCH for k in range(N + 1)]
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = edges
    mesh.y_grid = edges
    mesh.z_grid = z_grid
    return mesh


def axial_mesh(half: float) -> openmc.RectilinearMesh:
    """集合体の全体を表 P3-5 の 49 層に分けるメッシュ。"""
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [-half, half]
    mesh.y_grid = [-half, half]
    mesh.z_grid = list(LEVELS)
    return mesh


def read_run(source: Path) -> dict[str, Any]:
    return common_read_run(
        source, list(TEMPERATURE), OPTIONAL | openmc_tallies.OPTIONAL
    )


def solve(source: Path, work: Path, options: Options) -> None:
    run = read_run(source)
    problem = run["problem"]
    work.mkdir(parents=True, exist_ok=True)
    m = Materials(problem)
    geometry, half = geometry_of(problem, m)
    materials = openmc.Materials(geometry.get_all_materials().values())
    materials.cross_sections = str(run["cross-sections"])

    settings = box_settings(
        run["particles"],
        run["batches"],
        run["inactive"],
        run["seed"],
        (-half, -half, FUEL_BOTTOM),
        (half, half, FUEL_TOP),
    )
    # 軸方向の源の収束を世代ごとに見るためのエントロピー（燃料の範囲を 4×4×48 に分ける）
    entropy = openmc.RegularMesh()
    entropy.lower_left = (-half, -half, FUEL_BOTTOM)
    entropy.upper_right = (half, half, FUEL_TOP)
    entropy.dimension = (4, 4, 48)
    settings.entropy_mesh = entropy
    set_statepoints(settings, run)
    settings.temperature = temperature_mode(TEMPERATURE[problem])

    # 層ごとのピン、全長で積分したピン、層ごとの集合体全体。積分した値も別の tally で数えるのは、ピンや層の間の世代
    # ごとの相関を含めた統計誤差を、OpenMC に世代の値から出させるため（層ごとのピンの誤差の二乗和は相関を落とす）
    meshes = {
        "fission": pin_mesh(list(LEVELS)),
        "pins": pin_mesh([FUEL_BOTTOM, FUEL_TOP]),
        "axial": axial_mesh(half),
    }
    tallies = openmc.Tallies()
    for name, mesh in meshes.items():
        tally = openmc.Tally(name=name)
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["fission"]
        tallies.append(tally)
    if openmc_tallies.enabled(run):
        # 層は、燃料の下と上の構造（炉心板、ノズル、端栓、プレナム）の境と、表 P3-5 の 49 層
        below = [PLATE_BOTTOM - WATER, PLATE_BOTTOM, NOZZLE_BOTTOM_TOP, ROD_BOTTOM]
        above = [PLENUM_TOP, ROD_TOP, NOZZLE_TOP_BOTTOM, NOZZLE_TOP_TOP, PLATE_TOP]
        tallies.extend(
            openmc_tallies.assembly_3d(
                list(materials),
                [m.fuel],
                # 上の tally と別のメッシュにする（同じメッシュのフィルタを OpenMC がまとめるときに ID の警告が出る）
                pin_mesh([FUEL_BOTTOM, FUEL_TOP]),
                half,
                HALF_GAP,
                [*below, *LEVELS, *above, PLATE_TOP + WATER],
            )
        )
        openmc_tallies.quiet(settings)
    model = openmc.Model(geometry, materials, settings, tallies)
    run_model(model, run, work, options)


def extract(source: Path, work: Path, results: Path) -> None:
    run = read_run(source)
    problem = run["problem"]
    results.mkdir(parents=True, exist_ok=True)
    levels = len(LEVELS) - 1
    with openmc.StatePoint(work / f"statepoint.{run['batches']}.h5") as sp:
        k = sp.keff
        rate, std = mesh_values(sp, "fission", (levels, N, N))
        pins, pins_std = mesh_values(sp, "pins", (1, N, N))
        axial, axial_std = mesh_values(sp, "axial", (levels, 1, 1))
    heights = np.diff(LEVELS)
    pins, pins_std = pins[0], pins_std[0]
    height = LEVELS[-1] - LEVELS[0]

    # 全長で積分したピン: 燃料のピンの平均で 1
    mean = fuel_mean(pins)
    write_pin_powers(results / "pins.msgpack", pins / mean)
    write_pin_powers(results / "pins-std.msgpack", pins_std / mean)

    # 線密度: 燃料の全長 H での平均で 1（層ごとのピンは、燃料のピンと全長での平均で 1）
    axial_scale = height / axial.sum() / heights
    layer_scale = height / fuel_mean(rate.sum(axis=0)) / heights[:, None, None]
    layers = {
        "format": "xlc-axial-pin-powers",
        "bounds": list(LEVELS),
        "axial": (axial.ravel() * axial_scale).tolist(),
        "axial-std": (axial_std.ravel() * axial_scale).tolist(),
        "layers": (rate * layer_scale).tolist(),
        "layers-std": (std * layer_scale).tolist(),
    }
    (results / "layers.msgpack").write_bytes(msgpack.packb(layers))
    write_k(results / "k.json", problem, k)


if __name__ == "__main__":
    run_or_extract(solve, extract, __doc__)
