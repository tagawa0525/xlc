"""VERA Problem 1B〜1D（ピンセル）と 2A〜2K、2Q（17×17 集合体）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py model <入力> <モデルのディレクトリ>
        python3 openmc_model.py run <入力> <モデルのディレクトリ> <作業ディレクトリ>
        python3 openmc_model.py extract <入力> <作業ディレクトリ> <結果のディレクトリ>

<入力> は p1/openmc-1B.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 問題（`problem`）、
openmc_library.py で作った `cross-sections`（cross_sections.xml）、1 世代の粒子数（`particles`）、世代数（`batches`）、
捨てる世代数（`inactive`）、乱数の種（`seed`）、スレッド数（`threads`）。パスは入力のファイルからの相対。
標準のタリーの組（openmc_tallies.py、ADR-0152）を入れ、`standard-tallies = false` と書いたときだけ入れない。
model は問題の物質、形状、設定、tally を <モデルのディレクトリ>/model.xml に書く。model.xml は断面積のライブラリの
パスもスレッド数も含まない（このスクリプトを直しても、モデルが同じなら model.xml はバイト単位で変わらず、OpenMC を
走り直さずに済む。ADR-0136）。run はそのモデルを `cross-sections` のライブラリとスレッド数で解き、作業ディレクトリに
statepoint を残す。extract は statepoint から、k を k.json に、下のピン出力と実効断面積を結果のディレクトリに書く
（OpenMC を解き直さずに抽出を変えられるよう分ける。ADR-0123）。

問題は 1B, 1C, 1D, 2A, 2B, 2C, 2D, 2E, 2E-empty（2E の Pyrex の棒を抜いたもの。吸収体と減速材の密度の効果を分けるため）, 2F,
2G, 2H, 2I, 2J, 2K, 2Q。
組成は仕様書（CASL-U-2012-0131-004）の表 P1-4, P2-4 の個数密度で、Hf（1e-6 /(barn cm) 未満）と
ギャップの He は ADR-0054 のとおり除く（ギャップは真空。Pyrex の棒と制御棒の中の He も同じ。ADR-0120）。燃料以外は、
2A と 2Q は 565 K、ほかは 600 K。2Q の支持格子は、各セルの外周の Zircaloy-4 の枠にする（仕様書の表 3 と問題 2 の
「Geometry」）。2 で始まる問題は lattice と同じ形式のピンごとの核分裂率を pins.msgpack に書く。
標準のタリーの組を入れたときは、その WIMS 69 群（xlc のライブラリの群構造）の燃料の核種ごとの実効断面積
（燃料の領域の反応率 / (個数密度 × 中性子束)）を、ピンセル（1 で始まる問題）では fuel-xs.msgpack に、集合体（2 で始まる
問題）ではピンごとに pin-xs.msgpack に書く（自己遮蔽の検証用。燃料のないセルは 0）。
"""

import os
import sys
from pathlib import Path
from typing import Any, NamedTuple

import msgpack
import numpy as np
import openmc
import openmc_tallies

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import STORED, eigenvalue_settings, material
from openmc_run import (
    fuel_mean,
    mesh_values,
    model_run_extract,
    write_k,
    write_pin_powers,
)
from openmc_run import read_run as common_read_run

SCORES = ["total", "elastic", "fission"]

PITCH = 1.26
N = 17
# 集合体の外周の水ギャップ（片側）
HALF_GAP = 0.04
PELLET = 0.4096
CLAD_INNER = 0.418
CLAD_OUTER = 0.475
GUIDE_INNER = 0.561
GUIDE_OUTER = 0.602
INSTRUMENT_INNER = 0.559
INSTRUMENT_OUTER = 0.605
# ライブラリにある温度そのものの問題（1B〜1D、2B〜2K）の被覆管と減速材の温度 (K)。温度のモードの選択に使う
COLD = 600.0
FUEL_TEMPERATURE = {
    "1B": 600.0,
    "1C": 900.0,
    "1D": 1200.0,
    "2A": 565.0,
    "2B": 600.0,
    "2C": 900.0,
    "2D": 1200.0,
    "2E": 600.0,
    "2E-empty": 600.0,
    "2F": 600.0,
    "2G": 600.0,
    "2H": 600.0,
    "2I": 600.0,
    "2J": 600.0,
    "2K": 600.0,
    "2Q": 565.0,
}
# 被覆管と減速材の温度 (K)。仕様書の表 P2-1（被覆管は減速材の温度）
MODERATOR_TEMPERATURE = {
    "1B": 600.0,
    "1C": 600.0,
    "1D": 600.0,
    "2A": 565.0,
    "2B": 600.0,
    "2C": 600.0,
    "2D": 600.0,
    "2E": 600.0,
    "2E-empty": 600.0,
    "2F": 600.0,
    "2G": 600.0,
    "2H": 600.0,
    "2I": 600.0,
    "2J": 600.0,
    "2K": 600.0,
    "2Q": 565.0,
}
# 支持格子の枠の厚さ (cm)。中間の支持格子 1 つの Zircaloy-4 875 g（仕様書の表 3、密度 6.56 g/cc、高さ 3.810 cm）を
# 289 のセルに等しく分けた面積 0.121138 cm² を、ピッチ p の正方形の外周の帯 4pt − 4t² で表す厚さ t
GRID_FRAME = {"2Q": 0.024512}

# 個数密度 (/(barn cm))
FUEL = {
    "U234": 6.11864e-06,
    "U235": 7.18132e-04,
    "U236": 3.29861e-06,
    "U238": 2.21546e-02,
    "O16": 4.57642e-02,
}
# 2K の濃縮度 3.6% の燃料（表 P2-4）
FUEL_36 = {
    "U234": 7.21203e-06,
    "U235": 8.33952e-04,
    "U236": 3.82913e-06,
    "U238": 2.20384e-02,
    "O16": 4.57669e-02,
}
ZIRCALOY = {
    "Cr50": 3.30121e-06,
    "Cr52": 6.36606e-05,
    "Cr53": 7.21860e-06,
    "Cr54": 1.79686e-06,
    "Fe54": 8.68307e-06,
    "Fe56": 1.36306e-04,
    "Fe57": 3.14789e-06,
    "Fe58": 4.18926e-07,
    "Zr90": 2.18865e-02,
    "Zr91": 4.77292e-03,
    "Zr92": 7.29551e-03,
    "Zr94": 7.39335e-03,
    "Zr96": 1.19110e-03,
    "Sn112": 4.68066e-06,
    "Sn114": 3.18478e-06,
    "Sn115": 1.64064e-06,
    "Sn116": 7.01616e-05,
    "Sn117": 3.70592e-05,
    "Sn118": 1.16872e-04,
    "Sn119": 4.14504e-05,
    "Sn120": 1.57212e-04,
    "Sn122": 2.23417e-05,
    "Sn124": 2.79392e-05,
}
# 1300 ppm のホウ素。0.661 g/cc（1B〜1D、2B〜2D）と 0.743 g/cc（2A、2E〜2K、2Q）。仕様書の表 P2-1、P2-4
MODERATOR = {
    "O16": 2.20729e-02,
    "H1": 4.41459e-02,
    "B10": 9.52537e-06,
    "B11": 3.83408e-05,
}
MODERATOR_DENSE = {
    "O16": 2.48112e-02,
    "H1": 4.96224e-02,
    "B10": 1.07070e-05,
    "B11": 4.30971e-05,
}
PYREX = {
    "B10": 9.63266e-04,
    "B11": 3.90172e-03,
    "O16": 4.67761e-02,
    "Si28": 1.81980e-02,
    "Si29": 9.24474e-04,
    "Si30": 6.10133e-04,
}
# 天然の C は天然の存在比 C-12 98.93%、C-13 1.07% で分ける（ADR-0079）
SS304 = {
    "C12": 3.17461e-04,
    "C13": 3.43358e-06,
    "Si28": 1.58197e-03,
    "Si29": 8.03653e-05,
    "Si30": 5.30394e-05,
    "P31": 6.99938e-05,
    "Cr50": 7.64915e-04,
    "Cr52": 1.47506e-02,
    "Cr53": 1.67260e-03,
    "Cr54": 4.16346e-04,
    "Mn55": 1.75387e-03,
    "Fe54": 3.44776e-03,
    "Fe56": 5.41225e-02,
    "Fe57": 1.24992e-03,
    "Fe58": 1.66342e-04,
    "Ni58": 5.30854e-03,
    "Ni60": 2.04484e-03,
    "Ni61": 8.88879e-05,
    "Ni62": 2.83413e-04,
    "Ni64": 7.21770e-05,
}
AIC = {
    "Ag107": 2.36159e-02,
    "Ag109": 2.19403e-02,
    "Cd106": 3.41523e-05,
    "Cd108": 2.43165e-05,
    "Cd110": 3.41250e-04,
    "Cd111": 3.49720e-04,
    "Cd112": 6.59276e-04,
    "Cd113": 3.33873e-04,
    "Cd114": 7.84957e-04,
    "Cd116": 2.04641e-04,
    "In113": 3.44262e-04,
    "In115": 7.68050e-03,
}
# 天然の C は SS304 と同じく C-12 と C-13 に分ける（ADR-0079）
B4C = {
    "B10": 1.52689e-02,
    "B11": 6.14591e-02,
    "C12": 1.89768e-02,
    "C13": 2.05247e-04,
}
# 問題ごとの減速材。表にない問題は解かない（密度を取り違えないように）
MODERATOR_OF = {
    "1B": MODERATOR,
    "1C": MODERATOR,
    "1D": MODERATOR,
    "2A": MODERATOR_DENSE,
    "2B": MODERATOR,
    "2C": MODERATOR,
    "2D": MODERATOR,
    "2E": MODERATOR_DENSE,
    "2E-empty": MODERATOR_DENSE,
    "2F": MODERATOR_DENSE,
    "2G": MODERATOR_DENSE,
    "2H": MODERATOR_DENSE,
    "2I": MODERATOR_DENSE,
    "2J": MODERATOR_DENSE,
    "2K": MODERATOR_DENSE,
    "2Q": MODERATOR_DENSE,
}
# Pyrex の棒（仕様書の表 6）: 内管の内と外、Pyrex の内と外、被覆管の内と外の半径。内管の中と間は真空
PYREX_RADII = (0.214, 0.231, 0.241, 0.427, 0.437, 0.484)
# 制御棒（仕様書の表 8）: SS304 の被覆管の内と外の半径と、吸収体の半径。吸収体と被覆管の間は真空
ROD_CLAD = (0.386, 0.484)
# 計装管のシンブル（仕様書の表 10）: SS304 の内と外の半径。中は真空
THIMBLE_RADII = (0.258, 0.382)

# 案内管の位置（行, 列。0 始まり）。計装管は中央
GUIDE_TUBES = [
    (2, 5),
    (2, 8),
    (2, 11),
    (3, 3),
    (3, 13),
    (5, 2),
    (5, 5),
    (5, 8),
    (5, 11),
    (5, 14),
    (8, 2),
    (8, 5),
    (8, 11),
    (8, 14),
    (11, 2),
    (11, 5),
    (11, 8),
    (11, 11),
    (11, 14),
    (13, 3),
    (13, 13),
    (14, 5),
    (14, 8),
    (14, 11),
]
INSTRUMENT_TUBE = (8, 8)
# 図 5 の Pyrex の棒の配置（本数ごとの、1/8 の中の位置。中央からの行と列の数）
PYREX_OCTANTS = {
    8: {(3, 0), (5, 5)},
    12: {(3, 0), (6, 3)},
    16: {(3, 0), (5, 5), (6, 3)},
    20: {(3, 0), (5, 5), (6, 0), (6, 3)},
    24: {(3, 0), (3, 3), (5, 5), (6, 0), (6, 3)},
}
# 天然の C の同位体の存在比（ADR-0079）
CARBON = {"C12": 0.9893, "C13": 0.0107}
# 2K の濃縮度 3.1% のピンの、1/8 の中の位置（中央からの行と列の数。仕様書の図 P2-2）。ほかの燃料ピンは 3.6%
LOW_ENRICHMENT_2K = {(1, 0), (7, 7), (8, 6), (8, 7), (8, 8)}


def octant_images(octant: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """1/8 の中の位置（中央からの行と列の数）を、8 つの対称操作で移した集合体の位置（行, 列。0 始まり）にする。"""
    centre = N // 2
    for a, b in octant:
        if not 0 <= b <= a <= centre:
            raise SystemExit(f"{(a, b)} is not in the octant")
    return {
        (centre + si * x, centre + sj * y)
        for a, b in octant
        for x, y in ((a, b), (b, a))
        for si in (1, -1)
        for sj in (1, -1)
    }


def with_carbon(densities: dict[str, float], carbon: float) -> dict[str, float]:
    """仕様書の天然の C（6000）の個数密度 `carbon` を、C-12 と C-13 に分けて足す（ADR-0079）。"""
    return {**densities, **{n: carbon * f for n, f in CARBON.items()}}


Layers = list[tuple[openmc.Material | None, float]]
# 支持格子の枠の物質と厚さ (cm)
Frame = tuple[openmc.Material, float]


def outside(
    region: openmc.Region, water: openmc.Material, frame: Frame | None
) -> list[openmc.Cell]:
    """セルの最も外の円の外（`region`）の減速材と、`frame` があればセルの外周の正方形の帯の枠のセル。枠の外側は
    格子のセルの境界で切り取られる。"""
    if frame is None:
        return [openmc.Cell(fill=water, region=region)]
    material, thickness = frame
    box = openmc.model.RectangularPrism(PITCH - 2 * thickness, PITCH - 2 * thickness)
    return [
        openmc.Cell(fill=water, region=region & -box),
        openmc.Cell(fill=material, region=+box),
    ]


def concentric(
    layers: Layers, water: openmc.Material, frame: Frame | None = None
) -> openmc.Universe:
    """中心から外へ (物質, 外側の半径) の並び（None は真空）の同心円筒と、その外の減速材（`frame` があれば支持格子の
    枠も）のセル。`layers` が空なら減速材だけ（枠は置けない）。"""
    if not layers:
        if frame is not None:
            raise SystemExit("a spacer grid frame around an empty cell")
        return openmc.Universe(cells=[openmc.Cell(fill=water)])
    cells = []
    surface = None
    for fill, radius in layers:
        boundary = openmc.ZCylinder(r=radius)
        region = -boundary if surface is None else +surface & -boundary
        cells.append(openmc.Cell(fill=fill, region=region))
        surface = boundary
    cells.extend(outside(+surface, water, frame))
    return openmc.Universe(cells=cells)


def tube(
    layers: Layers,
    inner: float,
    outer: float,
    clad: openmc.Material,
    water: openmc.Material,
    frame: Frame | None = None,
) -> openmc.Universe:
    """管（案内管・計装管、内と外の半径 `inner`、`outer`）のセル。`layers` は管に入れたもので、中心から外へ
    (物質, 外側の半径) の並び（None は真空）。その外から管の内側までと、管の外は減速材。`frame` は支持格子の枠。"""
    return concentric([*layers, (water, inner), (clad, outer)], water, frame)


def fuel_layers(fuel: openmc.Material, clad: openmc.Material) -> Layers:
    """燃料のピンの層（ペレット、真空のギャップ、被覆管）。"""
    return [(fuel, PELLET), (None, CLAD_INNER), (clad, CLAD_OUTER)]


def pin_lattice(
    universes: list[list[openmc.Universe]], water: openmc.Material
) -> openmc.RectLattice:
    """17×17 のセルの格子（中心が原点）。格子の外（集合体の外周の水のギャップ）は減速材。"""
    lattice = openmc.RectLattice()
    lattice.pitch = (PITCH, PITCH)
    lattice.lower_left = (-N * PITCH / 2, -N * PITCH / 2)
    lattice.universes = universes
    lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=water)])
    return lattice


def pyrex_layers(temperature: float) -> Layers:
    """Pyrex の棒（内管、Pyrex、被覆管。内管の中と間は真空）の層。`temperature` は物質の温度 (K)。"""
    pyrex = material("pyrex", PYREX, temperature)
    steel = material("ss304", SS304, temperature)
    return pyrex_rod_layers(pyrex, steel)


def pyrex_rod_layers(pyrex: openmc.Material, steel: openmc.Material) -> Layers:
    """Pyrex の棒の毒物の部分の層（物質 `pyrex` と、内管と被覆管の `steel`）。"""
    fills = [None, steel, None, pyrex, None, steel]
    return list(zip(fills, PYREX_RADII))


def rod_layers(
    name: str, nuclides: dict[str, float], radius: float, temperature: float
) -> Layers:
    """制御棒（半径 `radius` の吸収体と SS304 の被覆管。間は真空）の層。`temperature` は物質の温度 (K)。"""
    absorber = material(name, nuclides, temperature)
    steel = material("ss304", SS304, temperature)
    clad_in, clad_out = ROD_CLAD
    return [(absorber, radius), (None, clad_in), (steel, clad_out)]


def thimble_layers(temperature: float) -> Layers:
    """計装管のシンブル（中が真空の SS304 の管）の層。`temperature` は物質の温度 (K)。"""
    inner, outer = THIMBLE_RADII
    return [(None, inner), (material("ss304", SS304, temperature), outer)]


# 案内管に入れるものの層を作る関数（引数は物質の温度。減速材の温度にする）と位置。Pyrex の 12 本は仕様書の図 5 の配置
INSERTS = {
    "2E": (pyrex_layers, sorted(octant_images(PYREX_OCTANTS[12]))),
    # 24 本はすべての案内管（図 5）
    "2F": (pyrex_layers, GUIDE_TUBES),
    "2G": (lambda t: rod_layers("aic", AIC, 0.382, t), GUIDE_TUBES),
    "2H": (lambda t: rod_layers("b4c", B4C, 0.373, t), GUIDE_TUBES),
    "2J": (pyrex_layers, GUIDE_TUBES),
    "2K": (pyrex_layers, GUIDE_TUBES),
}
# 計装管に入れるものの層を作る関数
INSTRUMENT_INSERTS = {"2I": thimble_layers, "2J": thimble_layers}
# 2 種目の燃料（個数密度）と、元の燃料（FUEL）のまま残すピンの位置
ZONED = {"2K": (FUEL_36, octant_images(LOW_ENRICHMENT_2K))}


def geometry_of(
    problem: str,
    fuel: openmc.Material,
    clad: openmc.Material,
    water: openmc.Material,
    gap: openmc.Material | None = None,
    insert: openmc.Universe | None = None,
    instrument_insert: openmc.Universe | None = None,
    zoned_fuel: openmc.Material | None = None,
    frame: Frame | None = None,
) -> tuple[openmc.Geometry, float]:
    """問題の形状と、外周の箱の半分の幅 (cm)。外周は鏡面反射。`gap` はペレットと被覆管の間（None なら真空）。
    `insert` は案内管に入れるもののセルで、INSERTS にある問題で要る。`instrument_insert` は計装管に入れるものを
    入れた計装管のセルで、INSTRUMENT_INSERTS にある問題で要る。`zoned_fuel` は 2 種目の燃料で、ZONED にある問題で要る。
    `frame` は、燃料のピン、案内管、計装管のセルに付ける支持格子の枠で、挿入物のある問題では使えない。"""
    pellet, inner, outer = (
        openmc.ZCylinder(r=r) for r in (PELLET, CLAD_INNER, CLAD_OUTER)
    )

    def rod_of(pellet_fill: openmc.Material) -> openmc.Universe:
        return openmc.Universe(
            cells=[
                openmc.Cell(fill=pellet_fill, region=-pellet),
                openmc.Cell(fill=gap, region=+pellet & -inner),
                openmc.Cell(fill=clad, region=+inner & -outer),
                *outside(+outer, water, frame),
            ]
        )

    rod = rod_of(fuel)

    if frame is not None and (insert is not None or instrument_insert is not None):
        raise SystemExit("a spacer grid frame is not put on the cells of inserts")
    if problem.startswith("1"):
        half = PITCH / 2
        root = openmc.Cell(fill=rod)
    else:
        guide = tube([], GUIDE_INNER, GUIDE_OUTER, clad, water, frame)
        if problem in INSTRUMENT_INSERTS:
            if instrument_insert is None:
                raise SystemExit(
                    f"{problem} needs the cell of its instrument tube insert"
                )
            instrument = instrument_insert
        else:
            instrument = tube(
                [], INSTRUMENT_INNER, INSTRUMENT_OUTER, clad, water, frame
            )
        universes = [[rod] * N for _ in range(N)]
        if problem in ZONED:
            if zoned_fuel is None:
                raise SystemExit(f"{problem} needs its second fuel")
            zoned = rod_of(zoned_fuel)
            kept = ZONED[problem][1]
            tubes = {*GUIDE_TUBES, INSTRUMENT_TUBE} & kept
            if tubes:
                raise SystemExit(f"fuel kept at tube positions {sorted(tubes)}")
            for i in range(N):
                for j in range(N):
                    if (i, j) not in kept:
                        universes[i][j] = zoned
        for i, j in GUIDE_TUBES:
            universes[i][j] = guide
        universes[INSTRUMENT_TUBE[0]][INSTRUMENT_TUBE[1]] = instrument
        if problem in INSERTS:
            if insert is None:
                raise SystemExit(f"{problem} needs the cell of its guide tube inserts")
            for i, j in INSERTS[problem][1]:
                if (i, j) not in GUIDE_TUBES:
                    raise SystemExit(f"insert at {(i, j)} is not in a guide tube")
                universes[i][j] = insert
        lattice = pin_lattice(universes, water)
        half = N * PITCH / 2 + HALF_GAP
        root = openmc.Cell(fill=lattice)
    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    root.region = -box
    return openmc.Geometry([root]), half


def pin_mesh() -> openmc.RegularMesh:
    """集合体のピンのセルと一致するメッシュ。"""
    mesh = openmc.RegularMesh()
    mesh.dimension = (N, N)
    mesh.lower_left = (-N * PITCH / 2, -N * PITCH / 2)
    mesh.upper_right = (N * PITCH / 2, N * PITCH / 2)
    return mesh


def pin_rows(sp: openmc.StatePoint, value: str = "mean"):
    """名前が fission の tally（`pin_mesh` ごと）の `value`（"mean" か "std_dev"）を、lattice と同じく上の行から並べた
    17 × 17 にする。"""
    mean, std = mesh_values(sp, "fission", (N, N))
    return {"mean": mean, "std_dev": std}[value]


def write_pins(sp: openmc.StatePoint, path: Path) -> None:
    """名前が fission の tally を、lattice と同じ形式のピンの核分裂率（平均 1）にして書く。"""
    rates = pin_rows(sp)
    write_pin_powers(path, rates / fuel_mean(rates))


class Built(NamedTuple):
    """問題の物質と形状。`fuel` は（2K では 3.1% の）燃料、`materials` は形状に使うすべての物質。"""

    fuel: openmc.Material
    materials: openmc.Materials
    geometry: openmc.Geometry
    half: float


def build(problem: str) -> Built:
    """問題の物質と形状を作る（外周の箱の半分の幅 `half` とともに）。"""
    fuel = material("fuel", FUEL, FUEL_TEMPERATURE[problem])
    clad = material("zircaloy-4", ZIRCALOY, MODERATOR_TEMPERATURE[problem])
    water = material("moderator", MODERATOR_OF[problem], MODERATOR_TEMPERATURE[problem])
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])

    def add_materials(layers: Layers) -> None:
        for fill, _ in layers:
            if fill is not None and all(fill is not m for m in materials):
                materials.append(fill)

    zoned_fuel = None
    if problem in ZONED:
        zoned_fuel = material(
            "fuel-zoned", ZONED[problem][0], FUEL_TEMPERATURE[problem]
        )
        materials.append(zoned_fuel)

    insert = None
    if problem in INSERTS:
        layers = INSERTS[problem][0](MODERATOR_TEMPERATURE[problem])
        add_materials(layers)
        insert = tube(layers, GUIDE_INNER, GUIDE_OUTER, clad, water)
    instrument_insert = None
    if problem in INSTRUMENT_INSERTS:
        layers = INSTRUMENT_INSERTS[problem](MODERATOR_TEMPERATURE[problem])
        add_materials(layers)
        instrument_insert = tube(
            layers, INSTRUMENT_INNER, INSTRUMENT_OUTER, clad, water
        )

    # 支持格子は被覆管と同じ組成（仕様書の表 P2-4 の Cladding and Grid）
    frame = (clad, GRID_FRAME[problem]) if problem in GRID_FRAME else None
    geometry, half = geometry_of(
        problem,
        fuel,
        clad,
        water,
        insert=insert,
        instrument_insert=instrument_insert,
        zoned_fuel=zoned_fuel,
        frame=frame,
    )
    return Built(fuel, materials, geometry, half)


def temperature_settings(problem: str) -> dict:
    """問題の温度の扱い（`openmc.Settings.temperature`）。"""
    # OpenMC は温度を kT (eV) で持ち、K に直すと丸めでわずかにずれるので、1 K の幅を許す。600 K の問題は、
    # ライブラリにあるその温度の表だけを使う。565 K の問題は温度の内挿のモードで解き（ADR-0080）、H in H2O の
    # S(α,β) を 550 K と 573.6 K の表から kT の比で確率的に選ぶ。ほかの核種は 565 K の表を持つので、内挿の
    # モードでも 565 K の表を使う（OpenMC 0.15.2 の src/nuclide.cpp。#0017）
    return temperature_mode(MODERATOR_TEMPERATURE[problem])


def temperature_mode(temperature: float) -> dict:
    """被覆管と減速材の温度 `temperature` (K) での温度の扱い（`temperature_settings` の説明のとおり）。"""
    method = "nearest" if temperature == COLD else "interpolation"
    return {"method": method, "tolerance": 1.0}


def read_run(source: Path) -> dict[str, Any]:
    """入力の TOML を読み、問題を確かめ、パスを絶対パスにする。"""
    run = common_read_run(source, list(FUEL_TEMPERATURE), openmc_tallies.OPTIONAL)
    problem = run["problem"]
    if any(problem not in table for table in (MODERATOR_TEMPERATURE, MODERATOR_OF)):
        raise SystemExit(f"{problem} lacks its moderator in the tables")
    return run


def fuels_of(materials: openmc.Materials) -> list[openmc.Material]:
    """燃料の物質（2K は 2 種）。"""
    return [m for m in materials if m.name in ("fuel", "fuel-zoned")]


def standard_tallies(problem: str, materials: openmc.Materials) -> list[openmc.Tally]:
    """問題の種類（ピンセルか集合体）の標準のタリーの組（ADR-0152）。"""
    if problem.startswith("1"):
        return openmc_tallies.pincell(list(materials))
    half = N * PITCH / 2 + HALF_GAP
    return openmc_tallies.assembly(
        list(materials), fuels_of(materials), pin_mesh(), half, HALF_GAP
    )


def without_outside(groups: np.ndarray) -> np.ndarray:
    """標準のタリーの組の 71 のビン（昇順、範囲の外の 2 つを含む）を、xlc と同じく高いエネルギーから並べた 69 群にする。"""
    return groups[..., 1:-1][..., ::-1]


def fuel_xs(sp: openmc.StatePoint, fuel: openmc.Material) -> dict:
    """ピンセルの燃料の核種ごとの実効断面積（標準のタリーの組の `nuclides` と `material groups` から）。"""
    bins = {"filters": [openmc.MaterialFilter], "filter_bins": [(fuel.id,)]}
    flux = sp.get_tally(name=openmc_tallies.name("material groups"))
    phi = flux.get_values(scores=["flux"], **bins).ravel()
    rates = sp.get_tally(name=openmc_tallies.name("nuclides"))
    densities = fuel.get_nuclide_atom_densities()
    xs = {}
    for nuclide, density in densities.items():
        xs[nuclide] = {}
        for score in SCORES:
            rate = rates.get_values(
                scores=[STORED.get(score, score)], nuclides=[nuclide], **bins
            ).ravel()
            # 範囲の外のビン（捨てる）には中性子束が 0 のものがある
            value = np.divide(
                rate, density * phi, out=np.zeros_like(rate), where=phi > 0
            )
            xs[nuclide][score] = without_outside(value).tolist()
    return xs


def pin_xs(sp: openmc.StatePoint, fuels: list[openmc.Material]) -> dict:
    """集合体のピンごとの燃料の核種の実効断面積（標準のタリーの組の `pin nuclides` から）。

    1 つのピンには 1 種の燃料しかないので、燃料の物質についての和で割る（反応率の和 / Σ 個数密度 × 中性子束）。
    燃料のないセルは 0。"""

    # [メッシュ（x が速く、y は下から）][燃料][ビン] を [行（上から）][列][燃料][ビン] にする
    def pins(values: np.ndarray) -> np.ndarray:
        return values.reshape(N, N, len(fuels), -1)[::-1]

    phi = pins(
        sp.get_tally(name=openmc_tallies.name("pin nuclides flux")).get_values(
            scores=["flux"]
        )
    )
    rates = sp.get_tally(name=openmc_tallies.name("pin nuclides"))
    xs = {}
    for nuclide in FUEL:
        density = np.array([f.get_nuclide_atom_densities()[nuclide] for f in fuels])[
            :, None
        ]
        reacting = (density * phi).sum(axis=2)
        xs[nuclide] = {}
        for score in SCORES:
            rate = pins(
                rates.get_values(scores=[STORED.get(score, score)], nuclides=[nuclide])
            ).sum(axis=2)
            value = np.divide(
                rate, reacting, out=np.zeros_like(rate), where=reacting > 0
            )
            xs[nuclide][score] = without_outside(value).tolist()
    return xs


def model(source: Path, directory: Path) -> None:
    """問題の model.xml を `directory` に書く。"""
    run = read_run(source)
    problem = run["problem"]
    _, materials, geometry, half = build(problem)

    settings = eigenvalue_settings(
        run["particles"], run["batches"], run["inactive"], run["seed"], half
    )
    settings.temperature = temperature_settings(problem)

    tallies = openmc.Tallies()
    if problem.startswith("2"):
        mesh = pin_mesh()
        tally = openmc.Tally(name="fission")
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["fission"]
        tallies.append(tally)
    if openmc_tallies.enabled(run):
        tallies.extend(standard_tallies(problem, materials))
        openmc_tallies.quiet(settings)

    directory.mkdir(parents=True, exist_ok=True)
    openmc.Model(geometry, materials, settings, tallies).export_to_model_xml(
        directory / "model.xml"
    )


def solve(source: Path, directory: Path, work: Path) -> None:
    """`directory` の model.xml を、入力の `cross-sections` のライブラリとスレッド数で解く。"""
    run = read_run(source)
    work.mkdir(parents=True, exist_ok=True)
    # model.xml にライブラリのパスを書かないので、環境変数で渡す（OpenMC は materials の cross_sections がなければ読む）
    os.environ["OPENMC_CROSS_SECTIONS"] = str(run["cross-sections"])
    openmc.run(
        threads=run["threads"],
        output=False,
        cwd=work,
        path_input=str((directory / "model.xml").resolve()),
    )


def extract(source: Path, work: Path, results: Path) -> None:
    run = read_run(source)
    problem = run["problem"]
    results.mkdir(parents=True, exist_ok=True)
    with openmc.StatePoint(work / f"statepoint.{run['batches']}.h5") as sp:
        k = sp.keff
        if problem.startswith("2"):
            write_pins(sp, results / "pins.msgpack")
        if openmc_tallies.enabled(run):
            # 物質は summary.h5（statepoint が読む）から、run の物質と同じ ID で取る
            materials = list(sp.summary.materials)
            if problem.startswith("1"):
                (fuel,) = fuels_of(materials)
                path, form = results / "fuel-xs.msgpack", "xlc-fuel-xs"
                xs = fuel_xs(sp, fuel)
            else:
                path, form = results / "pin-xs.msgpack", "xlc-pin-xs"
                xs = pin_xs(sp, fuels_of(materials))
            path.write_bytes(msgpack.packb({"format": form, "nuclides": xs}))
    write_k(results / "k.json", problem, k)


if __name__ == "__main__":
    model_run_extract(model, solve, extract, __doc__)
