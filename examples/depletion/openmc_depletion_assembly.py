"""VERA Problem 2 の 17×17 集合体（毒物なし、Pyrex、IFBA、WABA、Gd と、履歴を変えた 2 本）を OpenMC の燃焼（openmc.deplete）で
燃焼させ、xlc の格子計算の燃焼（ADR-0173）と比べる参照解を作る（ADR-0178）。

使い方（references の devShell で）:

    python3 openmc_depletion_assembly.py run <入力> <作業ディレクトリ>
    python3 openmc_depletion_assembly.py extract <入力> <作業ディレクトリ> <結果（msgpack）>
    python3 openmc_depletion_assembly.py transport <入力> <作業ディレクトリ>

<入力> は p2/openmc-2A.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 問題（`problem`。下の CASES）、
openmc_library.py で作った `cross-sections`、燃焼チェーン（`chain`）、初めの重金属の出力密度（`power`、W/g）、燃焼度の点（`burnup`、
MWd/kg、0 より大きく増える列）、1 世代の粒子数（`particles`）、世代数（`batches`）、非アクティブバッチ数（`inactive`）、乱数のシード（`seed`）、
スレッド数（`threads`）。Gd の問題には Gd 燃料のリングの数（`gd-rings`。ほかの問題には書けない）。時間の倍率を測る入力にだけ、
燃料の初めの組成を別の燃焼の結果からコピーする `start-from`（openmc_depletion.py の extract の msgpack）と `start-point`（その燃焼度の
点の番号）。パスは入力のファイルからの相対。

方式は openmc_depletion.py と同じ（CE/CM、CRAM48、一定の出力密度、チェーンの核分裂の Q、連続エネルギーの反応率、収率は 112.5 eV で
熱と速。ADR-0173）。違いは燃焼する物質: 集合体の 289 のセルを、中央からの距離で決まる対称操作の組（`orbit`。xlc が自己遮蔽と MOC で
1 つにまとめるピンの組。ADR-0112）に分け、組と物質の対ごとに 1 つの燃焼する物質にする（ADR-0173 の案 R1）。燃料のほか、
Pyrex、IFBA の被覆、WABA の吸収体も燃焼する物質にする。Gd 燃料はリング（等しい面積）ごとの物質にする。制御棒（AIC）は燃焼しない。
物質の体積は、組のセルの数 × 1 つのセルの中の面積 (cm²)（2 次元なので高さ 1 cm あたり）。

run は作業ディレクトリに depletion_results.h5、燃焼点ごとのピン出力と時間（pins.json）を書く。ピン出力は、燃焼点の始まりの輸送
（`openmc_simulation_n<点>.h5`。OpenMC の `write_bos_data`）の核分裂率のメッシュの tally。extract は燃焼点ごとの k と統計誤差、領域ごとの
輸送する核種の数密度 (1/(barn cm))、ピン出力、ステップごとの輸送と CRAM の時間を msgpack に書く（compare_depletion_assembly.py が
読む）。transport は燃焼せずに初めの組成の 1 回の輸送だけ解き、時間とピン出力の統計誤差を書く（時間の倍率の基準）。
"""

import itertools
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import NamedTuple

import h5py
import msgpack
import numpy as np
import openmc
import openmc.deplete

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "vera"))
sys.path.insert(0, str(HERE))
from compare_depletion_assembly import compared_openmc_names, orbit_std_rms
from openmc_common import eigenvalue_settings, material
from openmc_depletion import TRANSPORT_SECONDS, YIELD_CUTOFF, TimedOperator
from openmc_input import input_path, read_input
from openmc_model import (
    AIC,
    CARBON,
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
    MODERATOR,
    MODERATOR_DENSE,
    PELLET,
    PYREX,
    PYREX_OCTANTS,
    PYREX_RADII,
    ROD_CLAD,
    SS304,
    ZIRCALOY,
    N,
    concentric,
    octant_images,
    pin_mesh,
    pin_rows,
    temperature_mode,
)
from openmc_run import fuel_mean

KEYS = {
    "problem",
    "cross-sections",
    "chain",
    "power",
    "burnup",
    "particles",
    "batches",
    "inactive",
    "seed",
    "threads",
}
OPTIONAL = {"gd-rings", "start-from", "start-point"}

# 個数密度 (/(barn cm))。仕様書の表 P2-4（Gd 燃料は 5% Gd2O3 と 95% UO2 の 1.8% 濃縮の燃料）
GD_FUEL = {
    "U234": 3.18096e-06,
    "U235": 3.90500e-04,
    "U236": 1.79300e-06,
    "U238": 2.10299e-02,
    "Gd152": 3.35960e-06,
    "Gd154": 3.66190e-05,
    "Gd155": 2.48606e-04,
    "Gd156": 3.43849e-04,
    "Gd157": 2.62884e-04,
    "Gd158": 4.17255e-04,
    "Gd160": 3.67198e-04,
    "O16": 4.53705e-02,
}
# IFBA の ZrB2（表 P2-4 の IFBA）。厚さ 10 μm（表 11）
IFBA = {
    "B10": 2.16410e-02,
    "B11": 1.96824e-02,
    "Zr90": 1.06304e-02,
    "Zr91": 2.31824e-03,
    "Zr92": 3.54348e-03,
    "Zr94": 3.59100e-03,
    "Zr96": 5.78528e-04,
}
IFBA_OUTER = 0.4106
# WABA の B4C-Al2O3（表 P2-4）。天然の C（6000 の 3.77001E-03）は C-12 と C-13 に分ける（ADR-0079）
WABA = {
    "B10": 2.98553e-03,
    "B11": 1.21192e-02,
    "C12": 3.77001e-03 * CARBON["C12"],
    "C13": 3.77001e-03 * CARBON["C13"],
    "O16": 5.85563e-02,
    "Al27": 3.90223e-02,
}
# WABA の棒（表 13）: 内管の内と外、吸収体の内と外、被覆管の内と外の半径。内管の中は減速材、内管と吸収体の間、吸収体と被覆管の間は
# He（ADR-0120 に従い真空）
WABA_RADII = (0.286, 0.339, 0.353, 0.404, 0.418, 0.484)

ALL_GUIDE_TUBES = frozenset(GUIDE_TUBES)
# 図 P2-2 の 1/8 の中の位置（中央から行と列の数。列 ≤ 行）。2M は 128 本、2N は 104 本の IFBA と 20 本の WABA（案内管の 24 のうち計装管の
# 隣の (3, 3) を除く）、2O と 2P は Gd の棒
IFBA_128 = {
    (1, 0),
    (2, 0),
    (3, 1),
    (3, 2),
    (4, 0),
    (4, 3),
    (5, 0),
    (5, 3),
    (5, 4),
    (6, 1),
    (6, 2),
    (6, 4),
    (6, 5),
    (7, 0),
    (7, 3),
    (7, 6),
    (8, 1),
    (8, 4),
    (8, 8),
}
IFBA_104 = IFBA_128 - {(7, 6), (8, 1), (8, 4)}
WABA_20 = {(3, 0), (5, 5), (6, 0), (6, 3)}
GD_12 = {(4, 2), (6, 6)}
GD_24 = {(2, 2), (4, 4), (5, 0), (6, 6), (7, 3)}


class Case(NamedTuple):
    """問題の状態と配置。温度は燃料、被覆管、減速材、構造材で同じ（仕様書の表 P2-1）。"""

    temperature: float
    moderator: dict[str, float]
    # 案内管に入れるもの（`pyrex`、`waba`、`aic`）の種類と、その案内管の位置（行, 列。0 始まり）
    insert: str | None
    insert_positions: frozenset[tuple[int, int]]
    ifba: frozenset[tuple[int, int]]
    gd: frozenset[tuple[int, int]]


def case(
    temperature: float,
    moderator: dict[str, float],
    insert: str | None = None,
    insert_octant: set[tuple[int, int]] | None = None,
    ifba: set[tuple[int, int]] | None = None,
    gd: set[tuple[int, int]] | None = None,
) -> Case:
    positions = frozenset(octant_images(insert_octant or set()))
    return Case(
        temperature,
        moderator,
        insert,
        positions,
        frozenset(octant_images(ifba or set())),
        frozenset(octant_images(gd or set())),
    )


# 問題の表。表にない問題は解かない。2A と履歴の 2 本は 565 K（温度の内挿のモード。ADR-0080）、ほかは 600 K（表 P2-1）。履歴の 2 本は
# 2A の集合体の 1 つの条件だけを変える: 2A-density は減速材の密度を 0.743 から 0.661 g/cc（仕様書の 2B〜2D の値。温度は 565 K のまま）、
# 2A-rods は 24 本の案内管すべてに AIC の制御棒を入れる（2G と同じ配置）
CASES = {
    "2A": case(565.0, MODERATOR_DENSE),
    "2A-density": case(565.0, MODERATOR),
    "2A-rods": case(565.0, MODERATOR_DENSE, "aic", PYREX_OCTANTS[24]),
    "2F": case(600.0, MODERATOR_DENSE, "pyrex", PYREX_OCTANTS[24]),
    "2M": case(600.0, MODERATOR_DENSE, ifba=IFBA_128),
    "2N": case(600.0, MODERATOR_DENSE, "waba", WABA_20, ifba=IFBA_104),
    "2O": case(600.0, MODERATOR_DENSE, gd=GD_12),
    "2P": case(600.0, MODERATOR_DENSE, gd=GD_24),
}
for _name, _case in CASES.items():
    if not _case.insert_positions <= ALL_GUIDE_TUBES:
        raise SystemExit(f"{_name}: an insert outside the guide tubes")
    if (_case.ifba | _case.gd) & (ALL_GUIDE_TUBES | {INSTRUMENT_TUBE}):
        raise SystemExit(f"{_name}: fuel poison in a tube position")
    if _case.ifba & _case.gd:
        raise SystemExit(f"{_name}: IFBA and gadolinia in the same pin")

# 燃焼度の点 (MWd/kg)。0 より大きく増え、最後は 60（docs/plans/m4-depletion.md の「OpenMC の燃焼の参照解の物理コア時間」の 60 GWd/t）。
# Xe-135 と Sm-149 の立ち上がりの細かいステップから始める（1 MWd/kg は出力密度 38.41 W/g で約 26 日）。70 点は、IFBA の B-10 と
# Gd が燃える 16 MWd/kg までを 40 点の約 2 倍の細かさにする
BURNUP_40 = [
    0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
    14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
    32.0, 34.0, 36.0, 38.0, 40.0, 44.0, 48.0, 52.0, 56.0, 60.0,
]  # fmt: skip
BURNUP_70 = [
    *(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0),
    *(1.0 + 0.25 * i for i in range(1, 13)),
    *(4.0 + 0.5 * i for i in range(1, 25)),
    *(16.0 + 1.0 * i for i in range(1, 9)),
    *(24.0 + 2.0 * i for i in range(1, 19)),
]


def orbit(i: int, j: int) -> tuple[int, int]:
    """セル（行 i、列 j。0 始まり）が属する対称操作の組。中央からの行と列の距離を大きい順に並べた対 (a, b)（b ≤ a）。"""
    centre = N // 2
    a, b = sorted((abs(i - centre), abs(j - centre)), reverse=True)
    return a, b


def orbit_size(o: tuple[int, int]) -> int:
    """組のセルの数: 中央 1、軸上と対角線上 4、ほかは 8。"""
    a, b = o
    if a == 0:
        return 1
    return 4 if b in (0, a) else 8


def cell_kind(problem: str, i: int, j: int) -> str:
    """セル（行 i、列 j）の種類: `fuel`、`ifba`、`gd`（燃料のピン）、`gt`（案内管）、`it`（計装管）、`pyrex`、`waba`、`aic`
    （案内管に入れたもの）。"""
    if problem not in CASES:
        raise SystemExit(f"unknown problem {problem}; expected one of {list(CASES)}")
    c = CASES[problem]
    if (i, j) == INSTRUMENT_TUBE:
        return "it"
    if (i, j) in ALL_GUIDE_TUBES:
        return c.insert if (i, j) in c.insert_positions and c.insert else "gt"
    if (i, j) in c.ifba:
        return "ifba"
    if (i, j) in c.gd:
        return "gd"
    return "fuel"


def gd_radii(rings: int) -> list[float]:
    """ペレットを等しい面積に分けた Gd 燃料のリングの外半径 (cm)。"""
    return [PELLET * math.sqrt((k + 1) / rings) for k in range(rings)]


class Zone(NamedTuple):
    """セルの中の同心円の帯。外半径 `radius`（None は最も外）、MOC の分割数 `rings`（xlc の入力が使う）。"""

    material: str
    radius: float | None
    rings: int


# 燃料のピンの外側（ギャップ、被覆管、減速材）。lattice の入力の VERA の 2A と同じ分割
def _outside() -> list[Zone]:
    return [
        Zone("void", CLAD_INNER, 1),
        Zone("zircaloy", CLAD_OUTER, 1),
        Zone("water", 0.63, 3),
        Zone("water", 0.75, 2),
        Zone("water", None, 1),
    ]


def _tube_outside() -> list[Zone]:
    return [Zone("water", 0.75, 2), Zone("water", None, 1)]


def cell_zones(kind: str, gd_rings: int | None = None) -> list[Zone]:
    """セルの種類の同心円の帯（内側から）。OpenMC のモデルと xlc の入力の形状はこれから作る。"""
    match kind:
        case "fuel":
            return [Zone("fuel", PELLET, 5), *_outside()]
        case "ifba":
            return [Zone("fuel", PELLET, 5), Zone("ifba", IFBA_OUTER, 1), *_outside()]
        case "gd":
            if gd_rings is None:
                raise SystemExit("a gadolinia pin needs gd-rings")
            fuel = [Zone(f"gd{k}", r, 1) for k, r in enumerate(gd_radii(gd_rings))]
            return [*fuel, *_outside()]
        case "gt":
            return [
                Zone("water", GUIDE_INNER, 3),
                Zone("zircaloy", GUIDE_OUTER, 1),
                *_tube_outside(),
            ]
        case "it":
            return [
                Zone("water", INSTRUMENT_INNER, 3),
                Zone("zircaloy", INSTRUMENT_OUTER, 1),
                *_tube_outside(),
            ]
        case "pyrex":
            fills = ["void", "ss304", "void", "pyrex", "void", "ss304"]
            rings = [1, 1, 1, 2, 1, 1]
            inner = [Zone(f, r, n) for f, r, n in zip(fills, PYREX_RADII, rings)]
            return [
                *inner,
                Zone("water", GUIDE_INNER, 1),
                Zone("zircaloy", GUIDE_OUTER, 1),
                *_tube_outside(),
            ]
        case "waba":
            fills = ["water", "zircaloy", "void", "waba", "void", "zircaloy"]
            rings = [2, 1, 1, 2, 1, 1]
            inner = [Zone(f, r, n) for f, r, n in zip(fills, WABA_RADII, rings)]
            return [
                *inner,
                Zone("water", GUIDE_INNER, 1),
                Zone("zircaloy", GUIDE_OUTER, 1),
                *_tube_outside(),
            ]
        case "aic":
            clad_in, clad_out = ROD_CLAD
            return [
                Zone("aic", 0.382, 5),
                Zone("void", clad_in, 1),
                Zone("ss304", clad_out, 1),
                Zone("water", GUIDE_INNER, 1),
                Zone("zircaloy", GUIDE_OUTER, 1),
                *_tube_outside(),
            ]
    raise SystemExit(f"unknown cell kind {kind}")


def is_burnable(material_name: str) -> bool:
    """燃焼する物質の名前（燃料、Gd 燃料のリング、IFBA の被覆、Pyrex、WABA の吸収体）。"""
    return material_name in ("fuel", "ifba", "pyrex", "waba") or (
        material_name.startswith("gd") and material_name[2:].isdigit()
    )


class Region(NamedTuple):
    """燃焼する物質: 種類（xlc の入力の物質の名前）、対称操作の組、組のセルの数、体積の和 (cm²)。"""

    kind: str
    orbit: tuple[int, int]
    count: int
    volume: float

    @property
    def name(self) -> str:
        return f"{self.kind}@{self.orbit[0]}-{self.orbit[1]}"


def burnable_areas(zones: list[Zone]) -> dict[str, float]:
    """セルの帯の、燃焼する物質ごとの面積 (cm²)（現れる順）。"""
    areas: dict[str, float] = {}
    inner = 0.0
    for z in zones:
        if z.radius is None:
            break
        if is_burnable(z.material):
            area = math.pi * (z.radius**2 - inner**2)
            areas[z.material] = areas.get(z.material, 0.0) + area
        inner = z.radius
    return areas


def regions(run: dict) -> list[Region]:
    """燃焼する物質の領域。集合体を上の行から左の列へ走り、(組, 物質) の最初に現れる順（xlc の燃焼の領域の順）。"""
    problem = run["problem"]
    found: dict[tuple[str, tuple[int, int]], Region] = {}
    for i, j in itertools.product(range(N), range(N)):
        o = orbit(i, j)
        areas = burnable_areas(
            cell_zones(cell_kind(problem, i, j), run.get("gd-rings"))
        )
        for name, area in areas.items():
            if (name, o) not in found:
                found[(name, o)] = Region(name, o, orbit_size(o), orbit_size(o) * area)
    return list(found.values())


def composition(name: str, c: Case) -> dict[str, float]:
    """物質の名前の組成（個数密度）。"""
    if name == "fuel":
        return FUEL
    if name.startswith("gd") and name[2:].isdigit():
        return GD_FUEL
    return {
        "ifba": IFBA,
        "pyrex": PYREX,
        "waba": WABA,
        "aic": AIC,
        "ss304": SS304,
        "zircaloy": ZIRCALOY,
        "water": c.moderator,
    }[name]


def library_nuclides_at(cross_sections: Path, temperature: float) -> set[str]:
    """ライブラリが温度 `temperature` (K) のデータを持つ核種。OpenMC は温度を kT で持ち、K に直すと丸めでずれるので 1 K の幅を許す
    （openmc_depletion.py の TimedOperator と同じ規則）。"""
    found = set()
    for library in openmc.data.DataLibrary.from_xml(cross_sections).libraries:
        if library["type"] != "neutron":
            continue
        with h5py.File(library["path"], "r") as f:
            for name in library["materials"]:
                temperatures = [float(k.rstrip("K")) for k in f[name]["kTs"]]
                if any(abs(t - temperature) <= 1.0 for t in temperatures):
                    found.add(name)
    return found


def burned_fuel(run: dict, c: Case) -> dict[str, float]:
    """`start-from` の燃焼の結果（ピンセルの燃料）の `start-point` の点の組成のうち、ライブラリが燃料の温度で持つ核種の
    個数密度（複数のリングは体積の重みで平均する）。"""
    result = msgpack.unpackb(run["start-from"].read_bytes())
    rings = result["rings"]
    point = run["start-point"]
    if not 0 <= point < len(result["burnup"]):
        raise SystemExit(
            f"start-point {point} is outside the {len(result['burnup'])} burnup points of {run['start-from']}"
        )
    total = sum(r["volume"] for r in rings)
    mean: dict[str, float] = {}
    for ring in rings:
        for nuclide, values in ring["densities"].items():
            mean[nuclide] = (
                mean.get(nuclide, 0.0) + values[point] * ring["volume"] / total
            )
    kept = library_nuclides_at(run["cross-sections"], c.temperature)
    burned = {n: d for n, d in mean.items() if d > 0 and n in kept}
    if not burned:
        raise SystemExit(
            f"{run['start-from']} point {point} has no nuclide at {c.temperature} K"
        )
    return burned


def build(run: dict) -> tuple[openmc.Model, list[Region]]:
    """問題のモデル（燃焼する物質は depletable で体積を持つ）と、燃焼する物質の領域（`regions` と同じ順）。"""
    problem = run["problem"]
    c = CASES[problem]
    temperature = c.temperature
    water = material("moderator", composition("water", c), temperature)
    water.add_s_alpha_beta("c_H_in_H2O")
    fixed = {"water": water}

    def fixed_material(name: str) -> openmc.Material | None:
        if name == "void":
            return None
        if name not in fixed:
            label = {"zircaloy": "zircaloy-4"}.get(name, name)
            fixed[name] = material(label, composition(name, c), temperature)
        return fixed[name]

    start = burned_fuel(run, c) if "start-from" in run else None
    region_list = regions(run)
    burn: dict[tuple[str, tuple[int, int]], openmc.Material] = {}
    for r in region_list:
        densities = (
            start if start is not None and r.kind == "fuel" else composition(r.kind, c)
        )
        m = material(r.name, densities, temperature)
        m.depletable = True
        m.volume = r.volume
        burn[(r.kind, r.orbit)] = m
    if start is not None and {r.kind for r in region_list} != {"fuel"}:
        raise SystemExit("start-from is for the problems that burn only fuel")

    universes: dict[tuple[str, tuple[int, int]], openmc.Universe] = {}
    rows = []
    for i in range(N):
        row = []
        for j in range(N):
            kind, o = cell_kind(problem, i, j), orbit(i, j)
            if (kind, o) not in universes:
                zones = cell_zones(kind, run.get("gd-rings"))
                while zones[-1].material == "water":
                    zones.pop()
                layers = [
                    (
                        burn[(z.material, o)]
                        if is_burnable(z.material)
                        else fixed_material(z.material),
                        z.radius,
                    )
                    for z in zones
                ]
                universes[(kind, o)] = concentric(layers, water)
            row.append(universes[(kind, o)])
        rows.append(row)

    lattice = openmc.RectLattice()
    lattice.pitch = (1.26, 1.26)
    lattice.lower_left = (-N * 1.26 / 2, -N * 1.26 / 2)
    lattice.universes = rows
    lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=water)])
    half = N * 1.26 / 2 + HALF_GAP
    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    geometry = openmc.Geometry([openmc.Cell(fill=lattice, region=-box)])

    settings = eigenvalue_settings(
        run["particles"], run["batches"], run["inactive"], run["seed"], half
    )
    settings.temperature = temperature_mode(temperature)
    settings.output = {"summary": False, "tallies": False}
    tally = openmc.Tally(name="fission")
    tally.filters = [openmc.MeshFilter(pin_mesh())]
    tally.scores = ["fission"]
    materials = openmc.Materials(
        [*burn.values(), *(m for m in fixed.values() if m is not None)]
    )
    model = openmc.Model(geometry, materials, settings, openmc.Tallies([tally]))
    return model, region_list


def read(source: Path) -> dict:
    run = read_input(source, KEYS, OPTIONAL)
    problem = run["problem"]
    if problem not in CASES:
        raise SystemExit(f"unknown problem {problem}; expected one of {list(CASES)}")
    has_gd = bool(CASES[problem].gd)
    if has_gd != ("gd-rings" in run):
        raise SystemExit(
            f"{problem}: gd-rings is required exactly for gadolinia problems"
        )
    if "gd-rings" in run and not (
        isinstance(run["gd-rings"], int) and run["gd-rings"] > 0
    ):
        raise SystemExit(f"gd-rings must be a positive integer, got {run['gd-rings']}")
    if ("start-from" in run) != ("start-point" in run):
        raise SystemExit("start-from and start-point go together")
    if "start-point" in run and not (
        isinstance(run["start-point"], int) and run["start-point"] >= 0
    ):
        raise SystemExit(
            f"start-point must be a non-negative integer, got {run['start-point']}"
        )
    burnup = run["burnup"]
    if (
        not burnup
        or burnup[0] <= 0
        or any(b <= a for a, b in itertools.pairwise(burnup))
    ):
        raise SystemExit(f"burnup must be positive and increasing, got {burnup}")
    if not run["power"] > 0:
        raise SystemExit(f"power must be positive, got {run['power']}")
    run["cross-sections"] = input_path(source, run["cross-sections"])
    run["chain"] = input_path(source, run["chain"])
    if "start-from" in run:
        run["start-from"] = input_path(source, run["start-from"])
    return run


def fuel_temperature(run: dict) -> float:
    """問題の温度 (K)。燃料のデータのある核種だけを輸送する規則（TimedOperator）が使う。"""
    return CASES[run["problem"]].temperature


def operator(run: dict, model: openmc.Model, prev=None) -> TimedOperator:
    TimedOperator.fuel_temperature = fuel_temperature(run)
    return TimedOperator(
        model,
        chain_file=run["chain"],
        prev_results=prev,
        normalization_mode="fission-q",
        fission_yield_mode="cutoff",
        fission_yield_opts={"cutoff": YIELD_CUTOFF},
        reaction_rate_mode="direct",
    )


def bos_points(directory: Path) -> list[dict]:
    """`directory` の `openmc_simulation_n<番号>.h5`（燃焼点の始まりの輸送の statepoint）の、番号順の、ピン出力
    （燃料のピンの平均で 1。[行（上から）][列]）の平均と統計誤差、k、時間（秒）。"""
    files = sorted(
        directory.glob("openmc_simulation_n*.h5"),
        key=lambda p: int(p.stem.rsplit("_n", 1)[1]),
    )
    points = []
    for path in files:
        with openmc.StatePoint(path) as sp:
            mean, std = pin_rows(sp, "mean"), pin_rows(sp, "std_dev")
            scale = fuel_mean(mean)
            runtime = {k: float(v) for k, v in sp.runtime.items()}
            points.append(
                {
                    "mean": (mean / scale).tolist(),
                    "std": (std / scale).tolist(),
                    "k": sp.keff.nominal_value,
                    "k-std": sp.keff.std_dev,
                    "runtime": runtime,
                }
            )
    return points


def write_pins(directory: Path, points: list[dict]) -> None:
    (directory / "pins.json").write_text(json.dumps({"points": points}) + "\n")


def read_pins(directory: Path) -> list[dict]:
    return json.loads((directory / "pins.json").read_text())["points"]


# 区間に分けた実行（openmc_depletion_resume.py）が extract に渡すファイル
RESULT_FILES = ("depletion_results.h5", "transport-seconds.json", "pins.json")


def collect_interval(step: Path, k: int, previous: Path | None) -> None:
    """区間 k（`step`）の終わりに、燃焼点のピン出力を、前の区間 `previous` までの分に続けて `step/pins.json` に書く。区間 k は点 k − 1 の
    始まりの輸送（と、最後の区間は点 k の輸送）を持つ。"""
    before = read_pins(previous) if previous is not None else []
    write_pins(step, [*before[: k - 1], *bos_points(step)])


def clean(directory: Path) -> None:
    """輸送ごとに上書きされる statepoint（粒子の源を含み、数十 MB）を消す。"""
    for statepoint in directory.glob("statepoint.*.h5"):
        statepoint.unlink()


def solve(source: Path, work: Path) -> None:
    run = read(source)
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    os.environ["OMP_NUM_THREADS"] = str(run["threads"])
    openmc.config["cross_sections"] = run["cross-sections"]
    model, _ = build(run)
    op = operator(run, model)
    burnup = run["burnup"]
    steps = [b - a for a, b in itertools.pairwise([0.0, *burnup])]
    integrator = openmc.deplete.CECMIntegrator(
        op,
        steps,
        power_density=run["power"],
        timestep_units="MWd/kg",
        solver="cram48",
    )
    start = time.perf_counter()
    integrator.integrate(final_step=True)
    total = time.perf_counter() - start
    (work / "transport-seconds.json").write_text(
        json.dumps(
            {
                "transport": TRANSPORT_SECONDS,
                "total": total,
                "transportable": TimedOperator.transportable,
            }
        )
        + "\n"
    )
    write_pins(work, bos_points(work))
    clean(work)
    print(
        f"transport {sum(TRANSPORT_SECONDS):.1f} s in {len(TRANSPORT_SECONDS)} solves, "
        f"total {total:.1f} s"
    )


def extract(source: Path, work: Path, output: Path) -> None:
    run = read(source)
    results = openmc.deplete.Results(work / "depletion_results.h5")
    _, region_list = build_for_ids(run)
    times, k = results.get_keff(time_units="s")
    seconds = json.loads((work / "transport-seconds.json").read_text())
    keep = set(seconds["transportable"]) | COMPARED
    first = results[0]
    nuclides = [n for n in first.index_nuc if n in keep]
    if missing := COMPARED - set(nuclides):
        raise SystemExit(
            f"the chain results lack the compared nuclides {sorted(missing)}"
        )
    regions_out = []
    for r, mat_id in region_list:
        densities = {name: [] for name in nuclides}
        for step in results:
            atoms = step.data[0][step.index_mat[str(mat_id)]]
            for name in nuclides:
                # 原子数 / 体積 (cm³) × 1e-24 = 1/(barn cm)
                densities[name].append(
                    float(atoms[step.index_nuc[name]]) / r.volume * 1e-24
                )
        regions_out.append(
            {
                "material": r.kind,
                "orbit": list(r.orbit),
                "count": r.count,
                "volume": r.volume,
                "densities": densities,
            }
        )
    pins = read_pins(work)
    output.write_bytes(
        msgpack.packb(
            {
                "format": "openmc-depletion-assembly",
                "problem": run["problem"],
                "burnup": [0.0, *run["burnup"]],
                "time": np.asarray(times).tolist(),
                "k": k[:, 0].tolist(),
                "k-std": k[:, 1].tolist(),
                "regions": regions_out,
                "pins": [p["mean"] for p in pins],
                "pins-std": [p["std"] for p in pins],
                "runtime": [p["runtime"] for p in pins],
                "transport-seconds": seconds["transport"],
                "total-seconds": seconds["total"],
                "cram-seconds": [float(r.proc_time) for r in results],
                "transportable": seconds["transportable"],
            }
        )
    )
    for b, kk in zip([0.0, *run["burnup"]], k):
        print(f"burnup {b} k {kk[0]:.6f} +- {kk[1]:.6f}")


# 輸送する核種でなくても、比べる核種（compare_depletion_assembly.py が名前を引く）
COMPARED = compared_openmc_names()


def build_for_ids(run: dict) -> tuple[openmc.Model, list[tuple[Region, int]]]:
    """`build` の領域と、その物質の ID（depletion_results.h5 の物質の名前）。ID は作る順で決まる（同じ入力から同じ）。"""
    model, region_list = build(run)
    ids = {m.name: m.id for m in model.materials}
    return model, [(r, ids[r.name]) for r in region_list]


def transport(source: Path, work: Path) -> None:
    """燃焼せずに、初めの組成で 1 回だけ輸送を解く（時間の倍率の基準。燃焼の tally を入れない）。"""
    run = read(source)
    work.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = str(run["threads"])
    openmc.config["cross_sections"] = run["cross-sections"]
    model, _ = build(run)
    model.run(cwd=work, threads=run["threads"], output=False)
    with openmc.StatePoint(work / f"statepoint.{run['batches']}.h5") as sp:
        mean, std = pin_rows(sp, "mean"), pin_rows(sp, "std_dev")
        relative = std[mean > 0] / mean[mean > 0]
        runtime = {k: float(v) for k, v in sp.runtime.items()}
        summary = {
            "problem": run["problem"],
            "particles": run["particles"],
            "batches": run["batches"],
            "inactive": run["inactive"],
            "threads": run["threads"],
            "k": sp.keff.nominal_value,
            "k-std": sp.keff.std_dev,
            "pin-rms-std": float(np.sqrt(np.mean(relative**2))),
            "orbit-rms-std": orbit_std_rms(mean.tolist(), std.tolist()),
            "runtime": runtime,
        }
    (work / "transport.json").write_text(json.dumps(summary, indent=1) + "\n")
    clean(work)
    seconds = runtime["simulation"]
    print(
        f"{run['problem']}: simulation {seconds:.1f} s for {run['batches']} batches "
        f"({seconds / run['batches']:.3f} s per batch), "
        f"k {summary['k']:.6f} +- {summary['k-std']:.6f}, "
        f"pin rms std {summary['pin-rms-std'] * 100:.3f}%, "
        f"orbit-mean rms std {summary['orbit-rms-std'] * 100:.3f}%"
    )


if __name__ == "__main__":
    match sys.argv[1:]:
        case ["run", source, work]:
            solve(Path(source).resolve(), Path(work).resolve())
        case ["extract", source, work, output]:
            extract(
                Path(source).resolve(), Path(work).resolve(), Path(output).resolve()
            )
        case ["transport", source, work]:
            transport(Path(source).resolve(), Path(work).resolve())
        case _:
            raise SystemExit(__doc__)
