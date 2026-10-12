"""VERA Problem 5A-2D、5B-2D、5C-2D（炉心の高さの中央の 2D の 1/4 炉心、零出力）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_core_2d.py run [--restart] [--threads N] <入力> <作業ディレクトリ>
        python3 openmc_core_2d.py extract <入力> <作業ディレクトリ> <結果のディレクトリ>

<入力> は p5/openmc-5A-2D.toml の形の TOML で、openmc_model.py の入力と同じキー（`problem`、`cross-sections`、
`particles`、`batches`、`inactive`、`seed`、`threads`）をすべて書く（既定値はない。ADR-0130）。途中の statepoint を書く間隔は、任意のキー `statepoint-interval`（世代数）。
`--restart` と `--threads` は openmc_resume.py。run と extract の分け方も
openmc_model.py と同じ（ADR-0123）。

集合体の径方向の形（ピン、案内管、計装管、Pyrex の棒と制御棒の寸法）と、問題 2 と同じ組成は openmc_model.py から
import する。openmc_model.py の形状は集合体 1 体で外周が鏡面反射なので、炉心の並べ方と炉心の外の構造はここに書く。

形状（仕様書の問題 5-2D、1.12、1.13 節、表 16、図 M3-1、M3-2）:
- 図 M3-1 の 1/4 炉心（H〜A 列、8〜15 行）を、炉心の中心を原点に x ≥ 0、y ≤ 0 に置く（図と同じ向きで、A 列が +x、
  15 行が −y）。x = 0 と y = 0 は鏡面反射（H 列と 8 行の集合体はその面で半分に切れる）。z は一様（2D）
- 集合体は燃料の濃縮度（2.11%、2.619%、3.10%）と Pyrex の棒の本数（図 5 の配置）で決まる。D バンクの位置（H8、D8、
  H12、D12）の 24 本の案内管には、5B-2D で AIC、5C-2D で B4C の制御棒を入れる（5A-2D は空）。計装管のシンブルと
  シンブルプラグは入れない（仕様書のとおり）
- バッフル（SS304、厚さ 2.85 cm）は、集合体の外周（ピッチ 21.50 cm の正方形の和）から水の隙間 0.19 cm をおいて、
  炉心の階段の形に沿う。集合体の外周を正方形で d だけ太らせた形（ミンコフスキー和）の d = 0.19 と 0.19 + 2.85 の
  間とする（平らな板が凸の角でも凹の角でも直角に会う形）
- その外は水、バレル（SS304、半径 187.96〜193.68）、水、中性子パッド（SS304、半径 194.64〜201.63、45° を中心に
  32°）、水、圧力容器の内張り（SS304、219.15〜219.71）、圧力容器（炭素鋼、219.71〜241.70）。その外は真空
  （付録 J の Case 1。仕様書の参照解）

Hf とギャップの He は除いて真空とし（ADR-0054）、Pyrex の棒と制御棒の中の He も同じとする（ADR-0120）。すべて 565 K、
減速材は 0.743 g/cc、ホウ素 1300 ppm（表 M3-1）。

extract は結果のディレクトリに次を書く。ピンと集合体の出力は、鏡面反射で 1/4 から全炉心（15×15 の集合体）に
広げ、行は上（1 行）から、列は左（R 列）から並べる（core の `pin-powers` と同じく集合体の間の余白を除く）:
- k.json: k と統計誤差
- pins.msgpack: ピンの核分裂率（`xlc-pin-powers`。255 × 255、燃料のないセルは 0、燃料のピンの平均で 1）。対称面で
  半分（中心の集合体の中央の行と列）に切れたピンは、モデルの中の部分の面積で割って全体のピンに直す
- pins-std.msgpack: その統計誤差（同じ形式と規格化。規格化の定数の誤差は含めない）
- assemblies.msgpack: 集合体の出力（`xlc-pin-powers`。15 × 15、炉心の外は 0、集合体の平均で 1）
- assemblies-std.msgpack: その統計誤差
"""

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import openmc
import openmc_tallies

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import material
from openmc_model import (
    AIC,
    B4C,
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
    PYREX_OCTANTS,
    SS304,
    ZIRCALOY,
    N,
    concentric,
    fuel_layers,
    octant_images,
    pin_lattice,
    pyrex_layers,
    rod_layers,
    temperature_mode,
    tube,
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

TEMPERATURE = 565.0
ASSEMBLY_PITCH = N * PITCH + 2 * HALF_GAP

# 表 M3-4 の個数密度 (/(barn cm))。3.10% は openmc_model.py の FUEL と同じ
FUEL_211 = {
    "U234": 4.04814e-06,
    "U235": 4.88801e-04,
    "U236": 2.23756e-06,
    "U238": 2.23844e-02,
    "O16": 4.57591e-02,
}
FUEL_2619 = {
    "U234": 5.09503e-06,
    "U235": 6.06709e-04,
    "U236": 2.76809e-06,
    "U238": 2.22663e-02,
    "O16": 4.57617e-02,
}
# 炭素鋼（CS508）。天然の C は C-12 と C-13 に分ける（ADR-0079）
CARBON_STEEL = with_carbon(
    {
        "Fe54": 4.89841e-03,
        "Fe56": 7.68945e-02,
        "Fe57": 1.77583e-03,
        "Fe58": 2.36330e-04,
    },
    3.93598e-03,
)
FUELS = {"2.1": FUEL_211, "2.6": FUEL_2619, "3.1": FUEL}

# 図 M3-1 の 1/4 炉心: 行（8〜15）ごとに H 列から（濃縮度, Pyrex の棒の本数）
LOADING = (
    (("2.1", 0), ("2.6", 20), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("3.1", 12)),
    (("2.6", 20), ("2.1", 0), ("2.6", 24), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("3.1", 24), ("3.1", 0)),
    (("2.1", 0), ("2.6", 24), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("2.6", 16), ("2.1", 0), ("3.1", 8)),
    (("2.6", 20), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("2.6", 20), ("2.1", 0), ("3.1", 16), ("3.1", 0)),
    (("2.1", 0), ("2.6", 20), ("2.1", 0), ("2.6", 20), ("2.6", 0), ("2.6", 24), ("3.1", 0)),
    (("2.6", 20), ("2.1", 0), ("2.6", 16), ("2.1", 0), ("2.6", 24), ("3.1", 12), ("3.1", 0)),
    (("2.1", 0), ("3.1", 24), ("2.1", 0), ("3.1", 16), ("3.1", 0), ("3.1", 0)),
    (("3.1", 12), ("3.1", 0), ("3.1", 8), ("3.1", 0)),
)  # fmt: skip
# D バンクの位置（1/4 炉心の行, 列。0 始まり、H8 が (0, 0)。図 M3-1 の右）
BANK_D = {(0, 0), (0, 4), (4, 0), (4, 4)}
# 制御棒の吸収体（表 8）
RODS = {"5A-2D": None, "5B-2D": ("aic", AIC, 0.382), "5C-2D": ("b4c", B4C, 0.373)}

# 炉心の外の構造の半径 (cm)（表 16）
BAFFLE_GAP = 0.19
BAFFLE_THICKNESS = 2.85
BARREL = (187.96, 193.68)
PAD = (194.64, 201.63)
PAD_ANGLES = (45.0 - 16.0, 45.0 + 16.0)
LINER = 219.15
VESSEL = (219.71, 241.70)


def assembly(
    fuel: openmc.Material,
    clad: openmc.Material,
    water: openmc.Material,
    pyrex: int,
    rod: openmc.Universe | None,
    pyrex_cell: openmc.Universe,
) -> openmc.Universe:
    """集合体 1 体のユニバース（中心が原点、外周の水のギャップを含む）。`pyrex` は Pyrex の棒の本数（0 はなし）、
    `rod` は 24 本の案内管に入れる制御棒の案内管のセル。"""
    pin = concentric(fuel_layers(fuel, clad), water)
    universes = [[pin] * N for _ in range(N)]
    empty = tube([], GUIDE_INNER, GUIDE_OUTER, clad, water)
    poisoned = octant_images(PYREX_OCTANTS[pyrex]) if pyrex else set()
    if len(poisoned) != pyrex or poisoned - set(GUIDE_TUBES):
        raise SystemExit(f"bad Pyrex pattern for {pyrex} rods")
    if rod is not None and pyrex:
        raise SystemExit("Pyrex rods in a rodded assembly")
    for i, j in GUIDE_TUBES:
        universes[i][j] = rod or (pyrex_cell if (i, j) in poisoned else empty)
    i, j = INSTRUMENT_TUBE
    universes[i][j] = tube([], INSTRUMENT_INNER, INSTRUMENT_OUTER, clad, water)
    return openmc.Universe(cells=[openmc.Cell(fill=pin_lattice(universes, water))])


def row_extents() -> list[tuple[float, float, float]]:
    """1/4 炉心の行ごとの (下端の y, 上端の y, 右端の x)。H 列と 8 行は対称面で切れる。"""
    half = ASSEMBLY_PITCH / 2
    out = []
    for r, row in enumerate(LOADING):
        top = 0.0 if r == 0 else -(r * ASSEMBLY_PITCH - half)
        bottom = -(r * ASSEMBLY_PITCH + half)
        out.append((bottom, top, (len(row) - 1) * ASSEMBLY_PITCH + half))
    return out


def grown(d: float) -> openmc.Region:
    """1/4 炉心の集合体の和を、正方形で d だけ太らせた形（対称面の側は太らせない。対称面の外は形状の外）。"""
    boxes = []
    for bottom, top, right in row_extents():
        boxes.append(
            +openmc.YPlane(bottom - d)
            & -openmc.YPlane(top + d if top < 0 else 1.0)
            & -openmc.XPlane(right + d)
        )
    return openmc.Union(boxes)


def geometry_of(problem: str, m: dict[str, Any]) -> openmc.Geometry:
    """問題の形状。`m` は物質。"""
    water, clad, steel = m["water"], m["clad"], m["steel"]
    pyrex_cell = tube(pyrex_layers(TEMPERATURE), GUIDE_INNER, GUIDE_OUTER, clad, water)
    rod = None
    if RODS[problem] is not None:
        name, nuclides, radius = RODS[problem]
        rod = tube(
            rod_layers(name, nuclides, radius, TEMPERATURE),
            GUIDE_INNER,
            GUIDE_OUTER,
            clad,
            water,
        )
    outer_water = openmc.Universe(cells=[openmc.Cell(fill=water)])
    cache: dict[tuple[str, int, bool], openmc.Universe] = {}
    universes = []
    for r, row in enumerate(LOADING):
        line = []
        for c, (enrichment, pyrex) in enumerate(row):
            rodded = (r, c) in BANK_D
            key = (enrichment, pyrex, rodded)
            if key not in cache:
                cache[key] = assembly(
                    m[f"fuel {enrichment}"],
                    clad,
                    water,
                    pyrex,
                    rod if rodded else None,
                    pyrex_cell,
                )
            line.append(cache[key])
        universes.append(line + [outer_water] * (len(LOADING) - len(line)))
    core = openmc.RectLattice()
    core.pitch = (ASSEMBLY_PITCH, ASSEMBLY_PITCH)
    # H8 の中心が原点。8 行が上
    core.lower_left = (
        -ASSEMBLY_PITCH / 2,
        -(len(LOADING) - 0.5) * ASSEMBLY_PITCH,
    )
    core.universes = universes
    core.outer = outer_water

    x0 = openmc.XPlane(0.0, boundary_type="reflective")
    y0 = openmc.YPlane(0.0, boundary_type="reflective")
    quadrant = +x0 & -y0
    fuel = grown(0.0)
    inner = grown(BAFFLE_GAP)
    baffle = grown(BAFFLE_GAP + BAFFLE_THICKNESS)
    cylinders = [openmc.ZCylinder(r=r) for r in (*BARREL, *PAD, LINER, *VESSEL)]
    barrel_in, barrel_out, pad_in, pad_out, liner, vessel_in, vessel_out = cylinders
    vessel_out.boundary_type = "vacuum"
    # 原点を通る角度 θ の直線 −x sinθ + y cosθ = 0 の正の側は、角度が θ より大きい側。y ≤ 0 なので、パッドは角度
    # −61°〜−29°
    near_x, near_y = (
        openmc.Plane(a=-math.sin(math.radians(-a)), b=math.cos(math.radians(-a)))
        for a in PAD_ANGLES
    )
    pad = +pad_in & -pad_out & -near_x & +near_y
    cells = [
        openmc.Cell(fill=core, region=quadrant & fuel),
        openmc.Cell(fill=water, region=quadrant & inner & ~fuel),
        openmc.Cell(fill=steel, region=quadrant & baffle & ~inner),
        openmc.Cell(fill=water, region=quadrant & -barrel_in & ~baffle),
        openmc.Cell(fill=steel, region=quadrant & +barrel_in & -barrel_out),
        openmc.Cell(fill=water, region=quadrant & +barrel_out & -liner & ~pad),
        openmc.Cell(fill=steel, region=quadrant & pad),
        openmc.Cell(fill=steel, region=quadrant & +liner & -vessel_in),
        openmc.Cell(fill=m["carbon steel"], region=quadrant & +vessel_in & -vessel_out),
    ]
    return openmc.Geometry(cells)


def pin_edges() -> tuple[list[float], list[int]]:
    """1/4 炉心の x ≥ 0 の、ピンのセルと集合体の間の余白の境界（昇順）と、ピンのセルの番号（余白を除く）。
    H 列の中央の列は x = 0 で半分に切れる。"""
    edges = [0.0]
    pins = []
    for c in range(len(LOADING)):
        centre = c * ASSEMBLY_PITCH
        for j in range(N):
            low = centre + (j - N / 2) * PITCH
            high = low + PITCH
            if high <= 0.0:
                continue
            if low > edges[-1] + 1e-9:
                # 隣の集合体との間の水のギャップ（0.08 cm）のセル
                edges.append(low)
            pins.append(len(edges) - 1)
            edges.append(high)
    return edges, pins


def pin_mesh() -> openmc.RectilinearMesh:
    edges, _ = pin_edges()
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = edges
    mesh.y_grid = [-e for e in reversed(edges)]
    # 形状は z に一様で z の境界がなく、核分裂の位置は世代とともに z に広がるので、z は十分に広くとる
    mesh.z_grid = [-1.0e6, 1.0e6]
    return mesh


def assembly_mesh() -> openmc.RegularMesh:
    """1/4 炉心の集合体の格子（H8 の中心が原点）に合わせたメッシュ。H 列と 8 行のビンは、半分が対称面の外にある。"""
    mesh = openmc.RegularMesh()
    n = len(LOADING)
    mesh.lower_left = (-ASSEMBLY_PITCH / 2, -(n - 0.5) * ASSEMBLY_PITCH)
    mesh.upper_right = ((n - 0.5) * ASSEMBLY_PITCH, ASSEMBLY_PITCH / 2)
    mesh.dimension = (n, n)
    return mesh


def read_run(source: Path) -> dict[str, Any]:
    return common_read_run(source, list(RODS), OPTIONAL | openmc_tallies.OPTIONAL)


def solve(source: Path, work: Path, options: Options) -> None:
    run = read_run(source)
    problem = run["problem"]
    work.mkdir(parents=True, exist_ok=True)
    m = {
        **{
            f"fuel {e}": material(f"fuel {e}", d, TEMPERATURE) for e, d in FUELS.items()
        },
        "clad": material("zircaloy-4", ZIRCALOY, TEMPERATURE),
        "water": material("moderator", MODERATOR_DENSE, TEMPERATURE),
        "steel": material("ss304", SS304, TEMPERATURE),
        "carbon steel": material("carbon steel", CARBON_STEEL, TEMPERATURE),
    }
    m["water"].add_s_alpha_beta("c_H_in_H2O")
    geometry = geometry_of(problem, m)
    materials = openmc.Materials(geometry.get_all_materials().values())
    materials.cross_sections = str(run["cross-sections"])

    size = len(LOADING) * ASSEMBLY_PITCH - ASSEMBLY_PITCH / 2
    settings = box_settings(
        run["particles"],
        run["batches"],
        run["inactive"],
        run["seed"],
        (0.0, -size, -1.0),
        (size, 0.0, 1.0),
    )
    # 源の収束を世代ごとに見るためのエントロピー（集合体ごと）
    settings.entropy_mesh = assembly_mesh()
    set_statepoints(settings, run)
    settings.temperature = temperature_mode(TEMPERATURE)

    # ピンと集合体。集合体も別の tally で数えるのは、集合体の中のピンの間の世代ごとの相関を含めた統計誤差を、
    # OpenMC に世代の値から出させるため（ピンの誤差の二乗和は相関を落とす）
    tallies = openmc.Tallies()
    for name, mesh in (("fission", pin_mesh()), ("assemblies", assembly_mesh())):
        tally = openmc.Tally(name=name)
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["fission"]
        tallies.append(tally)
    if openmc_tallies.enabled(run):
        tallies.extend(
            openmc_tallies.core_2d(list(materials), assembly_mesh(), HALF_GAP)
        )
        openmc_tallies.quiet(settings)
    model = openmc.Model(geometry, materials, settings, tallies)
    run_model(model, run, work, options)


def unfold(quarter: np.ndarray) -> np.ndarray:
    """右下の 1/4（[行（上から）][列（左から）]、先頭の行と列が対称面の上）を、鏡面反射で全体に広げる。"""
    rows = np.concatenate([quarter[:0:-1], quarter], axis=0)
    return np.concatenate([rows[:, :0:-1], rows], axis=1)


def extract(source: Path, work: Path, results: Path) -> None:
    run = read_run(source)
    problem = run["problem"]
    edges, pins = pin_edges()
    cells = len(edges) - 1
    path = work / f"statepoint.{run['batches']}.h5"
    with openmc.StatePoint(path) as sp:
        # run を段階にしない問題では、入力を変えても statepoint は古いまま残るので、走らせた設定と入力を突き合わせる
        ran = {
            "particles": int(sp.n_particles),
            "batches": int(sp.n_batches),
            "inactive": int(sp.n_inactive),
            "seed": int(sp.seed),
        }
        expected = {key: run[key] for key in ran}
        if ran != expected:
            raise SystemExit(f"{path} was run with {ran}, but {source} says {expected}")
        k = sp.keff
        rate, std = mesh_values(sp, "fission", (cells, cells))
        by, by_std = mesh_values(sp, "assemblies", (len(LOADING), len(LOADING)))
    results.mkdir(parents=True, exist_ok=True)
    # 余白を除き、対称面で切れたピン（先頭の行と列）を全体のピンに直す
    widths = np.diff(edges)[pins] / PITCH
    area = np.outer(widths, widths)
    rate = unfold(rate[np.ix_(pins, pins)] / area)
    std = unfold(std[np.ix_(pins, pins)] / area)
    mean = fuel_mean(rate)
    write_pin_powers(results / "pins.msgpack", rate / mean)
    write_pin_powers(results / "pins-std.msgpack", std / mean)

    # 集合体: 対称面の上の集合体は半分（中心の H8 は 1/4）だけ数えたので、全体に直す
    full = np.ones(len(LOADING))
    full[0] = 2.0
    by, by_std = (unfold(v * np.outer(full, full)) for v in (by, by_std))
    mean = fuel_mean(by)
    write_pin_powers(results / "assemblies.msgpack", by / mean)
    write_pin_powers(results / "assemblies-std.msgpack", by_std / mean)
    write_k(results / "k.json", problem, k)


if __name__ == "__main__":
    run_or_extract(solve, extract, __doc__)
