"""VERA Problem 1B〜1D（ピンセル）と 2B〜2D（17×17 集合体）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py <cross_sections.xml> <問題> <作業ディレクトリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種] [xlc のライブラリ]

<問題> は 1B, 1C, 1D, 2B, 2C, 2D。断面積は openmc_library.py で作ったもの。
組成は仕様書（CASL-U-2012-0131-004）の表 P1-4, P2-4 の個数密度で、Hf（1e-6 /(barn cm) 未満）と
ギャップの He は計画のとおり除く（ギャップは真空）。被覆管・案内管・計装管と減速材は 600 K。
2B〜2D は lattice と同じ形式のピンごとの核分裂率を pins.msgpack に書く。
xlc のライブラリ（msgpack）を渡すと、その群構造で燃料の核種ごとの実効断面積
（燃料の領域の反応率 / (個数密度 × 中性子束)）を、1B〜1D では fuel-xs.msgpack に、2B〜2D ではピンごとに
pin-xs.msgpack に書く（自己遮蔽の検証用）。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

# tally の score の、statepoint に保存される名前（書いた名前と違うもの。弾性散乱は "(n,elastic)" になる）
STORED = {"elastic": "(n,elastic)"}

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
# 被覆管と減速材の温度 (K)
COLD = 600.0
FUEL_TEMPERATURE = {
    "1B": 600.0,
    "1C": 900.0,
    "1D": 1200.0,
    "2B": 600.0,
    "2C": 900.0,
    "2D": 1200.0,
}

# 個数密度 (/(barn cm))
FUEL = {
    "U234": 6.11864e-06,
    "U235": 7.18132e-04,
    "U236": 3.29861e-06,
    "U238": 2.21546e-02,
    "O16": 4.57642e-02,
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
# 0.661 g/cc、1300 ppm のホウ素
MODERATOR = {
    "O16": 2.20729e-02,
    "H1": 4.41459e-02,
    "B10": 9.52537e-06,
    "B11": 3.83408e-05,
}

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


def material(
    name: str, densities: dict[str, float], temperature: float
) -> openmc.Material:
    m = openmc.Material(name=name, temperature=temperature)
    for nuclide, density in densities.items():
        m.add_nuclide(nuclide, density, "ao")
    m.set_density("sum")
    return m


def tube(inner: float, outer: float, clad: openmc.Material, water: openmc.Material):
    """中が減速材の管（案内管・計装管）のセル。"""
    r_in, r_out = openmc.ZCylinder(r=inner), openmc.ZCylinder(r=outer)
    return openmc.Universe(
        cells=[
            openmc.Cell(fill=water, region=-r_in),
            openmc.Cell(fill=clad, region=+r_in & -r_out),
            openmc.Cell(fill=water, region=+r_out),
        ]
    )


def main() -> None:
    cross_sections = Path(sys.argv[1]).resolve()
    problem = sys.argv[2]
    work = Path(sys.argv[3]).resolve()
    particles = int(sys.argv[4]) if len(sys.argv) > 4 else 100_000
    batches = int(sys.argv[5]) if len(sys.argv) > 5 else 300
    inactive = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    # 統計誤差の見積もりを確かめるには、種を変えて解き直す
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    library = Path(sys.argv[8]).resolve() if len(sys.argv) > 8 else None
    if problem not in FUEL_TEMPERATURE:
        raise SystemExit(
            f"unknown problem {problem}; expected one of {list(FUEL_TEMPERATURE)}"
        )
    work.mkdir(parents=True, exist_ok=True)

    fuel = material("fuel", FUEL, FUEL_TEMPERATURE[problem])
    clad = material("zircaloy-4", ZIRCALOY, COLD)
    water = material("moderator", MODERATOR, COLD)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(cross_sections)

    pellet, inner, outer = (
        openmc.ZCylinder(r=r) for r in (PELLET, CLAD_INNER, CLAD_OUTER)
    )
    rod = openmc.Universe(
        cells=[
            openmc.Cell(fill=fuel, region=-pellet),
            # ギャップは真空
            openmc.Cell(region=+pellet & -inner),
            openmc.Cell(fill=clad, region=+inner & -outer),
            openmc.Cell(fill=water, region=+outer),
        ]
    )

    if problem.startswith("1"):
        half = PITCH / 2
        root = openmc.Cell(fill=rod)
    else:
        guide = tube(GUIDE_INNER, GUIDE_OUTER, clad, water)
        instrument = tube(INSTRUMENT_INNER, INSTRUMENT_OUTER, clad, water)
        universes = [[rod] * N for _ in range(N)]
        for i, j in GUIDE_TUBES:
            universes[i][j] = guide
        universes[INSTRUMENT_TUBE[0]][INSTRUMENT_TUBE[1]] = instrument
        lattice = openmc.RectLattice()
        lattice.pitch = (PITCH, PITCH)
        lattice.lower_left = (-N * PITCH / 2, -N * PITCH / 2)
        lattice.universes = universes
        lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=water)])
        half = N * PITCH / 2 + HALF_GAP
        root = openmc.Cell(fill=lattice)
    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    root.region = -box
    geometry = openmc.Geometry([root])

    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    # ライブラリにある温度だけを使う。OpenMC は温度を kT (eV) で持ち、K に直すと丸めで
    # わずかにずれるので、1 K の幅を許す（ライブラリの温度は 300 K ずつ離れている）
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
        constraints={"fissionable": True},
    )

    tallies = openmc.Tallies()
    if problem.startswith("2"):
        mesh = openmc.RegularMesh()
        mesh.dimension = (N, N)
        mesh.lower_left = (-N * PITCH / 2, -N * PITCH / 2)
        mesh.upper_right = (N * PITCH / 2, N * PITCH / 2)
        tally = openmc.Tally(name="fission")
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["fission"]
        tallies.append(tally)
        # ピンごとの燃料の核種の実効断面積（メッシュはピンのセルと一致する）
        if library is not None:
            bounds = msgpack.unpackb(library.read_bytes())["group_bounds"]
            pin_filters = [
                openmc.MeshFilter(mesh),
                openmc.MaterialFilter([fuel]),
                openmc.EnergyFilter(sorted(bounds)),
            ]
            tally = openmc.Tally(name="pin fuel")
            tally.filters = pin_filters
            tally.nuclides = list(FUEL)
            tally.scores = ["total", "elastic", "fission"]
            flux = openmc.Tally(name="pin fuel flux")
            flux.filters = pin_filters
            flux.scores = ["flux"]
            tallies.extend([tally, flux])

    # 燃料の核種ごとの実効断面積（群構造は xlc のライブラリから。OpenMC の境界は昇順）
    if problem.startswith("1") and library is not None:
        bounds = msgpack.unpackb(library.read_bytes())["group_bounds"]
        energy = openmc.EnergyFilter(sorted(bounds))
        tally = openmc.Tally(name="fuel")
        tally.filters = [openmc.MaterialFilter([fuel]), energy]
        tally.nuclides = list(FUEL)
        tally.scores = ["total", "elastic", "fission"]
        flux = openmc.Tally(name="fuel flux")
        flux.filters = [openmc.MaterialFilter([fuel]), energy]
        flux.scores = ["flux"]
        tallies.extend([tally, flux])

    model = openmc.Model(geometry, materials, settings, tallies)
    statepoint = model.run(cwd=work, output=False)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        if problem.startswith("2"):
            values = sp.get_tally(name="fission").get_values(scores=["fission"]).ravel()
            # メッシュは x が速く変わり、y は下から上。lattice と同じく上の行から並べる
            rates = values.reshape(N, N)[::-1, :]
            rates = rates / rates[rates > 0].mean()
            (work / "pins.msgpack").write_bytes(
                msgpack.packb({"format": "xlc-pin-powers", "rows": rates.tolist()})
            )
            if library is not None:
                # [メッシュ（x が速く、y は下から）][群（昇順）] を [行（上から）][列][群（高いエネルギーから）] にする
                def pins(values):
                    groups = values.size // (N * N)
                    return values.reshape(N, N, groups)[::-1, :, ::-1]

                phi = pins(
                    sp.get_tally(name="pin fuel flux").get_values(scores=["flux"])
                )
                rates = sp.get_tally(name="pin fuel")
                xs = {}
                for nuclide, density in FUEL.items():
                    xs[nuclide] = {}
                    for score in ["total", "elastic", "fission"]:
                        rate = pins(
                            rates.get_values(
                                scores=[STORED.get(score, score)], nuclides=[nuclide]
                            )
                        )
                        # 燃料のないセルは 0
                        value = np.divide(
                            rate,
                            density * phi,
                            out=np.zeros_like(rate),
                            where=phi > 0,
                        )
                        xs[nuclide][score] = value.tolist()
                (work / "pin-xs.msgpack").write_bytes(
                    msgpack.packb({"format": "xlc-pin-xs", "nuclides": xs})
                )
        if problem.startswith("1") and library is not None:
            phi = sp.get_tally(name="fuel flux").get_values(scores=["flux"]).ravel()
            rates = sp.get_tally(name="fuel")
            # 群は xlc と同じく高いエネルギーから並べる
            xs = {
                nuclide: {
                    score: (
                        rates.get_values(
                            scores=[STORED.get(score, score)], nuclides=[nuclide]
                        ).ravel()
                        / (density * phi)
                    )[::-1].tolist()
                    for score in ["total", "elastic", "fission"]
                }
                for nuclide, density in FUEL.items()
            }
            (work / "fuel-xs.msgpack").write_bytes(
                msgpack.packb({"format": "xlc-fuel-xs", "nuclides": xs})
            )
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
