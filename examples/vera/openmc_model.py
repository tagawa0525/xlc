"""VERA Problem 1B〜1D（ピンセル）と 2B〜2H（17×17 集合体）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py <cross_sections.xml> <問題> <作業ディレクトリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種] [xlc のライブラリ]

<問題> は 1B, 1C, 1D, 2B, 2C, 2D, 2E, 2E-empty（2E の Pyrex の棒を抜いたもの。吸収体と減速材の密度の効果を分けるため）, 2F,
2G, 2H。断面積は openmc_library.py で作ったもの。
組成は仕様書（CASL-U-2012-0131-004）の表 P1-4, P2-4 の個数密度で、Hf（1e-6 /(barn cm) 未満）と
ギャップの He は ADR-0054 のとおり除く（ギャップは真空。Pyrex の棒と制御棒の中の He も同じ。ADR-0120）。燃料以外は 600 K。
2B〜2H は lattice と同じ形式のピンごとの核分裂率を pins.msgpack に書く。
xlc のライブラリ（msgpack）を渡すと、その群構造で燃料の核種ごとの実効断面積
（燃料の領域の反応率 / (個数密度 × 中性子束)）を、ピンセル（1 で始まる問題）では fuel-xs.msgpack に、集合体（2 で始まる
問題）ではピンごとに pin-xs.msgpack に書く（自己遮蔽の検証用）。
"""

import sys
from pathlib import Path

import msgpack
import openmc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import (
    effective_xs,
    eigenvalue_settings,
    fuel_tallies,
    group_filter,
    material,
)

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
# 被覆管と減速材の温度 (K)
COLD = 600.0
FUEL_TEMPERATURE = {
    "1B": 600.0,
    "1C": 900.0,
    "1D": 1200.0,
    "2B": 600.0,
    "2C": 900.0,
    "2D": 1200.0,
    "2E": 600.0,
    "2E-empty": 600.0,
    "2F": 600.0,
    "2G": 600.0,
    "2H": 600.0,
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
# 1300 ppm のホウ素。0.661 g/cc（1B〜1D、2B〜2D）と 0.743 g/cc（2E〜2H）。仕様書の表 P2-1、P2-4
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
    "2B": MODERATOR,
    "2C": MODERATOR,
    "2D": MODERATOR,
    "2E": MODERATOR_DENSE,
    "2E-empty": MODERATOR_DENSE,
    "2F": MODERATOR_DENSE,
    "2G": MODERATOR_DENSE,
    "2H": MODERATOR_DENSE,
}
# Pyrex の棒（仕様書の表 6）: 内管の内と外、Pyrex の内と外、被覆管の内と外の半径。内管の中と間は真空
PYREX_RADII = (0.214, 0.231, 0.241, 0.427, 0.437, 0.484)
# 制御棒（仕様書の表 8）: SS304 の被覆管の内と外の半径と、吸収体の半径。吸収体と被覆管の間は真空
ROD_CLAD = (0.386, 0.484)
ROD_ABSORBER = {"AIC": (AIC, 0.382), "B4C": (B4C, 0.373)}

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
# 案内管に入れるもの（Pyrex または制御棒の吸収体の名前）と位置。Pyrex の 12 本は仕様書の図 5 の配置
INSERTS = {
    "2E": (
        "Pyrex",
        [
            (2, 5),
            (2, 11),
            (5, 2),
            (5, 8),
            (5, 14),
            (8, 5),
            (8, 11),
            (11, 2),
            (11, 8),
            (11, 14),
            (14, 5),
            (14, 11),
        ],
    ),
    # 24 本はすべての案内管（図 5）
    "2F": ("Pyrex", GUIDE_TUBES),
    "2G": ("AIC", GUIDE_TUBES),
    "2H": ("B4C", GUIDE_TUBES),
}


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


def insert_tube(
    layers: list[tuple[openmc.Material | None, float]],
    clad: openmc.Material,
    water: openmc.Material,
) -> openmc.Universe:
    """案内管に入れたもののセル。`layers` は中心から外へ (物質, 外側の半径) の並び（None は真空）で、その外から案内管の
    内側までは減速材、案内管の被覆管、その外は減速材。"""
    cells = []
    inner = None
    for fill, radius in [*layers, (water, GUIDE_INNER), (clad, GUIDE_OUTER)]:
        outer = openmc.ZCylinder(r=radius)
        region = -outer if inner is None else +inner & -outer
        cells.append(openmc.Cell(fill=fill, region=region))
        inner = outer
    cells.append(openmc.Cell(fill=water, region=+inner))
    return openmc.Universe(cells=cells)


def pyrex_layers(
    pyrex: openmc.Material, steel: openmc.Material
) -> list[tuple[openmc.Material | None, float]]:
    """Pyrex の棒（内管、Pyrex、被覆管。内管の中と間は真空）の層。"""
    fills = [None, steel, None, pyrex, None, steel]
    return list(zip(fills, PYREX_RADII))


def rod_layers(
    absorber: openmc.Material, radius: float, steel: openmc.Material
) -> list[tuple[openmc.Material | None, float]]:
    """制御棒（半径 `radius` の吸収体と SS304 の被覆管。間は真空）の層。"""
    clad_in, clad_out = ROD_CLAD
    return [(absorber, radius), (None, clad_in), (steel, clad_out)]


def geometry_of(
    problem: str,
    fuel: openmc.Material,
    clad: openmc.Material,
    water: openmc.Material,
    gap: openmc.Material | None = None,
    insert: openmc.Universe | None = None,
) -> tuple[openmc.Geometry, float]:
    """問題の形状と、外周の箱の半分の幅 (cm)。外周は鏡面反射。`gap` はペレットと被覆管の間（None なら真空）。
    `insert` は案内管に入れるもののセル。INSERTS にある問題で要る。"""
    pellet, inner, outer = (
        openmc.ZCylinder(r=r) for r in (PELLET, CLAD_INNER, CLAD_OUTER)
    )
    rod = openmc.Universe(
        cells=[
            openmc.Cell(fill=fuel, region=-pellet),
            openmc.Cell(fill=gap, region=+pellet & -inner),
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
        if problem in INSERTS:
            if insert is None:
                raise SystemExit(f"{problem} needs the cell of its guide tube inserts")
            for i, j in INSERTS[problem][1]:
                if (i, j) not in GUIDE_TUBES:
                    raise SystemExit(f"insert at {(i, j)} is not in a guide tube")
                universes[i][j] = insert
        lattice = openmc.RectLattice()
        lattice.pitch = (PITCH, PITCH)
        lattice.lower_left = (-N * PITCH / 2, -N * PITCH / 2)
        lattice.universes = universes
        lattice.outer = openmc.Universe(cells=[openmc.Cell(fill=water)])
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


def write_pins(sp: openmc.StatePoint, path: Path) -> None:
    """名前が fission の tally（`pin_mesh` ごと）を、lattice と同じ形式のピンの核分裂率（平均 1）にして書く。"""
    values = sp.get_tally(name="fission").get_values(scores=["fission"]).ravel()
    # メッシュは x が速く変わり、y は下から上。lattice と同じく上の行から並べる
    rates = values.reshape(N, N)[::-1, :]
    rates = rates / rates[rates > 0].mean()
    path.write_bytes(
        msgpack.packb({"format": "xlc-pin-powers", "rows": rates.tolist()})
    )


class Run:
    """コマンドラインの 2〜7 番目の引数: <問題> <作業ディレクトリ> [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]。"""

    def __init__(self, argv: list[str]) -> None:
        self.problem = argv[2]
        if self.problem not in FUEL_TEMPERATURE or self.problem not in MODERATOR_OF:
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


def main() -> None:
    cross_sections = Path(sys.argv[1]).resolve()
    run = Run(sys.argv)
    problem, work = run.problem, run.work
    library = Path(sys.argv[8]).resolve() if len(sys.argv) > 8 else None
    work.mkdir(parents=True, exist_ok=True)

    fuel = material("fuel", FUEL, FUEL_TEMPERATURE[problem])
    clad = material("zircaloy-4", ZIRCALOY, COLD)
    water = material("moderator", MODERATOR_OF[problem], COLD)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    insert = None
    if problem in INSERTS:
        kind = INSERTS[problem][0]
        steel = material("ss304", SS304, COLD)
        if kind == "Pyrex":
            pyrex = material("pyrex", PYREX, COLD)
            materials.extend([pyrex, steel])
            insert = insert_tube(pyrex_layers(pyrex, steel), clad, water)
        else:
            nuclides, radius = ROD_ABSORBER[kind]
            absorber = material(kind.lower(), nuclides, COLD)
            materials.extend([absorber, steel])
            insert = insert_tube(rod_layers(absorber, radius, steel), clad, water)
    materials.cross_sections = str(cross_sections)

    geometry, half = geometry_of(problem, fuel, clad, water, insert=insert)

    settings = run.settings(half)
    # ライブラリにある温度だけを使う。OpenMC は温度を kT (eV) で持ち、K に直すと丸めで
    # わずかにずれるので、1 K の幅を許す（ライブラリの温度は 300 K ずつ離れている）
    settings.temperature = {"method": "nearest", "tolerance": 1.0}

    tallies = openmc.Tallies()
    if problem.startswith("2"):
        mesh = pin_mesh()
        tally = openmc.Tally(name="fission")
        tally.filters = [openmc.MeshFilter(mesh)]
        tally.scores = ["fission"]
        tallies.append(tally)
        # ピンごとの燃料の核種の実効断面積（メッシュはピンのセルと一致する）
        if library is not None:
            pin_filters = [
                openmc.MeshFilter(mesh),
                openmc.MaterialFilter([fuel]),
                group_filter(library),
            ]
            tallies.extend(fuel_tallies("pin fuel", pin_filters, list(FUEL), SCORES))

    # 燃料の核種ごとの実効断面積（群構造は xlc のライブラリから。OpenMC の境界は昇順）
    if problem.startswith("1") and library is not None:
        filters = [openmc.MaterialFilter([fuel]), group_filter(library)]
        tallies.extend(fuel_tallies("fuel", filters, list(FUEL), SCORES))

    model = openmc.Model(geometry, materials, settings, tallies)
    statepoint = model.run(cwd=work, output=False)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        if problem.startswith("2"):
            write_pins(sp, work / "pins.msgpack")
            if library is not None:
                # [メッシュ（x が速く、y は下から）][群（昇順）] を [行（上から）][列][群（高いエネルギーから）] にする
                def pins(values):
                    groups = values.size // (N * N)
                    return values.reshape(N, N, groups)[::-1, :, ::-1]

                # 燃料のないセルは 0
                xs = {
                    nuclide: {score: value.tolist() for score, value in by.items()}
                    for nuclide, by in effective_xs(
                        sp, "pin fuel", FUEL, SCORES, pins
                    ).items()
                }
                (work / "pin-xs.msgpack").write_bytes(
                    msgpack.packb({"format": "xlc-pin-xs", "nuclides": xs})
                )
        if problem.startswith("1") and library is not None:
            # 群は xlc と同じく高いエネルギーから並べる
            xs = {
                nuclide: {score: value.tolist() for score, value in by.items()}
                for nuclide, by in effective_xs(
                    sp, "fuel", FUEL, SCORES, lambda v: v.ravel()[::-1]
                ).items()
            }
            (work / "fuel-xs.msgpack").write_bytes(
                msgpack.packb({"format": "xlc-fuel-xs", "nuclides": xs})
            )
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
