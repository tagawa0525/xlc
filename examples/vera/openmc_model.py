"""VERA Problem 1B〜1D（ピンセル）と 2B〜2D（17×17 集合体）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py <cross_sections.xml> <問題> <作業ディレクトリ>
        [粒子数/世代] [世代数] [非アクティブバッチ数] [乱数のシード] [xlc のライブラリ]

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


def geometry_of(
    problem: str,
    fuel: openmc.Material,
    clad: openmc.Material,
    water: openmc.Material,
    gap: openmc.Material | None = None,
) -> tuple[openmc.Geometry, float]:
    """問題の形状と、外周の箱の半分の幅 (cm)。外周は鏡面反射。`gap` はペレットと被覆管の間（None なら真空）。"""
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
    """コマンドラインの 2〜7 番目の引数: <問題> <作業ディレクトリ> [粒子数/世代] [世代数] [非アクティブバッチ数] [乱数のシード]。"""

    def __init__(self, argv: list[str]) -> None:
        self.problem = argv[2]
        if self.problem not in FUEL_TEMPERATURE:
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
        """世代の数と乱数のシード、燃料に一様な初期の中性子源（`half` は外周の箱の半分の幅）。"""
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
    water = material("moderator", MODERATOR, COLD)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(cross_sections)

    geometry, half = geometry_of(problem, fuel, clad, water)

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
