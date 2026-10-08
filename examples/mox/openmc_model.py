"""LWR MOX ピンセル（MOX-1, MOX-2）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py <cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

<問題> は MOX-1, MOX-2。断面積は examples/vera/openmc_library.py で作ったもの。
組成は pin_cell.pdf の付録 2 の個数密度で、セルは同じ面積の正方形（ピッチ 1.26 cm、外周は鏡面反射）、
温度はすべて 600 K にする（verification.md）。xlc のライブラリ（msgpack）の群構造で、燃料の核種ごとの実効断面積
（燃料の領域の反応率 / (個数密度 × 中性子束)）を fuel-xs.msgpack に書く（examples/vera/openmc_model.py と
同じ形式に捕獲 `capture` を足したもの）。比べやすいよう、燃料の中性子束 `flux`（体積で割らない tally の値）と
その相対標準偏差 `flux-rsd`、反応率の相対標準偏差 `<score>-rsd`、同じ規格化でのセル全体の吸収率 `absorption` と
ν × 核分裂率 `nu-fission` も書く。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import (
    STORED,
    effective_xs,
    eigenvalue_settings,
    fuel_tallies,
    group_filter,
    material,
)

SCORES = ["total", "elastic", "fission", "capture"]

PITCH = 1.26
PELLET = 0.41
CLAD_OUTER = 0.475
TEMPERATURE = 600.0

# 個数密度 (/(barn cm))
FUEL = {
    "MOX-1": {
        "U235": 5.105e-5,
        "U238": 2.037e-2,
        "Pu238": 4.669e-5,
        "Pu239": 1.465e-3,
        "Pu240": 5.691e-4,
        "Pu241": 2.713e-4,
        "Pu242": 1.413e-4,
        "Am241": 3.028e-5,
        "O16": 4.588e-2,
    },
    "MOX-2": {
        "U235": 5.118e-5,
        "U238": 2.042e-2,
        "Pu238": 2.714e-5,
        "Pu239": 1.972e-3,
        "Pu240": 4.256e-4,
        "Pu241": 3.577e-5,
        "Pu242": 1.234e-5,
        "Am241": 1.406e-5,
        "O16": 4.588e-2,
    },
}
# 天然の Zr を同位体の存在比（at%）で分ける
ZIRCONIUM = {
    f"Zr{a}": 0.0388 * f / 100
    for a, f in {90: 51.45, 91: 11.22, 92: 17.15, 94: 17.38, 96: 2.80}.items()
}
WATER = {"H1": 0.04744, "O16": 0.02372}


def relative_std(tally: openmc.Tally, **selection) -> np.ndarray:
    """tally の値の相対標準偏差（値が 0 のところは 0）を、群の高いエネルギーから並べる。"""
    mean = tally.get_values(**selection).ravel()
    std = tally.get_values(value="std_dev", **selection).ravel()
    return np.divide(std, mean, out=np.zeros_like(mean), where=mean > 0)[::-1]


def write_fuel_xs(statepoint: Path, path: Path, densities: dict[str, float]) -> None:
    """statepoint から燃料の実効断面積を取り出して `path` に書く。"""
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        # 群は xlc と同じく高いエネルギーから並べる
        xs = effective_xs(sp, "fuel", densities, SCORES, lambda v: v.ravel()[::-1])
        rates = sp.get_tally(name="fuel")
        flux = sp.get_tally(name="fuel flux")
        balance = sp.get_tally(name="balance")
        nuclides = {}
        for nuclide, by in xs.items():
            nuclides[nuclide] = {score: value.tolist() for score, value in by.items()}
            for score in SCORES:
                rsd = relative_std(
                    rates, scores=[STORED.get(score, score)], nuclides=[nuclide]
                )
                nuclides[nuclide][f"{score}-rsd"] = rsd.tolist()
        data = {
            "format": "xlc-fuel-xs",
            "keff": [k.nominal_value, k.std_dev],
            "flux": flux.get_values(scores=["flux"]).ravel()[::-1].tolist(),
            "flux-rsd": relative_std(flux, scores=["flux"]).tolist(),
            **{
                score: float(balance.get_values(scores=[score]).ravel()[0])
                for score in ["absorption", "nu-fission"]
            },
            "densities": densities,
            "nuclides": nuclides,
        }
    path.write_bytes(msgpack.packb(data))


def main() -> None:
    cross_sections = Path(sys.argv[1]).resolve()
    problem = sys.argv[2]
    if problem not in FUEL:
        raise SystemExit(f"unknown problem {problem}; expected one of {list(FUEL)}")
    work = Path(sys.argv[3]).resolve()
    library = Path(sys.argv[4]).resolve()
    particles = int(sys.argv[5]) if len(sys.argv) > 5 else 100_000
    batches = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    inactive = int(sys.argv[7]) if len(sys.argv) > 7 else 50
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    work.mkdir(parents=True, exist_ok=True)

    densities = FUEL[problem]
    fuel = material("fuel", densities, TEMPERATURE)
    clad = material("zirconium", ZIRCONIUM, TEMPERATURE)
    water = material("water", WATER, TEMPERATURE)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(cross_sections)

    pellet, outer = openmc.ZCylinder(r=PELLET), openmc.ZCylinder(r=CLAD_OUTER)
    box = openmc.model.RectangularPrism(PITCH, PITCH, boundary_type="reflective")
    geometry = openmc.Geometry(
        [
            openmc.Cell(fill=fuel, region=-pellet),
            openmc.Cell(fill=clad, region=+pellet & -outer),
            openmc.Cell(fill=water, region=+outer & -box),
        ]
    )

    settings = eigenvalue_settings(particles, batches, inactive, seed, PITCH / 2)
    # ライブラリの温度だけを使う（kT から K に直すときの丸めのために 1 K の幅を許す）
    settings.temperature = {"method": "nearest", "tolerance": 1.0}

    # 燃料の核種ごとの実効断面積（群構造は xlc のライブラリから）
    filters = [openmc.MaterialFilter([fuel]), group_filter(library)]
    tallies = fuel_tallies("fuel", filters, list(densities), SCORES)

    # 核種の吸収率の差を k の差に換算するための、セル全体の吸収と生成
    balance = openmc.Tally(name="balance")
    balance.scores = ["absorption", "nu-fission"]
    tallies.append(balance)

    model = openmc.Model(geometry, materials, settings, openmc.Tallies(tallies))
    statepoint = model.run(cwd=work, output=False)
    write_fuel_xs(statepoint, work / "fuel-xs.msgpack", densities)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
