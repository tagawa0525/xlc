"""LWR MOX ピンセル（MOX-1, MOX-2）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py <cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [非アクティブバッチ数] [乱数のシード]

<問題> は MOX-1, MOX-2。断面積は examples/vera/openmc_library.py で作ったもの。
組成は pin_cell.pdf の付録 2 の個数密度で、セルは同じ面積の正方形（ピッチ 1.26 cm、外周は鏡面反射）、
温度はすべて 600 K にする（README）。xlc のライブラリ（msgpack）の群構造で、燃料の核種ごとの実効断面積
（燃料の領域の反応率 / (個数密度 × 中性子束)）を fuel-xs.msgpack に書く（examples/vera/openmc_model.py と
同じ形式）。比べやすいよう、燃料の中性子束 `flux`（体積で割らない tally の値）、反応率の相対標準偏差
`<score>-rsd`、同じ規格化でのセル全体の吸収率 `absorption` と ν × 核分裂率 `nu-fission` も書く。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

# 書き出す名前と tally の score の名前（statepoint では弾性散乱は "(n,elastic)" になる）。
# capture は統計誤差を見るためで、比べる捕獲断面積は全 − 弾性 − 核分裂から求める
STORED = {"elastic": "(n,elastic)", "capture": "(n,gamma)"}
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


def material(name: str, densities: dict[str, float]) -> openmc.Material:
    m = openmc.Material(name=name, temperature=TEMPERATURE)
    for nuclide, density in densities.items():
        m.add_nuclide(nuclide, density, "ao")
    m.set_density("sum")
    return m


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
    fuel = material("fuel", densities)
    clad = material("zirconium", ZIRCONIUM)
    water = material("water", WATER)
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

    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    half = PITCH / 2
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
        constraints={"fissionable": True},
    )
    # ライブラリの温度だけを使う（kT から K に直すときの丸めのために 1 K の幅を許す）
    settings.temperature = {"method": "nearest", "tolerance": 1.0}

    # 燃料の核種ごとの実効断面積（群構造は xlc のライブラリから。OpenMC の境界は昇順）
    bounds = msgpack.unpackb(library.read_bytes())["group_bounds"]
    filters = [openmc.MaterialFilter([fuel]), openmc.EnergyFilter(sorted(bounds))]
    rates = openmc.Tally(name="fuel")
    rates.filters = filters
    rates.nuclides = list(densities)
    rates.scores = [STORED.get(score, score) for score in SCORES]
    flux = openmc.Tally(name="fuel flux")
    flux.filters = filters
    flux.scores = ["flux"]

    # 核種の吸収率の差を k の差に換算するための、セル全体の吸収と生成
    balance = openmc.Tally(name="balance")
    balance.scores = ["absorption", "nu-fission"]

    model = openmc.Model(
        geometry, materials, settings, openmc.Tallies([rates, flux, balance])
    )
    statepoint = model.run(cwd=work, output=False)
    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        phi = sp.get_tally(name="fuel flux").get_values(scores=["flux"]).ravel()
        tally = sp.get_tally(name="fuel")
        xs = {}
        for nuclide, density in densities.items():
            xs[nuclide] = {}
            for score in SCORES:
                name = STORED.get(score, score)
                mean = tally.get_values(scores=[name], nuclides=[nuclide]).ravel()
                std = tally.get_values(
                    scores=[name], nuclides=[nuclide], value="std_dev"
                ).ravel()
                rsd = np.divide(std, mean, out=np.zeros_like(mean), where=mean > 0)
                # 群は xlc と同じく高いエネルギーから並べる
                xs[nuclide][score] = (mean / (density * phi))[::-1].tolist()
                xs[nuclide][f"{score}-rsd"] = rsd[::-1].tolist()
        (work / "fuel-xs.msgpack").write_bytes(
            msgpack.packb(
                {
                    "format": "xlc-fuel-xs",
                    "keff": [k.nominal_value, k.std_dev],
                    "flux": phi[::-1].tolist(),
                    **{
                        score: float(
                            sp.get_tally(name="balance")
                            .get_values(scores=[score])
                            .ravel()[0]
                        )
                        for score in ["absorption", "nu-fission"]
                    },
                    "densities": densities,
                    "nuclides": xs,
                }
            )
        )
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


if __name__ == "__main__":
    main()
