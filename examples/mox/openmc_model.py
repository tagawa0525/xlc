"""LWR MOX ピンセル（MOX-1, MOX-2）を OpenMC の連続エネルギーモードで解く。

使い方: python3 openmc_model.py run <入力> <作業ディレクトリ>
        python3 openmc_model.py extract <入力> <作業ディレクトリ> <結果のディレクトリ>

<入力> は openmc-mox1.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 問題（`problem`。MOX-1 か
MOX-2）、examples/vera/openmc_library.py で作った `cross-sections`（cross_sections.xml）、1 世代の粒子数（`particles`）、
世代数（`batches`）、捨てる世代数（`inactive`）、乱数の種（`seed`）、スレッド数（`threads`）、群構造を読む xlc の
ライブラリ（`effective-xs-groups`）。パスは入力のファイルからの相対。
run は OpenMC を解いて作業ディレクトリに statepoint を残し、extract は statepoint から k を k.json に、
燃料の実効断面積を fuel-xs.msgpack に結果のディレクトリへ書く（ADR-0123）。
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
    input_path,
    material,
    read_input,
)
from openmc_run import RUN_KEYS, run_or_extract, write_k

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


def read_run(source: Path) -> dict:
    """入力の TOML を読み、問題を確かめ、パスを絶対パスにする。"""
    run = read_input(source, RUN_KEYS | {"effective-xs-groups"}, set())
    if run["problem"] not in FUEL:
        raise SystemExit(
            f"unknown problem {run['problem']}; expected one of {list(FUEL)}"
        )
    run["cross-sections"] = input_path(source, run["cross-sections"])
    run["effective-xs-groups"] = input_path(source, run["effective-xs-groups"])
    return run


def solve(source: Path, work: Path) -> None:
    run = read_run(source)
    densities = FUEL[run["problem"]]
    work.mkdir(parents=True, exist_ok=True)

    fuel = material("fuel", densities, TEMPERATURE)
    clad = material("zirconium", ZIRCONIUM, TEMPERATURE)
    water = material("water", WATER, TEMPERATURE)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(run["cross-sections"])

    pellet, outer = openmc.ZCylinder(r=PELLET), openmc.ZCylinder(r=CLAD_OUTER)
    box = openmc.model.RectangularPrism(PITCH, PITCH, boundary_type="reflective")
    geometry = openmc.Geometry(
        [
            openmc.Cell(fill=fuel, region=-pellet),
            openmc.Cell(fill=clad, region=+pellet & -outer),
            openmc.Cell(fill=water, region=+outer & -box),
        ]
    )

    settings = eigenvalue_settings(
        run["particles"], run["batches"], run["inactive"], run["seed"], PITCH / 2
    )
    # ライブラリの温度だけを使う（kT から K に直すときの丸めのために 1 K の幅を許す）
    settings.temperature = {"method": "nearest", "tolerance": 1.0}

    # 燃料の核種ごとの実効断面積（群構造は xlc のライブラリから）
    filters = [
        openmc.MaterialFilter([fuel]),
        group_filter(run["effective-xs-groups"]),
    ]
    tallies = fuel_tallies("fuel", filters, list(densities), SCORES)

    # 核種の吸収率の差を k の差に換算するための、セル全体の吸収と生成
    balance = openmc.Tally(name="balance")
    balance.scores = ["absorption", "nu-fission"]
    tallies.append(balance)

    model = openmc.Model(geometry, materials, settings, openmc.Tallies(tallies))
    model.run(cwd=work, threads=run["threads"], output=False)


def extract(source: Path, work: Path, results: Path) -> None:
    run = read_run(source)
    results.mkdir(parents=True, exist_ok=True)
    statepoint = work / f"statepoint.{run['batches']}.h5"
    write_fuel_xs(statepoint, results / "fuel-xs.msgpack", FUEL[run["problem"]])
    with openmc.StatePoint(statepoint) as sp:
        write_k(results / "k.json", run["problem"], sp.keff)


if __name__ == "__main__":
    run_or_extract(solve, extract, __doc__)
