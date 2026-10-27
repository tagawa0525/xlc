"""LWR MOX ピンセル（examples/mox の MOX-1、MOX-2）を OpenMC の燃焼（openmc.deplete）で燃焼させ、xlc の格子計算の燃焼と比べる
参照解を作る（ADR-0206）。

使い方（references の devShell で）:

    python3 openmc_depletion_mox.py run <入力> <作業ディレクトリ>
    python3 openmc_depletion_mox.py extract <入力> <作業ディレクトリ> <結果（msgpack）>

<入力> は openmc-mox1-rings1.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 問題（`problem`、MOX-1 か
MOX-2）、燃料のリングの数（`rings`、等しい面積に分け、リングごとに 1 つの燃焼する物質にする）、燃料の温度（`fuel-temperature`、K）、
被覆管と減速材の温度（`moderator-temperature`、K）、openmc_library.py で作った `cross-sections`、燃焼チェーン（`chain`）、初めの
重金属の出力密度（`power`、W/g）、燃焼度の点（`burnup`、MWd/kg、0 より大きく増える列）、1 世代の粒子数（`particles`）、世代数
（`batches`）、非アクティブバッチ数（`inactive`）、乱数のシード（`seed`）、スレッド数（`threads`）。パスは入力のファイルからの相対。

形状と組成は examples/mox/openmc_model.py の問題のもの（ギャップのない被覆管、ホウ素のない減速材）。燃焼の方式、輸送する核種の
規則、run と extract の出力の形は openmc_depletion.py と同じ（ADR-0173、ADR-0183）で、extract の msgpack は compare_depletion.py が読む。
`read`、`ring_radii`、`solve`、`extract` は openmc_depletion.py の同じ名前の関数のコピー（問題、形状、温度の扱いだけが違う）。
共通の関数にすると、そちらを deps にする openmc-1C が古くなるのでコピーしている（#0132）。
"""

import importlib.util
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import msgpack
import numpy as np
import openmc
import openmc.deplete

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from openmc_common import eigenvalue_settings, material
from openmc_depletion import TRANSPORT_SECONDS, YIELD_CUTOFF, TimedOperator
from openmc_input import input_path, read_input


def load_mox_model():
    """examples/mox/openmc_model.py（形状と組成）。openmc_depletion.py が import する examples/vera/openmc_model.py と
    モジュールの名前が同じなので、sys.path でなく別の名前で読み込む。"""
    spec = importlib.util.spec_from_file_location(
        "mox_openmc_model", HERE.parent / "mox" / "openmc_model.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MOX = load_mox_model()

KEYS = {
    "problem",
    "rings",
    "fuel-temperature",
    "moderator-temperature",
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


def read(source: Path) -> dict:
    run = read_input(source, KEYS, set())
    if run["problem"] not in MOX.FUEL:
        raise SystemExit(
            f"unknown problem {run['problem']}; expected one of {list(MOX.FUEL)}"
        )
    if not (isinstance(run["rings"], int) and run["rings"] > 0):
        raise SystemExit(f"rings must be a positive integer, got {run['rings']}")
    for key in ("fuel-temperature", "moderator-temperature"):
        if not run[key] > 0:
            raise SystemExit(f"{key} must be positive, got {run[key]}")
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
    return run


def ring_radii(rings: int) -> list[float]:
    """ペレットを等しい面積に分けたリングの外半径 (cm)。"""
    return [MOX.PELLET * math.sqrt((i + 1) / rings) for i in range(rings)]


def build(run: dict) -> tuple[openmc.Model, list[openmc.Material]]:
    """ピンセルのモデルと、燃焼するリングの物質（内側から）。"""
    radii = ring_radii(run["rings"])
    hot, cold = float(run["fuel-temperature"]), float(run["moderator-temperature"])
    fuels = []
    inner = 0.0
    for i, r in enumerate(radii):
        fuel = material(f"fuel-{i}", MOX.FUEL[run["problem"]], hot)
        fuel.depletable = True
        fuel.volume = math.pi * (r * r - inner * inner)
        fuels.append(fuel)
        inner = r
    clad = material("zirconium", MOX.ZIRCONIUM, cold)
    water = material("water", MOX.WATER, cold)
    water.add_s_alpha_beta("c_H_in_H2O")
    cells = []
    surface = None
    for fuel, r in zip(fuels, radii):
        boundary = openmc.ZCylinder(r=r)
        region = -boundary if surface is None else +surface & -boundary
        cells.append(openmc.Cell(fill=fuel, region=region))
        surface = boundary
    outer = openmc.ZCylinder(r=MOX.CLAD_OUTER)
    box = openmc.model.RectangularPrism(
        MOX.PITCH, MOX.PITCH, boundary_type="reflective"
    )
    cells.append(openmc.Cell(fill=clad, region=+surface & -outer))
    cells.append(openmc.Cell(fill=water, region=+outer & -box))
    settings = eigenvalue_settings(
        run["particles"], run["batches"], run["inactive"], run["seed"], MOX.PITCH / 2
    )
    # ライブラリにある温度の表だけを使う（kT から K に直すときの丸めのために 1 K の幅を許す）。燃焼のライブラリは入力の温度
    # そのものの表を持ち、内挿の要る温度（VERA の 565 K の S(α,β)）を持たない
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    settings.output = {"summary": False, "tallies": False}
    model = openmc.Model(
        openmc.Geometry(cells), openmc.Materials([*fuels, clad, water]), settings
    )
    return model, fuels


def check_transportable(run: dict, transportable: list[str]) -> None:
    """初めの燃料の核種がすべて燃料で輸送される（ライブラリが燃料の温度で持つ）ことを確かめる。欠けると TimedOperator は
    その核種を崩壊だけにし、輸送を解き終えてから比較の検査で分かるので、輸送の前に止める。"""
    missing = sorted(set(MOX.FUEL[run["problem"]]) - set(transportable))
    if missing:
        raise SystemExit(
            f"the library has no data at {run['fuel-temperature']} K for the fuel "
            f"nuclides {missing}"
        )


def solve(source: Path, work: Path) -> None:
    run = read(source)
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    os.environ["OMP_NUM_THREADS"] = str(run["threads"])
    openmc.config["cross_sections"] = run["cross-sections"]
    model, _ = build(run)
    TimedOperator.fuel_temperature = float(run["fuel-temperature"])
    operator = TimedOperator(
        model,
        chain_file=run["chain"],
        normalization_mode="fission-q",
        fission_yield_mode="cutoff",
        fission_yield_opts={"cutoff": YIELD_CUTOFF},
        reaction_rate_mode="direct",
    )
    check_transportable(run, TimedOperator.transportable)
    burnup = run["burnup"]
    steps = [b - a for a, b in itertools.pairwise([0.0, *burnup])]
    integrator = openmc.deplete.CECMIntegrator(
        operator,
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
    print(
        f"transport {sum(TRANSPORT_SECONDS):.1f} s in {len(TRANSPORT_SECONDS)} solves, "
        f"total {total:.1f} s"
    )


def extract(source: Path, work: Path, output: Path) -> None:
    run = read(source)
    results = openmc.deplete.Results(work / "depletion_results.h5")
    _, fuels = build(run)
    times, k = results.get_keff(time_units="s")
    nuclides = list(results[0].index_nuc)
    rings = []
    for fuel in fuels:
        mat = str(fuel.id)
        densities = {
            name: results.get_atoms(mat, name, nuc_units="atom/b-cm")[1].tolist()
            for name in nuclides
        }
        rings.append({"name": fuel.name, "volume": fuel.volume, "densities": densities})
    seconds = json.loads((work / "transport-seconds.json").read_text())
    output.write_bytes(
        msgpack.packb(
            {
                "format": "openmc-depletion",
                "burnup": [0.0, *run["burnup"]],
                "time": np.asarray(times).tolist(),
                "k": k[:, 0].tolist(),
                "k-std": k[:, 1].tolist(),
                "rings": rings,
                "transport-seconds": seconds["transport"],
                "total-seconds": seconds["total"],
                "cram-seconds": [float(r.proc_time) for r in results],
                "transportable": seconds["transportable"],
            }
        )
    )
    for b, kk in zip([0.0, *run["burnup"]], k):
        print(f"burnup {b} k {kk[0]:.6f} +- {kk[1]:.6f}")


if __name__ == "__main__":
    match sys.argv[1:]:
        case ["run", source, work]:
            solve(Path(source).resolve(), Path(work).resolve())
        case ["extract", source, work, output]:
            extract(
                Path(source).resolve(), Path(work).resolve(), Path(output).resolve()
            )
        case _:
            raise SystemExit(__doc__)
