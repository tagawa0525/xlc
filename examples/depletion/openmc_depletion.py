"""VERA Problem 1 のピンセル（1B〜1D）を OpenMC の燃焼（openmc.deplete）で燃やし、xlc の格子計算の燃焼（ADR-0173）と比べる参照解を作る。

使い方（references の devShell で）:

    python3 openmc_depletion.py run <入力> <作業ディレクトリ>
    python3 openmc_depletion.py extract <入力> <作業ディレクトリ> <結果（msgpack）>

<入力> は openmc-1C-rings1.toml の形の TOML で、計算の条件をすべて書く（既定値はない。ADR-0130）: 問題（`problem`、1B〜1D）、
燃料のリングの数（`rings`、等しい面積に分け、リングごとに 1 つの燃焼する物質にする）、openmc_library.py で作った
`cross-sections`、燃焼鎖（`chain`。`Chain.from_endf` の XML）、初めの重金属の出力密度（`power`、W/g）、燃焼度の点（`burnup`、
MWd/kg、0 より大きく増える列）、1 世代の粒子数（`particles`）、世代数（`batches`）、捨てる世代数（`inactive`）、乱数の種
（`seed`）、スレッド数（`threads`）。パスは入力のファイルからの相対。

方式は ADR-0173 の表のとおり xlc と揃える: CE/CM（`CECMIntegrator`、CRAM48）、一定の出力密度（`timestep_units="MWd/kg"`）、
鎖の核分裂の Q で規格化（`normalization_mode="fission-q"`）、連続エネルギーの反応率（`reaction_rate_mode="direct"`）、
収率は 112.5 eV（0.0253 eV と 500 keV の幾何平均。ADR-0164 の区間の境界）で熱と速に分ける（`fission_yield_mode="cutoff"`）。
2 次元のピンセルなので、物質の体積は高さ 1 cm あたりの面積 (cm²) にする（出力と重金属の質量はどちらも 1 cm あたり）。

run は作業ディレクトリに depletion_results.h5 を書き、輸送（`CoupledOperator.__call__`）と CRAM の時間を表示する。extract は
燃焼点ごとの k と統計誤差、リングごとの全核種の数密度 (1/(barn cm))、ステップごとの輸送と CRAM の時間を msgpack に書く
（compare_depletion.py が読む）。
"""

import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import h5py
import msgpack
import numpy as np
import openmc
import openmc.deplete

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "vera"))
from openmc_common import eigenvalue_settings, material
from openmc_input import input_path, read_input
from openmc_model import (
    CLAD_INNER,
    CLAD_OUTER,
    FUEL,
    FUEL_TEMPERATURE,
    MODERATOR_OF,
    MODERATOR_TEMPERATURE,
    PELLET,
    PITCH,
    ZIRCALOY,
    temperature_settings,
)

KEYS = {
    "problem",
    "rings",
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
PROBLEMS = ["1B", "1C", "1D"]
# 収率の熱と速の境界 (eV): ADR-0164 の区間の境界（0.0253 eV と 500 keV の幾何平均）
YIELD_CUTOFF = math.sqrt(0.0253 * 5.0e5)
# 輸送の時間（CoupledOperator.__call__ の 1 回ごと、s）
TRANSPORT_SECONDS: list[float] = []


def read(source: Path) -> dict:
    run = read_input(source, KEYS, set())
    if run["problem"] not in PROBLEMS:
        raise SystemExit(
            f"unknown problem {run['problem']}; expected one of {PROBLEMS}"
        )
    if not (isinstance(run["rings"], int) and run["rings"] > 0):
        raise SystemExit(f"rings must be a positive integer, got {run['rings']}")
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
    return [PELLET * math.sqrt((i + 1) / rings) for i in range(rings)]


def build(run: dict) -> tuple[openmc.Model, list[openmc.Material]]:
    """ピンセルのモデルと、燃焼するリングの物質（内側から）。"""
    problem = run["problem"]
    radii = ring_radii(run["rings"])
    fuels = []
    inner = 0.0
    for i, r in enumerate(radii):
        fuel = material(f"fuel-{i}", FUEL, FUEL_TEMPERATURE[problem])
        fuel.depletable = True
        fuel.volume = math.pi * (r * r - inner * inner)
        fuels.append(fuel)
        inner = r
    clad = material("zircaloy-4", ZIRCALOY, MODERATOR_TEMPERATURE[problem])
    water = material("moderator", MODERATOR_OF[problem], MODERATOR_TEMPERATURE[problem])
    water.add_s_alpha_beta("c_H_in_H2O")
    cells = []
    surface = None
    for fuel, r in zip(fuels, radii):
        boundary = openmc.ZCylinder(r=r)
        region = -boundary if surface is None else +surface & -boundary
        cells.append(openmc.Cell(fill=fuel, region=region))
        surface = boundary
    gap, outer = openmc.ZCylinder(r=CLAD_INNER), openmc.ZCylinder(r=CLAD_OUTER)
    box = openmc.model.RectangularPrism(PITCH, PITCH, boundary_type="reflective")
    cells.append(openmc.Cell(fill=None, region=+surface & -gap))
    cells.append(openmc.Cell(fill=clad, region=+gap & -outer))
    cells.append(openmc.Cell(fill=water, region=+outer & -box))
    settings = eigenvalue_settings(
        run["particles"], run["batches"], run["inactive"], run["seed"], PITCH / 2
    )
    settings.temperature = temperature_settings(problem)
    settings.output = {"summary": False, "tallies": False}
    model = openmc.Model(
        openmc.Geometry(cells), openmc.Materials([*fuels, clad, water]), settings
    )
    return model, fuels


class TimedOperator(openmc.deplete.CoupledOperator):
    """輸送の 1 回ごとの時間を測り、燃料の温度 `fuel_temperature` (K) のデータのある核種だけを燃料で輸送する
    `CoupledOperator`。

    xlc は、ライブラリが燃焼の領域の温度で持たない鎖の核種を輸送しない（崩壊だけ。ADR-0173）。OpenMC は最も近い温度の
    モードでは別の温度のデータを使うので、燃焼する物質の組成を OpenMC に渡すとき（`_update_materials`）と反応率の tally の
    核種（`_get_reaction_nuclides`）に限って、`nuclides_with_data` を燃料の温度のデータのある核種に絞り、xlc と揃える。
    `nuclides_with_data` そのものを絞ると、燃焼しない物質（減速材の H-1 など）からも核種が落ちる。
    """

    fuel_temperature: float
    # 燃料で輸送する核種（記録に、xlc と同じ集合であることを書く）
    transportable: list[str]

    def _get_nuclides_with_data(self, cross_sections):
        everything = super()._get_nuclides_with_data(cross_sections)
        kept = set()
        for library in openmc.data.DataLibrary.from_xml(cross_sections).libraries:
            if library["type"] != "neutron":
                continue
            with h5py.File(library["path"], "r") as f:
                for name in library["materials"]:
                    temperatures = [float(k.rstrip("K")) for k in f[name]["kTs"]]
                    # OpenMC は温度を kT で持ち、K に直すと丸めでずれるので 1 K の幅を許す（openmc_model.py と同じ）
                    if any(abs(t - self.fuel_temperature) <= 1.0 for t in temperatures):
                        kept.add(name)
        print(
            f"decay only in the fuel (no data at {self.fuel_temperature} K): "
            f"{sorted(everything - kept)}"
        )
        TimedOperator.transportable = sorted(kept)
        return everything

    def _in_fuel(self, method):
        """`nuclides_with_data` を燃料の温度のある核種に替えて `method` を呼ぶ。"""
        everything = self.nuclides_with_data
        self.nuclides_with_data = set(TimedOperator.transportable)
        try:
            return method()
        finally:
            self.nuclides_with_data = everything

    def _update_materials(self):
        return self._in_fuel(super()._update_materials)

    def _get_reaction_nuclides(self):
        return self._in_fuel(super()._get_reaction_nuclides)

    def __call__(self, vec, source_rate):
        start = time.perf_counter()
        result = super().__call__(vec, source_rate)
        TRANSPORT_SECONDS.append(time.perf_counter() - start)
        return result


def solve(source: Path, work: Path) -> None:
    run = read(source)
    work.mkdir(parents=True, exist_ok=True)
    os.chdir(work)
    os.environ["OMP_NUM_THREADS"] = str(run["threads"])
    openmc.config["cross_sections"] = run["cross-sections"]
    model, _ = build(run)
    TimedOperator.fuel_temperature = FUEL_TEMPERATURE[run["problem"]]
    operator = TimedOperator(
        model,
        chain_file=run["chain"],
        normalization_mode="fission-q",
        fission_yield_mode="cutoff",
        fission_yield_opts={"cutoff": YIELD_CUTOFF},
        reaction_rate_mode="direct",
    )
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
