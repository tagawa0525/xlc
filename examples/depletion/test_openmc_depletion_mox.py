"""openmc_depletion_mox.py の、MOX のピンセルの燃焼の入力とモデルの組み立て（ADR-0206）。

OpenMC のある環境で実行する（モデルを組み立てるだけで、輸送は解かない）:
nix develop .#references -c python3 -m unittest discover -s examples/depletion -p 'test_openmc_depletion_mox.py'
"""

import math
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import openmc_depletion_mox as mox

PELLET_AREA = math.pi * 0.41**2
# examples/mox/openmc_model.py の MOX-1（pin_cell.pdf の付録 2）
MOX1 = {
    "U235": 5.105e-5,
    "U238": 2.037e-2,
    "Pu238": 4.669e-5,
    "Pu239": 1.465e-3,
    "Pu240": 5.691e-4,
    "Pu241": 2.713e-4,
    "Pu242": 1.413e-4,
    "Am241": 3.028e-5,
    "O16": 4.588e-2,
}


def run_of(problem: str = "MOX-1", **extra) -> dict:
    """入力を読んだ形（モデルの組み立てに要るものだけ）。"""
    return {
        "problem": problem,
        "rings": 1,
        "fuel-temperature": 900.0,
        "moderator-temperature": 600.0,
        "particles": 1000,
        "batches": 10,
        "inactive": 5,
        "seed": 1,
        **extra,
    }


INPUT = """
problem = "MOX-1"
rings = 1
fuel-temperature = 900
moderator-temperature = 600
cross-sections = "openmc-library/cross_sections.xml"
chain = "openmc-chain.xml"
power = 38.41
burnup = [0.1, 0.5]
particles = 1000
batches = 10
inactive = 5
seed = 1
threads = 1
"""


def read_text(text: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "openmc-mox1.toml"
        source.write_text(text)
        return mox.read(source)


class Read(unittest.TestCase):
    def test_the_example_input_reads(self):
        run = read_text(INPUT)
        self.assertEqual(run["problem"], "MOX-1")
        self.assertEqual(run["fuel-temperature"], 900)

    def test_an_unknown_problem_is_an_error(self):
        with self.assertRaises(SystemExit):
            read_text(INPUT.replace('"MOX-1"', '"1C"'))

    def test_the_temperatures_are_required(self):
        with self.assertRaises(SystemExit):
            read_text(INPUT.replace("fuel-temperature = 900\n", ""))

    def test_rings_must_be_positive(self):
        with self.assertRaises(SystemExit):
            read_text(INPUT.replace("rings = 1", "rings = 0"))


class Transportable(unittest.TestCase):
    def test_initial_fuel_nuclides_without_data_at_the_fuel_temperature_are_an_error(
        self,
    ):
        transportable = sorted(set(MOX1) - {"Pu242"})
        with self.assertRaises(SystemExit):
            mox.check_transportable(run_of(), transportable)

    def test_initial_fuel_nuclides_with_data_pass(self):
        mox.check_transportable(run_of(), sorted(MOX1) + ["Xe135"])


class Build(unittest.TestCase):
    def test_the_rings_split_the_pellet_into_equal_areas(self):
        _, fuels = mox.build(run_of(rings=4))
        self.assertEqual(len(fuels), 4)
        for fuel in fuels:
            self.assertAlmostEqual(fuel.volume, PELLET_AREA / 4, places=12)

    def test_the_fuel_has_the_mox1_composition_at_the_fuel_temperature(self):
        _, fuels = mox.build(run_of())
        (fuel,) = fuels
        self.assertTrue(fuel.depletable)
        self.assertEqual(fuel.temperature, 900.0)
        self.assertEqual({n.name: n.percent for n in fuel.nuclides}, MOX1)

    def test_the_clad_touches_the_pellet_and_the_water_has_h_in_h2o(self):
        model, _ = mox.build(run_of())
        cells = list(model.geometry.get_all_cells().values())
        filled = [c.fill.name for c in cells if c.fill is not None]
        # ギャップはない（pin_cell.pdf の付録 2）: 空のセルを作らない
        self.assertEqual(filled, ["fuel-0", "zirconium", "water"])
        water = next(c.fill for c in cells if c.fill and c.fill.name == "water")
        self.assertEqual(water.temperature, 600.0)
        self.assertEqual([s[0] for s in water._sab], ["c_H_in_H2O"])

    def test_the_temperature_mode_uses_only_the_library_temperatures(self):
        model, _ = mox.build(run_of(**{"moderator-temperature": 565.0}))
        self.assertEqual(
            model.settings.temperature, {"method": "nearest", "tolerance": 1.0}
        )


if __name__ == "__main__":
    unittest.main()
