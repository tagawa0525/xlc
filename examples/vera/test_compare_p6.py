"""compare_p6.py の、core の標準出力の読み方と、Problem 6 の確認の量。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_p6.py'
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_p6

HFP = """state hfp
k-eff 1.16000000
leakage per fission neutron 0.0123
axial fission fractions (bottom first: z from, z to (cm), fraction)
0.000000 50.000000 0.00000000e0
50.000000 150.000000 4.00000000e-1
150.000000 250.000000 6.00000000e-1
250.000000 300.000000 0.00000000e0
assembly fission rates (relative to the assembly mean, top row is +y)
(`.` is a reflector site, `-` a position without a node)
1.0000
thermal-hydraulics iterations 7 change 3.000e-6
fuel temperature volume average 900.0000 K max 1000.0000 K
moderator density interpolation error max 2.000e-4 g/cc
channel outlets (assembly position x y from the -x -y corner: rise (kJ/kg), outlet (K))
0 0 205.540000 600.500000
thermal layers (bottom first: z from, z to (cm), moderator (K), fuel (K) or - without fission)
0.000000 50.000000 565.0000 -
50.000000 150.000000 575.0000 950.0000
150.000000 250.000000 595.0000 850.0000
250.000000 300.000000 600.5000 -
"""

LOW = HFP.replace("k-eff 1.16000000", "k-eff 1.20000100").replace(
    "iterations 7", "iterations 2"
)

ISOTHERMAL = """state inlet
k-eff 1.20000000
leakage per fission neutron 0.0123
axial fission fractions (bottom first: z from, z to (cm), fraction)
0.000000 300.000000 1.00000000e0
assembly fission rates (relative to the assembly mean, top row is +y)
(`.` is a reflector site, `-` a position without a node)
1.0000
"""

LATTICE = """[branches]
boron = { ppm = [1300], materials = ["water"] }
moderator = { temperature = [565, 570, 580, 590, 600, 610], density = [0.742713, 0.732936, 0.711901, 0.688335, 0.661165, 0.628290], materials = ["water"] }
fuel = { temperature = [565, 600], materials = ["fuel"] }
"""

CORE = """cycle = "c.toml"
[thermal-hydraulics]
power = 17.67
flow = 85.9684
pressure = 15.513204
tolerance = 1.0e-05
[state.hfp]
inlet = 565
"""


class ParseCore(unittest.TestCase):
    def test_reads_the_eigenvalue_the_thermal_summary_and_the_layers(self):
        c = compare_p6.parse_core(HFP)
        self.assertEqual(c["k"], 1.16)
        t = c["thermal"]
        self.assertEqual(t["iterations"], 7)
        self.assertEqual(t["change"], 3e-6)
        self.assertEqual((t["fuel_average"], t["fuel_max"]), (900.0, 1000.0))
        self.assertEqual(t["density_error"], 2e-4)
        self.assertEqual(
            t["channels"], [{"position": [0, 0], "rise": 205.54, "outlet": 600.5}]
        )
        self.assertEqual(
            t["layers"][:2],
            [
                {"z": [0.0, 50.0], "moderator": 565.0, "fuel": None},
                {"z": [50.0, 150.0], "moderator": 575.0, "fuel": 950.0},
            ],
        )
        self.assertEqual(c["fractions"][1], {"z": [50.0, 150.0], "fraction": 0.4})

    def test_a_solve_without_thermal_hydraulics_has_no_thermal_summary(self):
        self.assertIsNone(compare_p6.parse_core(ISOTHERMAL)["thermal"])

    def test_missing_eigenvalue_is_an_error(self):
        with self.assertRaises(SystemExit):
            compare_p6.parse_core(HFP.replace("k-eff 1.16000000\n", ""))


class TableDensity(unittest.TestCase):
    GRID = ((570.0, 0.732936), (580.0, 0.711901))

    def test_interpolates_the_grid_densities_linearly_in_temperature(self):
        self.assertAlmostEqual(
            compare_p6.table_density(self.GRID, 575.0),
            (0.732936 + 0.711901) / 2,
            places=12,
        )
        self.assertEqual(compare_p6.table_density(self.GRID, 580.0), 0.711901)

    def test_temperatures_outside_the_grid_are_errors(self):
        for t in (569.9, 580.1):
            with self.assertRaises(SystemExit):
                compare_p6.table_density(self.GRID, t)


class Summary(unittest.TestCase):
    def summary(self, hfp: str = HFP) -> dict:
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for name, text in [
                ("lattice.toml", LATTICE),
                ("core.toml", CORE),
                ("hfp.out", hfp),
                ("low.out", LOW),
                ("iso.out", ISOTHERMAL),
            ]:
                p = Path(d) / name
                p.write_text(text)
                paths.append(str(p))
            return json.loads(json.dumps(compare_p6.compare(paths)))

    def test_compares_the_rise_with_the_power_over_the_flow(self):
        s = self.summary()
        self.assertAlmostEqual(s["channel"]["expected_rise"], 17.67e3 / 85.9684, 12)
        self.assertAlmostEqual(
            s["channel"]["rise_relative_diff"], 205.54 / (17.67e3 / 85.9684) - 1, 12
        )
        self.assertEqual(s["channel"]["outlet"], 600.5)

    def test_reports_convergence_against_the_input_tolerance(self):
        s = self.summary()
        self.assertEqual(s["thermal"]["tolerance"], 1e-5)
        self.assertTrue(s["thermal"]["converged"])
        s = self.summary(HFP.replace("change 3.000e-6", "change 2.000e-5"))
        self.assertFalse(s["thermal"]["converged"])

    def test_compares_the_low_power_solve_with_the_isothermal_solve_in_pcm(self):
        s = self.summary()
        self.assertAlmostEqual(s["low_power_vs_isothermal_pcm"], 0.1, places=9)
        self.assertEqual(s["low_power"]["iterations"], 2)

    def test_layers_carry_the_table_density_and_the_fission_fraction(self):
        layer = self.summary()["layers"][1]
        self.assertAlmostEqual(layer["density"], (0.732936 + 0.711901) / 2, 12)
        self.assertEqual(layer["fuel"], 950.0)
        self.assertEqual(layer["fraction"], 0.4)

    def test_records_the_inputs_with_their_md5(self):
        s = self.summary()
        self.assertEqual(len(s["inputs"]), 5)
        self.assertEqual(len(s["inputs"][0]["md5"]), 32)

    def test_report_names_the_inputs_first(self):
        lines = compare_p6.report(self.summary())
        self.assertTrue(lines[0].startswith("inputs: "))


if __name__ == "__main__":
    unittest.main()
