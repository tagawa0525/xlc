"""added_effect.py（全核種を輸送した燃焼と試行の燃焼の k∞ の差を、核種の吸収の変化に分ける）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_added_effect.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import added_effect


def capture(sigma):
    return {"fission": 0.0, "reactions": [{"mt": 102, "channels": [sigma]}]}


# 2 群の状態: 散乱が無く、核分裂の中性子はすべて速い群に出るので、中性子束は [1, 0]、P = νΣf,1 = 2、A = Σa,1 = 1
STATE = {
    "k_inf": 2.0,
    "total": [1.0, 1.0],
    "scatter": [[0.0, 0.0], [0.0, 0.0]],
    "chi": [1.0, 0.0],
    "absorption": [1.0, 0.0],
    "nu_fission": [2.0, 0.0],
}
NUCLIDES = ["U-235", "Sm-149", "Pm-148"]


def history(densities, transportable):
    """燃焼度 0 と 10 の点を持つ、体積 0.5 の燃焼の領域 1 つの履歴"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": NUCLIDES,
        "transportable": transportable,
        "regions": [{"volume": 0.5}],
        "points": [
            {"burnup": 0.0, "densities": [[0.0, 0.0, 0.0]]},
            {"burnup": 10.0, "densities": [densities]},
        ],
    }


def table(k10):
    state10 = {**STATE, "k_inf": k10}
    return {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "points": [[0], [1]],
        "states": [STATE, state10],
        "cooling": {
            "nuclides": NUCLIDES,
            "points": [
                [capture(0.0), capture(0.0), capture(0.0)],
                [capture(10.0), capture(100.0), capture(1000.0)],
            ],
        },
    }


class Contributions(unittest.TestCase):
    def setUp(self):
        self.before = history([1e-3, 1e-6, 2e-7], ["U-235", "Sm-149"])
        self.after = history([1e-3, 3e-6, 1e-7], ["U-235", "Sm-149", "Pm-148"])
        self.points = added_effect.contributions(
            self.before, self.after, table(1.5), table(1.4), 1.0, ["Pm-148"]
        )

    def test_the_burnup_zero_point_is_left_out(self):
        self.assertEqual([p["burnup"] for p in self.points], [10.0])

    def test_the_direct_difference_comes_from_the_two_tables(self):
        self.assertAlmostEqual(self.points[0]["direct_pcm"], 1e5 * (1 / 1.5 - 1 / 1.4))

    def test_a_nuclide_of_both_burnups_counts_its_density_change(self):
        # 均質化した数密度は体積の比 0.5 / 1.0 を掛けたもの。Δρ = −ΔN̄ σa φ̄ / P
        sm = self.points[0]["nuclides"]["Sm-149"]
        self.assertAlmostEqual(sm["pcm"], -1e5 * (3e-6 - 1e-6) * 0.5 * 100.0 / 2.0)
        self.assertAlmostEqual(sm["density_ratio"], 3.0)

    def test_an_added_nuclide_counts_its_whole_absorption(self):
        # 試行の燃焼は足した核種を輸送しないので、その吸収は 0 から数える
        pm = self.points[0]["nuclides"]["Pm-148"]
        self.assertAlmostEqual(pm["pcm"], -1e5 * 1e-7 * 0.5 * 1000.0 / 2.0)
        self.assertTrue(pm["added"])

    def test_the_terms_are_summed_by_kind_and_the_rest_is_the_remainder(self):
        p = self.points[0]
        nuclides = p["nuclides"]
        self.assertAlmostEqual(p["added_pcm"], nuclides["Pm-148"]["pcm"])
        self.assertAlmostEqual(p["carried_pcm"], nuclides["Sm-149"]["pcm"])
        self.assertAlmostEqual(p["heavy_pcm"], nuclides["U-235"]["pcm"])
        self.assertAlmostEqual(
            p["remainder_pcm"],
            p["direct_pcm"] - p["added_pcm"] - p["carried_pcm"] - p["heavy_pcm"],
        )

    def test_the_added_list_must_match_the_gained_transported_nuclides(self):
        with self.assertRaises(SystemExit):
            added_effect.contributions(
                self.before, self.after, table(1.5), table(1.4), 1.0, []
            )

    def test_absorption_beyond_the_lattice_total_is_an_error(self):
        # 格子の体積を小さくしすぎると、燃料の核種の吸収が格子全体の吸収を超える
        with self.assertRaises(SystemExit):
            added_effect.contributions(
                self.before, self.after, table(1.5), table(1.4), 1e-9, ["Pm-148"]
            )


class Shares(unittest.TestCase):
    def test_the_added_nuclides_get_their_lattice_absorption_fraction(self):
        # 吸収の割合の JSON の形（dilution_shielding.py が読む）。割合は N̄ σa φ̄ / A（A = Σa,1 φ1 = 1）、全核種の燃焼が
        # 自分の吸収で減らした数密度なので、数密度の比は 1
        before = history([1e-3, 1e-6, 2e-7], ["U-235", "Sm-149"])
        after = history([1e-3, 3e-6, 1e-7], ["U-235", "Sm-149", "Pm-148"])
        points = added_effect.contributions(
            before, after, table(1.5), table(1.4), 1.0, ["Pm-148"]
        )
        shares = added_effect.shares(points)
        (point,) = shares["points"]
        self.assertEqual(point["burnup"], 10.0)
        self.assertEqual(point["k_inf"], 1.4)
        self.assertEqual(list(point["untransported_all"]), ["Pm-148"])
        pm = point["untransported_all"]["Pm-148"]
        self.assertAlmostEqual(pm["fraction"], 1e-7 * 0.5 * 1000.0 / 1.0)
        self.assertEqual(pm["captured_ratio"], 1.0)


class Heavy(unittest.TestCase):
    def test_actinides_are_heavy_and_fission_products_are_not(self):
        self.assertTrue(added_effect.heavy("Pu-239"))
        self.assertTrue(added_effect.heavy("Am-242m"))
        self.assertFalse(added_effect.heavy("Pm-148"))


class LatticeVolume(unittest.TestCase):
    def test_the_volume_is_the_pitch_squared_per_cell(self):
        text = 'pitch = 1.26\nrows = ["pin pin pin", "pin gd  pin", "pin pin pin"]\n'
        self.assertAlmostEqual(added_effect.lattice_volume(text), 9 * 1.26**2)


class Report(unittest.TestCase):
    def test_the_table_has_a_row_per_fuel_and_burnup(self):
        summary = {
            "fuel": "u",
            "points": [
                {
                    "burnup": 60.0,
                    "direct_pcm": -188.0,
                    "added_pcm": -40.0,
                    "carried_pcm": -150.0,
                    "heavy_pcm": 1.0,
                    "remainder_pcm": 1.0,
                    "nuclides": {},
                }
            ],
        }
        rows = added_effect.terms_table([summary]).splitlines()[2:]
        self.assertEqual(len(rows), 1)
        self.assertIn("U", rows[0])
        self.assertIn("−188.0", rows[0])


if __name__ == "__main__":
    unittest.main()
