"""burned_composition.py の体積の平均（#0274）。

実行: cd examples/depletion/issue-0274/fine-pin && python3 -m unittest test_burned_composition
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import burned_composition as bc


def history() -> dict:
    return {
        "format": "xlc-depletion-history",
        "nuclides": ["O-16", "U-235", "Pu-239", "B-10"],
        "transportable": ["O-16", "U-235", "Pu-239", "B-10"],
        "regions": [
            {"material": "fuel-0", "volume": 1.0},
            {"material": "fuel-1", "volume": 3.0},
            {"material": "water", "volume": 10.0},
        ],
        "points": [
            {"burnup": 0.0, "densities": [[1, 2, 0, 0], [1, 2, 0, 0], [0, 0, 0, 5]]},
            {"burnup": 60.0, "densities": [[1, 1, 4, 0], [1, 2, 0, 0], [0, 0, 0, 5]]},
        ],
    }


class AverageTest(unittest.TestCase):
    def test_densities_are_weighted_by_the_volume_of_the_matching_regions(self):
        got = bc.average(history(), 60.0, "fuel", ["U-235", "Pu-239"])
        self.assertAlmostEqual(got["U-235"], (1 * 1 + 3 * 2) / 4)
        self.assertAlmostEqual(got["Pu-239"], (1 * 4) / 4)

    def test_regions_of_other_materials_are_left_out(self):
        got = bc.average(history(), 60.0, "fuel", ["B-10"])
        self.assertEqual(got["B-10"], 0.0)

    def test_a_burnup_off_the_history_points_is_an_error(self):
        with self.assertRaises(SystemExit):
            bc.average(history(), 30.0, "fuel", ["U-235"])

    def test_a_nuclide_the_history_does_not_transport_is_an_error(self):
        with self.assertRaises(SystemExit):
            bc.average(history(), 60.0, "fuel", ["Am-241"])

    def test_a_prefix_matching_no_region_is_an_error(self):
        with self.assertRaises(SystemExit):
            bc.average(history(), 60.0, "pellet", ["U-235"])


if __name__ == "__main__":
    unittest.main()
