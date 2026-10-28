"""make_hybrid.py（ミクロ燃焼で追う元素だけを条件を変えた履歴の値にした履歴）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_make_hybrid.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_hybrid

NUCLIDES = ["O-16", "U-235", "Am-242m", "Xe-135", "Cs-133", "Gd-155", "Ag-109"]


def history(scale: float, moc_regions: list[int]) -> dict:
    """燃焼点 2 つ × 領域 2 つ。数密度は (点, 領域, 核種) と scale から決まる"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": NUCLIDES,
        "regions": [
            {"material": "fuel", "cell": 20, "regions": moc_regions, "volume": 2.0},
            {"material": "fuel", "cell": 21, "regions": moc_regions, "volume": 4.0},
        ],
        "transportable": NUCLIDES,
        "points": [
            {
                "burnup": burnup,
                "k_inf": 1.0 + scale,
                "densities": [
                    [scale * (100 * p + 10 * r + n) for n in range(len(NUCLIDES))]
                    for r in range(2)
                ],
            }
            for p, burnup in enumerate([0.0, 10.0])
        ],
    }


class Hybrid(unittest.TestCase):
    def test_tracked_elements_take_the_varied_densities_and_the_rest_stay_nominal(self):
        nominal, varied = history(1.0, [1, 2]), history(2.0, [1, 2])
        result = make_hybrid.hybrid(copy.deepcopy(nominal), varied)
        tracked = {"U-235", "Am-242m", "Xe-135", "Gd-155"}
        for point, n, v in zip(result["points"], nominal["points"], varied["points"]):
            for densities, n_row, v_row in zip(
                point["densities"], n["densities"], v["densities"]
            ):
                for nuclide, value, n_value, v_value in zip(
                    NUCLIDES, densities, n_row, v_row
                ):
                    self.assertEqual(
                        value, v_value if nuclide in tracked else n_value, nuclide
                    )

    def test_fields_other_than_the_densities_stay_nominal(self):
        result = make_hybrid.hybrid(history(1.0, [1, 2]), history(2.0, [1, 2]))
        self.assertEqual([p["k_inf"] for p in result["points"]], [2.0, 2.0])

    def test_moc_region_numbers_may_differ_between_lattices(self):
        result = make_hybrid.hybrid(history(1.0, [1, 2]), history(2.0, [5, 6, 7]))
        self.assertEqual(result["regions"][0]["regions"], [1, 2])

    def test_histories_with_other_depletion_regions_are_errors(self):
        varied = history(2.0, [1, 2])
        varied["regions"][1]["volume"] = 4.5
        with self.assertRaises(SystemExit):
            make_hybrid.hybrid(history(1.0, [1, 2]), varied)

    def test_histories_with_other_burnup_points_are_errors(self):
        varied = history(2.0, [1, 2])
        varied["points"][1]["burnup"] = 12.0
        with self.assertRaises(SystemExit):
            make_hybrid.hybrid(history(1.0, [1, 2]), varied)

    def test_histories_with_other_nuclides_are_errors(self):
        varied = history(2.0, [1, 2])
        varied["nuclides"] = NUCLIDES[:-1] + ["Ag-107"]
        with self.assertRaises(SystemExit):
            make_hybrid.hybrid(history(1.0, [1, 2]), varied)


if __name__ == "__main__":
    unittest.main()
