"""compare_steps.py（#0293 の項目 2 の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_compare_steps.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_steps


def history(points):
    """燃焼の領域 2 つ（fuel と gd0、体積 1 と 3）。`points` は (燃焼度, k, [[領域ごとの核種の数密度]])"""
    return {
        "nuclides": ["U-235", "Gd-155"],
        "regions": [
            {"material": "fuel", "volume": 1.0},
            {"material": "gd0", "volume": 3.0},
        ],
        "points": [{"burnup": b, "k_inf": k, "densities": d} for b, k, d in points],
    }


class Compare(unittest.TestCase):
    def test_the_fine_history_is_compared_at_the_coarse_points_only(self):
        coarse = history(
            [(0.0, 1.2, [[1.0, 0.0], [1.0, 2.0]]), (2.0, 1.1, [[1.0, 0.0], [1.0, 2.0]])]
        )
        fine = history(
            [
                (0.0, 1.2, [[1.0, 0.0], [1.0, 2.0]]),
                (1.0, 1.15, [[9.0, 0.0], [9.0, 9.0]]),
                (2.0, 1.1001, [[1.1, 0.0], [1.0, 2.2]]),
            ]
        )
        rows = compare_steps.compare(
            coarse, fine, {"fuel": ["U-235"], "gd": ["Gd-155"]}
        )
        self.assertEqual([r["burnup"] for r in rows], [0.0, 2.0])
        self.assertAlmostEqual(rows[1]["dk_pcm"], 10.0, places=6)
        self.assertAlmostEqual(rows[1]["densities"]["fuel:U-235"], 0.1)
        self.assertAlmostEqual(rows[1]["densities"]["gd:Gd-155"], 0.1)

    def test_densities_are_averaged_over_the_regions_by_volume(self):
        h = history([(0.0, 1.0, [[1.0, 0.0], [2.0, 0.0]])])
        h["regions"][0]["material"] = "gd7"
        self.assertAlmostEqual(
            compare_steps.averaged(h, h["points"][0], "gd", "U-235"), 1.75
        )

    def test_a_missing_point_in_the_fine_history_is_an_error(self):
        coarse = history([(2.0, 1.1, [[1.0, 0.0], [1.0, 2.0]])])
        fine = history([(1.0, 1.1, [[1.0, 0.0], [1.0, 2.0]])])
        with self.assertRaises(SystemExit):
            compare_steps.compare(coarse, fine, {"fuel": ["U-235"]})


if __name__ == "__main__":
    unittest.main()
