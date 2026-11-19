"""estimates.py（#0293 の項目 6、7、3 の見積もり）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_estimates.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimates


class Flux(unittest.TestCase):
    def test_two_group_infinite_medium_gives_the_analytic_k(self):
        # 1 群目から 2 群目へ 0.02、吸収 0.01 と 0.1、生成は 2 群目だけ 0.15: k = 0.15 × 0.02 / (0.03 × 0.1)
        state = {
            "total": [0.53, 0.5],
            "scatter": [[0.5, 0.02], [0.0, 0.4]],
            "chi": [1.0, 0.0],
            "nu_fission": [0.0, 0.15],
        }
        phi, k = estimates.flux(state)
        self.assertAlmostEqual(k, 0.15 * 0.02 / (0.03 * 0.1))
        self.assertAlmostEqual(phi[1] / phi[0], 0.02 / 0.1)


class ShieldingFactor(unittest.TestCase):
    def test_the_table_is_interpolated_in_log_sigma0_and_divided_by_infinite_dilution(
        self,
    ):
        sigma0 = [1e10, 1e4, 1e2]
        table = [10.0, 8.0, 6.0]
        self.assertAlmostEqual(estimates.shielding_factor(table, sigma0, 1e3), 0.7)
        self.assertEqual(estimates.shielding_factor(table, sigma0, 1e11), 1.0)
        self.assertEqual(estimates.shielding_factor([0.0, 0.0, 0.0], sigma0, 1e3), 1.0)


class UrrGroups(unittest.TestCase):
    def test_only_groups_inside_both_the_ufg_and_the_unresolved_range_count(self):
        bounds = [1e7, 24780.0, 15030.0, 9118.0, 5530.0, 3519.1, 2239.45, 1425.1]
        self.assertEqual(estimates.urr_groups(bounds, (2250.0, 25000.0)), [2, 3, 4])
        self.assertEqual(
            estimates.urr_groups(bounds, (1400.0, 25000.0)), [2, 3, 4, 5, 6]
        )
        self.assertEqual(estimates.urr_groups(bounds, (20000.0, 1.5e5)), [])


class CumulativeYields(unittest.TestCase):
    def test_independent_yields_follow_the_decay_branches(self):
        chain = {
            "nuclides": [
                {"decay": {"modes": [{"ratio": 0.8, "products": [{"nuclide": 1, "count": 1}]},
                                     {"ratio": 0.2, "products": [{"nuclide": 2, "count": 1}]}]}},
                {"decay": {"modes": [{"ratio": 1.0, "products": [{"nuclide": 2, "count": 1}]}]}},
                {"decay": None},
            ]
        }  # fmt: skip
        cumulative = estimates.cumulative_yields(chain, {0: 0.1, 1: 0.01, 2: 0.001})
        self.assertAlmostEqual(cumulative[1], 0.01 + 0.08)
        self.assertAlmostEqual(cumulative[2], 0.001 + 0.02 + 0.09)


class YieldPoints(unittest.TestCase):
    def test_nuclides_without_both_points_use_one_yield(self):
        self.assertEqual(
            estimates.yield_points({"energies": [0.0253, 5e5, 1.4e7]}), (0, 1)
        )
        self.assertIsNone(estimates.yield_points({"energies": [5e5, 1.4e7]}))


if __name__ == "__main__":
    unittest.main()
