"""fission_production.py（無限希釈で数えた核分裂する核種の、落とした核分裂の生成の 1 次の見積もり）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_fission_production.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fission_production as fp


class DroppedProduction(unittest.TestCase):
    def test_the_first_order_estimate_matches_the_exact_change_of_one_over_k(self):
        # 無限媒質 k = P / A。核分裂 F を吸収と生成（ν F）の両方から外したときの Δρ = 1/k_with − 1/k_without
        absorption, k, nu = 1.0, 0.8, 3.5
        production = k * absorption
        fission = 1e-5
        exact = (
            absorption / production
            - (absorption - fission) / (production - nu * fission)
        ) * 1e5
        estimate = fp.dropped_production_pcm(fission / absorption, nu, k)
        self.assertLess(estimate, 0.0)
        self.assertAlmostEqual(estimate / exact, 1.0, places=4)

    def test_a_nuclide_with_nu_equal_to_k_changes_nothing(self):
        self.assertEqual(fp.dropped_production_pcm(1e-3, 1.1, 1.1), 0.0)


class Nubar(unittest.TestCase):
    class Tabulated:
        def __init__(self, x, y):
            self.x, self.y = x, y

        def __call__(self, e):
            for (x0, y0), (x1, y1) in zip(
                zip(self.x, self.y), zip(self.x[1:], self.y[1:])
            ):
                if x0 <= e <= x1:
                    return y0 + (y1 - y0) * (e - x0) / (x1 - x0)
            raise ValueError(e)

    def test_tabulated_nu_gives_the_thermal_value_and_the_largest_value(self):
        nu = self.Tabulated([1e-5, 1.0, 2e7], [2.4, 2.4, 5.0])
        thermal, largest = fp.nubar({"LNU": 2, "nu": nu})
        self.assertAlmostEqual(thermal, 2.4)
        self.assertAlmostEqual(largest, 5.0)

    def test_polynomial_nu_gives_the_thermal_value_and_the_largest_value(self):
        thermal, largest = fp.nubar({"LNU": 1, "C": [2.0, 1e-7]})
        self.assertAlmostEqual(thermal, 2.0 + 1e-7 * 0.0253)
        self.assertAlmostEqual(largest, 2.0 + 1e-7 * 2e7)


def table(fission_sigma):
    """燃焼度 0 と 10 の 2 点、核種 U-238 と X のテーブル。2 群の中性子束は (10, 5)、格子の吸収は 0.1 × 10 + 0.2 × 5 = 2"""
    state = {
        "total": [0.5, 1.0],
        "scatter": [[0.4, 0.05], [0.0, 0.9]],
        "chi": [1.0, 0.0],
        "absorption": [0.1, 0.2],
        "k_inf": 1.1,
    }
    return {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "points": [[0], [1]],
        "states": [state, state],
        "micro": {"nuclides": ["U-238"], "densities": [[0.5], [0.5]]},
        "cooling": {
            "nuclides": ["U-238", "X", "Y"],
            "points": [
                [
                    {"fission": 0.0, "reactions": []},
                    {"fission": s, "reactions": []},
                    {"fission": 0.0, "reactions": []},
                ]
                for s in fission_sigma
            ],
        },
    }


def history(x_density, k_inf):
    """燃焼の領域 1 つ（体積 2）。U-238 の数密度 1 で、テーブルの均質化した 0.5 から格子の体積は 4"""
    nuclides = ["U-238", "X", "Y"]
    return {
        "format": "xlc-depletion-history",
        "nuclides": nuclides,
        "transportable": ["U-238"],
        "dilute_absorbers": ["X", "Y"],
        "regions": [{"volume": 2.0}],
        "points": [
            {"burnup": b, "k_inf": k_inf, "densities": [[1.0, d, 1.0]]}
            for b, d in [(0.0, 0.0), (10.0, x_density)]
        ],
    }


class FissionFractions(unittest.TestCase):
    def test_the_fraction_is_the_homogenized_fission_rate_over_the_lattice_absorption(
        self,
    ):
        h = history(0.2, 0.9)
        (point,) = fp.fission_fractions(h, table([3.0, 3.0]), h)
        self.assertEqual(point["burnup"], 10.0)
        self.assertEqual(point["k_inf"], 0.9)
        # N̄ = 0.2 × 2 / 4 = 0.1、σf φ̄ = 3 × 15、格子の吸収 2
        self.assertAlmostEqual(point["fractions"]["X"], 0.1 * 3.0 * 15.0 / 2.0)
        # 核分裂しない吸収体は入らない
        self.assertEqual(list(point["fractions"]), ["X"])

    def test_densities_come_from_the_history_and_the_volume_from_the_one_before(
        self,
    ):
        before = history(0.0, 1.1)
        dilute = history(0.4, 0.9)
        (point,) = fp.fission_fractions(before, table([3.0, 3.0]), dilute)
        self.assertAlmostEqual(point["fractions"]["X"], 0.2 * 3.0 * 15.0 / 2.0)

    def test_a_history_without_dilute_absorbers_is_an_error(self):
        h = history(0.2, 0.9)
        del h["dilute_absorbers"]
        with self.assertRaises(SystemExit):
            fp.fission_fractions(h, table([3.0, 3.0]), h)


class Select(unittest.TestCase):
    def test_nuclides_are_added_until_the_rest_is_at_most_the_threshold(self):
        # 核種ごとに、燃料と燃焼度を通じた最大の絶対値で並べ、残りの和が閾値以下になるまで足す
        estimates = {
            "A": [("u", 10.0, -3.0), ("mox", 60.0, -40.0)],
            "B": [("u", 60.0, -4.0)],
            "C": [("gd", 60.0, -1.5)],
            "D": [("u", 60.0, -0.5)],
        }
        result = fp.select(estimates, 5.0)
        self.assertEqual(
            [r["nuclide"] for r in result["ranking"]], ["A", "B", "C", "D"]
        )
        self.assertEqual(result["ranking"][0]["fuel"], "mox")
        self.assertEqual(result["ranking"][0]["burnup"], 60.0)
        self.assertAlmostEqual(result["ranking"][0]["rest_pcm"], 6.0)
        self.assertAlmostEqual(result["ranking"][1]["rest_pcm"], 2.0)
        self.assertEqual(result["selected"], ["A", "B"])

    def test_nothing_is_added_when_the_whole_sum_is_at_most_the_threshold(self):
        result = fp.select({"A": [("u", 60.0, -2.0)], "B": [("u", 60.0, -3.0)]}, 5.0)
        self.assertEqual(result["selected"], [])


class Estimate(unittest.TestCase):
    def test_the_estimate_sums_the_nuclides_and_sits_beside_the_measured_difference(
        self,
    ):
        points = [
            {
                "burnup": 10.0,
                "k_inf": 0.9,
                "fractions": {"X": 1e-4, "Y": 2e-4, "Z": 1.0},
            }
        ]
        nus = {"X": (3.0, 5.0), "Y": (2.5, 5.0), "Z": (3.0, 5.0)}
        measured = [{"burnup": 10.0, "dilute_unshielded_pcm": -30.0}]
        (row,) = fp.estimate(points, nus, ["X", "Y"], measured)
        x = fp.dropped_production_pcm(1e-4, 3.0, 0.9)
        y = fp.dropped_production_pcm(2e-4, 2.5, 0.9)
        self.assertAlmostEqual(row["nuclides"]["X"], x)
        self.assertAlmostEqual(row["estimate_pcm"], x + y)
        self.assertEqual(row["measured_pcm"], -30.0)

    def test_a_burnup_missing_from_the_measured_comparison_is_an_error(self):
        points = [{"burnup": 10.0, "k_inf": 0.9, "fractions": {"X": 1e-4}}]
        with self.assertRaises(SystemExit):
            fp.estimate(points, {"X": (3.0, 5.0)}, ["X"], [{"burnup": 20.0}])


if __name__ == "__main__":
    unittest.main()
