"""dilution_shielding.py（組成の希釈での足す核種の自己遮蔽の目安）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_dilution_shielding.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dilution_shielding as ds


def nuclide(tables):
    """1 温度（565 K）の、σ0 ごとの群の捕獲の表を持つ xs のファイル"""
    return {
        "nuclides": [
            {
                "name": "X",
                "temperatures": [
                    {
                        "temperature": 565.0,
                        "shielded": [
                            {"capture": t, "fission": [0.0] * len(t)} for t in tables
                        ],
                    }
                ],
            }
        ]
    }


class Interpolate(unittest.TestCase):
    def test_values_interpolate_in_the_logarithm_of_the_background(self):
        sigma0 = [1e4, 1e2, 1.0]
        self.assertAlmostEqual(ds.interpolate(sigma0, [10.0, 6.0, 2.0], 1e3), 8.0)
        self.assertEqual(ds.interpolate(sigma0, [10.0, 6.0, 2.0], 1e6), 10.0)
        self.assertEqual(ds.interpolate(sigma0, [10.0, 6.0, 2.0], 0.1), 2.0)


class Weights(unittest.TestCase):
    def test_each_range_sums_to_its_two_group_flux(self):
        w = ds.weights([1e3, 10.0, 1.0, 0.1, 1e-5], 565.0, 3.0, 2.0)
        self.assertAlmostEqual(w[0] + w[1], 3.0)
        self.assertAlmostEqual(w[2] + w[3], 2.0)
        # 熱外は単位のレサジーあたり一定
        self.assertAlmostEqual(w[0] / w[1], 2.0)


class Ratio(unittest.TestCase):
    def test_a_dilute_nuclide_keeps_its_infinite_dilution_cross_section(self):
        x = nuclide([[100.0, 1.0], [50.0, 1.0]])
        r, one = ds.ratio(x, [1.0, 1.0], 1e-12, [1e10, 1.0], [1.0, 1.0], 565.0)
        self.assertAlmostEqual(r, 1.0, places=6)
        self.assertAlmostEqual(one, 50.5)

    def test_a_concentrated_nuclide_takes_the_shielded_table(self):
        x = nuclide([[100.0, 1.0], [50.0, 1.0]])
        r, _ = ds.ratio(x, [1.0, 1.0], 1.0, [1e10, 1.0], [1.0, 1.0], 565.0)
        self.assertAlmostEqual(r, 51.0 / 101.0)


class Shares(unittest.TestCase):
    def test_a_dropped_nuclide_counts_as_zero_and_an_unknown_one_is_an_error(self):
        point = {"untransported_all": {"A": {"fraction": 1e-3, "captured_ratio": 0.5}}}
        pruned = {"P"}
        self.assertEqual(ds.share(point, "A", pruned)["fraction"], 1e-3)
        self.assertEqual(ds.share(point, "P", pruned)["fraction"], 0.0)
        with self.assertRaises(SystemExit):
            ds.share(point, "Z", pruned)


class Grids(unittest.TestCase):
    def test_files_with_another_sigma0_grid_or_groups_are_errors(self):
        a = {"sigma0": [1e10, 1.0], "group_bounds": [10.0, 1.0, 0.1]}
        ds.check_grid(a, a, "X")
        with self.assertRaises(SystemExit):
            ds.check_grid(a, {**a, "sigma0": [1e10, 10.0]}, "X")
        with self.assertRaises(SystemExit):
            ds.check_grid(a, {**a, "group_bounds": [10.0, 0.1]}, "X")


if __name__ == "__main__":
    unittest.main()
