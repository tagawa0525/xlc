"""scatter_shares.py の散乱源の割合（#0180）。

実行: cd examples/depletion/issue-0180 && python3 -m unittest test_scatter_shares
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scatter_shares import shares


class Shares(unittest.TestCase):
    def test_constant_shares_are_the_ratio_of_the_scattering(self):
        # 一定の断面積では、どの割合も N σs の比になる
        grid = np.geomspace(4.0, 1000.0, 50)
        scatter = {"a": np.full_like(grid, 3.0), "b": np.full_like(grid, 1.0)}
        total = np.full_like(grid, 8.0)
        out = shares(grid, scatter, total, 4.0, 1000.0)
        for value in out["a"]:
            self.assertAlmostEqual(value, 0.75, places=12)
        for value in out["b"]:
            self.assertAlmostEqual(value, 0.25, places=12)

    def test_narrow_resonance_weight_lowers_the_share_where_the_total_is_large(self):
        # 全断面積の大きい半分（レサジーで）の散乱は、NR の中性子束 1/Σt で軽くなり、平均より NR の割合が小さい
        grid = np.geomspace(4.0, 1000.0, 2001)
        upper = grid > np.sqrt(4.0 * 1000.0)
        scatter = {"a": np.where(upper, 1.0, 0.0), "b": np.ones_like(grid)}
        total = np.where(upper, 100.0, 2.0)
        mean, nr, peak = shares(grid, scatter, total, 4.0, 1000.0)["a"]
        self.assertTrue(0.2 < mean < 0.3)
        self.assertLess(nr, 0.1 * mean)
        self.assertEqual(peak, 0.5)


if __name__ == "__main__":
    unittest.main()
