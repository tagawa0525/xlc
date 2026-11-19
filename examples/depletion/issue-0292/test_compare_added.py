"""compare_added.py（輸送する核種を足して燃焼させ直した k∞ の差と 1 次の目安）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_added.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_added


def history(ks):
    return {
        "format": "xlc-depletion-history",
        "points": [{"burnup": b, "k_inf": k} for b, k in ks],
    }


def absorption(points):
    return {
        "points": [
            {
                "burnup": b,
                "k_inf": k,
                "untransported_all": {
                    n: {"fraction": f, "captured_ratio": r}
                    for n, (f, r) in shares.items()
                },
            }
            for b, k, shares in points
        ]
    }


class Compare(unittest.TestCase):
    def test_the_direct_difference_and_both_estimates_sum_the_added_nuclides(self):
        before = history([(0.0, 1.2), (10.0, 1.0), (20.0, 0.9)])
        after = history([(0.0, 1.2), (10.0, 0.99), (20.0, 0.88)])
        shares = absorption(
            [(10.0, 1.0, {"A": (0.01, 0.5), "B": (0.002, 1.0), "C": (0.5, 1.0)})]
        )
        (row,) = compare_added.compare(before, after, shares, ["A", "B"])["points"]
        self.assertAlmostEqual(row["direct_pcm"], (1 / 1.0 - 1 / 0.99) * 1e5)
        self.assertAlmostEqual(row["decay_only_pcm"], -1e5 * 0.012)
        self.assertAlmostEqual(row["captured_pcm"], -1e5 * 0.007)

    def test_added_nuclides_missing_from_the_absorption_shares_are_errors(self):
        h = history([(10.0, 1.0)])
        with self.assertRaises(SystemExit):
            compare_added.compare(
                h, h, absorption([(10.0, 1.0, {"A": (0.1, 1.0)})]), ["Z"]
            )

    def test_histories_with_other_burnup_points_are_errors(self):
        with self.assertRaises(SystemExit):
            compare_added.compare(
                history([(10.0, 1.0)]),
                history([(20.0, 1.0)]),
                absorption([(10.0, 1.0, {})]),
                [],
            )


if __name__ == "__main__":
    unittest.main()
