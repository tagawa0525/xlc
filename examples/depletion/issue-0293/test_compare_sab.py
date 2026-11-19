"""compare_sab.py（#0293 の項目 8 の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_compare_sab.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_sab


def table(burnup_indices, ks):
    return {
        "axes": [
            {"name": "burnup", "points": [0.0, 20.0]},
            {"name": "boron", "points": [1300.0]},
        ],
        "points": [[b, 0] for b in burnup_indices],
        "states": [{"k_inf": k} for k in ks],
    }


class Compare(unittest.TestCase):
    def test_the_difference_is_interpolated_minus_evaluated_per_burnup(self):
        rows = compare_sab.compare(
            table([0, 1], [1.2, 1.0]), table([0, 1], [1.2, 1.00002])
        )
        self.assertEqual([r["burnup"] for r in rows], [0.0, 20.0])
        self.assertAlmostEqual(rows[1]["dk_pcm"], 2.0)

    def test_tables_of_different_points_are_an_error(self):
        with self.assertRaises(SystemExit):
            compare_sab.compare(table([0, 1], [1.2, 1.0]), table([0], [1.2]))


if __name__ == "__main__":
    unittest.main()
