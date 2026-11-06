"""compare_core_acceleration.py の要約（固有値計算の加速だけを変えた core の燃焼の差と、固有値計算ごとの反復の回数と時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_acceleration.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_acceleration
from test_compare_core_guess import core, log


class Compare(unittest.TestCase):
    def test_differences_are_the_cmfd_minus_the_unaccelerated_solves(self):
        plain = core()
        cmfd = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4)
        c = compare_core_acceleration.compare(
            plain,
            log([4000, 1200, 900, 3000], midpoints=((1.0, 900.0, 1),), total=9.0),
            cmfd,
            log([60, 20, 15, 40], midpoints=((1.0, 900.0 + 2e-4, 1),), total=1.5),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(sorted(c["runs"]), ["cmfd", "none"])
        self.assertEqual(
            c["runs"]["cmfd"],
            {
                "solves": 4,
                "iterations": 135,
                "eigenvalue_seconds": 1.0,
                "seconds": 1.5,
                "per_solve": {"iterations": [15, 60], "seconds": [0.25, 0.25]},
            },
        )
        self.assertEqual(c["runs"]["none"]["per_solve"]["iterations"], [900, 4000])

    def test_report_gives_the_ratios_and_the_range_per_solve(self):
        c = compare_core_acceleration.compare(
            core(),
            log([4000, 1200, 900, 3000], total=9.0),
            core(),
            log([60, 20, 15, 40], total=1.5),
        )
        text = "\n".join(compare_core_acceleration.report({"cases": {"small": c}}))
        self.assertIn("外側の反復 9,100 → 135 回（0.015 倍）", text)
        self.assertIn("1 回の外側の反復 900〜4,000 → 15〜60 回", text)


if __name__ == "__main__":
    unittest.main()
