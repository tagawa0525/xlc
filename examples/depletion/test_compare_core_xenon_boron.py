"""compare_core_xenon_boron.py の要約（平衡 Xe と臨界ホウ素の収束のさせ方だけを変えた core の燃焼の差と、固有値計算の数、反復の回数、時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_xenon_boron.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_xenon_boron
from test_compare_core_guess import core, log


class Compare(unittest.TestCase):
    def test_differences_are_the_joint_minus_the_nested_iteration(self):
        nested = core(solves=(1, 6))
        joint = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4, solves=(1, 3))
        c = compare_core_xenon_boron.compare(
            nested,
            log([100] * 8, midpoints=((1.0, 900.0, 1),), total=9.0),
            joint,
            log([100, 40, 30, 20, 10], midpoints=((1.0, 900.0 + 2e-4, 1),), total=4.0),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(sorted(c["runs"]), ["joint", "nested"])
        self.assertEqual(c["runs"]["nested"]["solves"], 8)
        self.assertEqual(
            c["runs"]["joint"],
            {
                "solves": 5,
                "iterations": 200,
                "eigenvalue_seconds": 1.25,
                "seconds": 4.0,
            },
        )

    def test_report_gives_the_solves_and_their_ratio(self):
        c = compare_core_xenon_boron.compare(
            core(solves=(1, 6)),
            log([100] * 8, total=9.0),
            core(solves=(1, 3)),
            log([50] * 5, total=4.0),
        )
        text = "\n".join(compare_core_xenon_boron.report({"cases": {"small": c}}))
        self.assertIn("固有値計算 8 → 5 回（0.625 倍）", text)


if __name__ == "__main__":
    unittest.main()
