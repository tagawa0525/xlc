"""fission_energy.py（#0293 の項目 3 の、核分裂あたりのエネルギーの定義の比べ）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_fission_energy.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fission_energy


class FissionRowTest(unittest.TestCase):
    def test_compares_er_and_the_local_fission_heating_with_qi(self):
        row = fission_energy.fission_row(qi=200.0e6, er=202.0e6, local=190.0e6)
        self.assertAlmostEqual(row["er_over_qi"], 0.01)
        self.assertAlmostEqual(row["local_over_qi"], -0.05)


class KappaRowTest(unittest.TestCase):
    def test_converts_the_kappa_joules_to_ev_before_comparing(self):
        joule = fission_energy.EV_IN_JOULES
        row = fission_energy.kappa_row(
            qi=200.0e6, et=210.0e6, kappa=212.1e6 * joule, kappa0=200.02e6 * joule
        )
        self.assertAlmostEqual(row["kappa0_over_qi"], 1.0e-4)
        self.assertAlmostEqual(row["kappa_over_et"], 0.01)
        self.assertAlmostEqual(row["kappa_over_kappa0"], 212.1 / 200.02)


if __name__ == "__main__":
    unittest.main()
