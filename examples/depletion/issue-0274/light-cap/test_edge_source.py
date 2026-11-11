"""edge_source.py の散乱源（自由ガスの核を入射エネルギーで静止した標的に切り替えたときの形。ADR-0246）。

実行: cd examples/depletion/issue-0274/light-cap && python3 -m unittest test_edge_source
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import edge_source as es

KT_FUEL = es.BOLTZMANN * 900.0
ALPHA = ((es.O16 - 1) / (es.O16 + 1)) ** 2


def ratio(e: float, cap: float) -> float:
    return (
        es.source(e, cap, es.O16, KT_FUEL) / es.source(e, np.inf, es.O16, KT_FUEL) - 1
    )


class EdgeSourceTest(unittest.TestCase):
    def test_the_free_gas_kernel_keeps_the_static_mean_in_the_bulk(self):
        # E ≫ kT の 1 回の下方散乱の出射の平均は、静止した標的の (1 + α)/2 E に近い
        e = 60.0
        grid = np.linspace(ALPHA * e * 0.95, e, 4001)
        values = np.array([es.kernel(e, x, es.O16, KT_FUEL) for x in grid])
        total = np.trapezoid(values, grid)
        mean = np.trapezoid(values * grid, grid) / total
        self.assertAlmostEqual(mean / e, (1 + ALPHA) / 2, places=3)

    def test_the_cap_digs_the_source_at_alpha_times_the_cap(self):
        cap = 27.14
        self.assertLess(ratio(ALPHA * cap, cap), -0.04)
        # 下端の窓の外では、上限によらない
        self.assertLess(abs(ratio(ALPHA * cap - 2.5, cap)), 1e-3)

    def test_the_cap_raises_the_source_between_alpha_times_the_cap_and_the_cap(self):
        # 上方散乱の側（自由ガス）の得を、下端の側（静止した標的）の損が打ち消さない
        cap = 27.14
        self.assertGreater(ratio(ALPHA * cap + 2.0, cap), 0.01)

    def test_free_hydrogen_lands_below_a_thermal_cut_by_half_a_kt_less(self):
        # 0.625 eV より下に落ちる割合は、静止した標的の 0.625/E から kT/2 だけ減る（E ≫ 0.625 eV ≫ kT）
        kt = es.BOLTZMANN * 580.0
        for e in [20.0, 200.0]:
            f = es.thermal_fraction(e, 0.625, es.H1, kt)
            self.assertAlmostEqual(f * e, 0.625 - kt / 2, delta=2e-4)


if __name__ == "__main__":
    unittest.main()
