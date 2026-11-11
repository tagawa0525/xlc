"""doppler_exponents.py の J 関数と局所のべき（#0274）。

実行: cd examples/depletion/issue-0274 && python3 -m unittest test_doppler_exponents
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import doppler_exponents as de


class JFunctionTest(unittest.TestCase):
    def test_a_resonance_much_wider_than_the_doppler_width_keeps_the_lorentzian_integral(
        self,
    ):
        # ψ → 1/(1 + x²) のとき ∫ 1/(1 + β(1 + x²)) dx = π / √(β(1 + β))
        for beta in (1e-3, 0.1, 10.0):
            self.assertAlmostEqual(
                de.j_function(1e3, beta) / (math.pi / math.sqrt(beta * (1 + beta))),
                1.0,
                places=5,
            )

    def test_the_doppler_width_shrinks_as_the_square_root_of_the_temperature(self):
        self.assertAlmostEqual(
            de.xi_at(6.67, 0.024, 236.0, 400.0) / de.xi_at(6.67, 0.024, 236.0, 1600.0),
            2.0,
            places=12,
        )


class ExponentTest(unittest.TestCase):
    def test_a_dilute_resonance_far_wider_in_doppler_than_natural_falls_as_t_to_minus_three_halves(
        self,
    ):
        n = de.local_exponent(1e5, 0.024, 236.0, 1e3, 900.0)
        self.assertAlmostEqual(n, 1.5, delta=0.01)

    def test_a_strongly_shielded_resonance_grows_linearly_in_temperature(self):
        n = de.local_exponent(6.67, 0.024, 236.0, 1e-4, 900.0)
        self.assertLess(abs(n), 0.1)

    def test_shielding_between_the_limits_gives_exponents_between_them(self):
        ns = [
            de.local_exponent(6.67, 0.024, 236.0, beta, 900.0)
            for beta in (1e-3, 1e-2, 1e-1, 1.0)
        ]
        self.assertEqual(ns, sorted(ns))
        self.assertTrue(0.0 < ns[1] < 0.5 < ns[2] < 1.5)


if __name__ == "__main__":
    unittest.main()
