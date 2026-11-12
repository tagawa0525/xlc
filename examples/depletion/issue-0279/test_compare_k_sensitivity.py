"""compare_k_sensitivity.py の比較（ADR-0247 の段 2。#0279）。"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_k_sensitivity as cks


def report(
    adjoint: float,
    central: tuple[float, float],
    epsilon: float = 0.1,
    tau: float = 1e-12,
) -> dict:
    coarse, fine = central
    return {
        "epsilon": epsilon,
        "k_tolerance": tau,
        "nuclides": [
            {
                "name": "Pu-240",
                "relative": adjoint,
                "direct": {
                    "central": [coarse, fine],
                    "richardson": fine + (fine - coarse) / 3.0,
                },
            }
        ],
    }


class Compare(unittest.TestCase):
    def test_agreement_within_the_step_error_passes(self):
        (row,) = cks.compare(report(-0.05926612, (-0.05925, -0.059266)))
        self.assertAlmostEqual(row["direct"], -0.059266 + (-0.059266 + 0.05925) / 3.0)
        self.assertEqual(row["limited_by"], "step")
        self.assertTrue(row["ok"])

    def test_difference_beyond_every_limit_fails(self):
        (row,) = cks.compare(report(-0.0600, (-0.059266, -0.059266)))
        self.assertEqual(row["limited_by"], "relative")
        self.assertFalse(row["ok"])

    def test_small_sensitivities_are_limited_by_the_k_tolerance(self):
        # k の変化が判定の誤差に近い核種: 3 τ / ε = 3e-9 の差まで許す
        (row,) = cks.compare(report(-5.0e-9, (-7.0e-9, -7.0e-9), tau=1e-10))
        self.assertEqual(row["limited_by"], "k tolerance")
        self.assertTrue(row["ok"])
        (row,) = cks.compare(report(-5.0e-9, (-9.0e-9, -9.0e-9), tau=1e-10))
        self.assertFalse(row["ok"])


if __name__ == "__main__":
    unittest.main()
