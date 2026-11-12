"""compare.py の反応度の変化（#0180、ADR-0249）。

実行: cd examples/depletion/issue-0180/fuels && python3 -m unittest test_compare
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import reactivity


def table(nu_fission_derivative: float, density: float) -> dict:
    """2 群の 1 状態のテーブル。Pu-240 の微分は νΣf の高速群だけ"""
    zero = [0.0, 0.0]
    return {
        "states": [
            {
                "total": [0.5, 1.0],
                "nu_fission": [0.01, 0.2],
                "scatter": [[0.47, 0.02], [0.0, 0.9]],
                "chi": [1.0, 0.0],
                "k_inf": 1.0,
            }
        ],
        "micro": {
            "nuclides": ["U-235", "Pu-240"],
            "densities": [[1e-3, density]],
            "states": [
                [
                    {},
                    {
                        "derivative": {
                            "total": zero,
                            "nu_fission": [nu_fission_derivative, 0.0],
                            "scatter": [zero, zero],
                        }
                    },
                ]
            ],
        },
    }


class Reactivity(unittest.TestCase):
    def test_a_zero_derivative_changes_nothing(self):
        self.assertEqual(reactivity(table(0.0, 1e-4), 0)["drho_pcm"], 0.0)

    def test_the_change_is_the_derivative_times_a_tenth_of_the_density(self):
        # νΣf を増やせば k∞ は増え、反応度の変化は正。変化は微分 × 0.1 N̄ で、微分を倍にしても N̄ を倍にしても同じ
        a = reactivity(table(10.0, 1e-4), 0)["drho_pcm"]
        self.assertGreater(a, 0.0)
        self.assertAlmostEqual(
            reactivity(table(20.0, 0.5e-4), 0)["drho_pcm"], a, places=9
        )
        self.assertAlmostEqual(reactivity(table(5.0, 2e-4), 0)["drho_pcm"], a, places=9)


if __name__ == "__main__":
    unittest.main()
