"""compare.py の反応度の変化（#0180、ADR-0249）。

実行: cd examples/depletion/issue-0180/fuels && python3 -m unittest test_compare
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import reactivity, same_nominal, share


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


class SameNominal(unittest.TestCase):
    def test_tables_that_solve_the_same_states_have_the_same_nominal(self):
        # 微分だけが違うテーブル
        self.assertTrue(same_nominal(table(10.0, 1e-4), table(20.0, 1e-4)))

    def test_a_different_state_is_not_the_same_nominal(self):
        # 摂動の自己遮蔽のほかで状態が違えば（別の入力、別の条件）、差に状態の差が混ざる
        other = table(10.0, 1e-4)
        other["states"][0]["k_inf"] = 1.00001
        self.assertFalse(same_nominal(table(10.0, 1e-4), other))
        other = table(10.0, 1e-4)
        other["states"][0]["total"][1] = 1.0000001
        self.assertFalse(same_nominal(table(10.0, 1e-4), other))
        # 2 群の定数のほかの項目（χ、拡散係数、不連続因子など）も状態の一部
        other = table(10.0, 1e-4)
        other["states"][0]["chi"] = [0.999, 0.001]
        self.assertFalse(same_nominal(table(10.0, 1e-4), other))
        # 追う核種の均質化した数密度は微分に掛けるので、これも同じであること
        other = table(10.0, 1e-4)
        other["micro"]["densities"][0][1] = 2e-4
        self.assertFalse(same_nominal(table(10.0, 1e-4), other))


class Share(unittest.TestCase):
    def test_the_share_is_the_difference_over_the_solved_value(self):
        self.assertAlmostEqual(share(-1.0, -200.0), 0.005, places=15)

    def test_a_zero_solved_value_has_no_share(self):
        # Pu-240 の数密度が 0 の燃焼点（燃焼度 0）は反応度の変化が 0 で、割合を持たない
        self.assertIsNone(share(0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
