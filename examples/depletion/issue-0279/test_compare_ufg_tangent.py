"""compare_ufg_tangent.py の比較（ADR-0247 の段 I2-2。#0279）。"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_ufg_tangent as cut


def report(
    flux: list[list[float]], macroscopic: list[float], implicit: list[list[float]]
) -> dict:
    """1 本のピン、燃料 2 リング、2 群、反応 1 つ（全断面積）の結果。`implicit` は [向き][リング] の、群によらない暗黙の分"""
    return {
        "ufg_tangent": {
            "nuclides": [f"N-{d}" for d in range(len(implicit))],
            "reactions": ["total"],
            "pins": [
                {
                    "first_group": 3,
                    "flux": flux,
                    "macroscopic": [[[m, m]] for m in macroscopic],
                    "implicit": [
                        [[[t, t]] for t in direction] for direction in implicit
                    ],
                }
            ],
        }
    }


FLUX = [[1.0, 1.0], [1.0, 1.0]]


class Compare(unittest.TestCase):
    def test_the_same_tangent_has_no_error(self):
        # 状態が変わっても暗黙の分が同じなら、公称の係数のずれは 0。相対の暗黙の分は Σ φ T / Σ φ Σ
        nominal = report(FLUX, [2.0, 2.0], [[-0.1, -0.3], [0.01, 0.0]])
        out = cut.compare(
            nominal, {"hot": report(FLUX, [4.0, 4.0], [[-0.1, -0.3], [0.01, 0.0]])}
        )
        row = out["nuclides"]["N-0"]["total"]
        self.assertAlmostEqual(row["nominal"]["relative"], -0.1)
        self.assertAlmostEqual(row["hot"]["relative"], -0.05)
        self.assertEqual(row["hot"]["error"], 0.0)
        # 大きさは反応ごとの最大の向きとの比
        self.assertEqual(row["size"], 1.0)
        self.assertAlmostEqual(
            out["nuclides"]["N-1"]["total"]["size"], 0.01 / (0.1**2 + 0.3**2) ** 0.5
        )

    def test_a_scaled_tangent_has_the_relative_error_of_the_scale(self):
        # 状態の暗黙の分が公称の 1.1 倍なら、ずれは 0.1 / 1.1
        nominal = report(FLUX, [2.0, 2.0], [[-0.1, -0.3]])
        state = report(FLUX, [2.0, 2.0], [[-0.11, -0.33]])
        error = cut.compare(nominal, {"s": state})["nuclides"]["N-0"]["total"]["s"][
            "error"
        ]
        self.assertAlmostEqual(error, 0.1 / 1.1)

    def test_states_of_other_rings_stop(self):
        nominal = report(FLUX, [2.0, 2.0], [[-0.1, -0.3]])
        state = report([[1.0, 1.0]], [2.0], [[-0.1]])
        with self.assertRaises(SystemExit):
            cut.compare(nominal, {"s": state})

    def test_states_of_other_directions_stop(self):
        nominal = report(FLUX, [2.0, 2.0], [[-0.1, -0.3]])
        with self.assertRaises(SystemExit):
            cut.compare(
                nominal, {"s": report(FLUX, [2.0, 2.0], [[-0.1, -0.3], [0.0, 0.0]])}
            )


if __name__ == "__main__":
    unittest.main()
