"""compare_branch_shielding.py の要約（分岐の自己遮蔽を公称の減速材で解いたテーブルの、状態ごとに解き直したテーブルとの差。
ADR-0201）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_branch_shielding.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_branch_shielding as cbs


def state(k: float, scale: float = 1.0) -> dict:
    return {
        "k_inf": k,
        "total": [0.5, 1.0],
        "absorption": [0.01 * scale, 0.1],
        "nu_fission": [0.005, 0.15 * scale],
        "fission": [0.002, 0.06],
        "chi": [1.0, 0.0],
        "diffusion": [1.4, 0.4 * scale],
        "scatter": [[0.47, 0.02 * scale], [0.0, 0.9]],
        "adf": [[1.0, 1.0 * scale]] + [[1.0, 1.0]] * 5,
        "vacuum_ratio": [None] * 6,
        "normalized_adf": [False] * 6,
        "height": None,
        "pins": None,
    }


def table(states: list[dict]) -> dict:
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "axes": [
            {"name": "burnup", "points": [0.0, 20.0]},
            {"name": "boron", "points": [900.0, 1300.0], "b10": 0.199},
            {"name": "moderator", "points": [565.0], "density": [0.743]},
            {"name": "fuel", "points": [565.0]},
        ],
        "split": 1,
        "heavy_metal": 1.0,
        "states": states,
    }


class CompareTest(unittest.TestCase):
    def test_states_carry_their_axis_values_and_the_k_difference_in_pcm(self):
        every = table([state(1.1), state(1.2), state(1.0), state(0.9)])
        nominal = table([state(1.1), state(1.20012), state(1.0), state(0.89991)])
        case = cbs.compare("2A", every, nominal)
        self.assertEqual(case["name"], "2A")
        self.assertEqual(len(case["states"]), 4)
        s = case["states"][1]
        self.assertEqual(
            (
                s["burnup"],
                s["boron"],
                s["moderator_temperature"],
                s["density"],
                s["fuel"],
            ),
            (0.0, 1300.0, 565.0, 0.743, 565.0),
        )
        self.assertAlmostEqual(s["dk_pcm"], 12.0, places=6)
        self.assertAlmostEqual(case["states"][3]["dk_pcm"], -9.0, places=6)
        self.assertEqual(case["states"][0]["dk_pcm"], 0.0)

    def test_constants_report_the_largest_relative_difference_per_item(self):
        every = table([state(1.0)] * 4)
        nominal = table([state(1.0), state(1.0, 1.001), state(1.0), state(1.0, 0.998)])
        case = cbs.compare("2A", every, nominal)
        worst = case["states"][3]["relative"]
        for item in ["absorption", "nu_fission", "diffusion", "scatter", "adf"]:
            self.assertAlmostEqual(worst[item], 0.002, places=12, msg=item)
        # 相対差の分母は 2 つの値の絶対値の大きい方
        self.assertAlmostEqual(
            case["states"][1]["relative"]["absorption"], 0.001 / 1.001, places=12
        )
        self.assertEqual(case["states"][0]["relative"]["adf"], 0.0)

    def test_tables_on_different_axes_are_errors(self):
        every = table([state(1.0)] * 4)
        other = table([state(1.0)] * 4)
        other["axes"][1]["points"] = [900.0, 1700.0]
        with self.assertRaises(SystemExit):
            cbs.compare("2A", every, other)
        short = table([state(1.0)] * 3)
        with self.assertRaises(SystemExit):
            cbs.compare("2A", every, short)

    def test_split_tables_are_errors(self):
        # 部分ノードに分けたテーブルは、分岐の状態ごとに split² の状態を持つ
        every = table([state(1.0)] * 16)
        every["split"] = 2
        with self.assertRaises(SystemExit):
            cbs.compare("2A", every, every)

    def test_values_that_are_not_finite_are_errors(self):
        every = table([state(1.0)] * 4)
        nominal = table([state(1.0)] * 3 + [state(1.0, float("nan"))])
        with self.assertRaises(SystemExit):
            cbs.compare("2A", every, nominal)
        nominal = table([state(1.0)] * 3 + [state(float("nan"))])
        with self.assertRaises(SystemExit):
            cbs.compare("2A", every, nominal)


if __name__ == "__main__":
    unittest.main()
