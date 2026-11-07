"""compare_core_energy.py の要約（ノードの出力の割合を核分裂率と κΣf で求めた core の燃焼の差）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_energy.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_energy


def core(powers, burnups, k=1.0, boron=600.0):
    """集合体 3 体（3 体目は反射体で、核分裂がない）・軸方向のノード 2 つの燃焼の履歴。点は 0 と 10 EFPD で、`powers` と `burnups` は
    点ごとの [集合体][ノード]"""
    return {
        "specific_power": 40.0,
        "states": {"aro": {}},
        "axial_mesh": [0.0, 10.0, 20.0],
        "assemblies": ["A", "B", "R"],
        "heavy_metal": [[1.0, 1.0], [1.0, 1.0], [0.0, 0.0]],
        "points": [
            {
                "efpd": efpd,
                "days": efpd,
                "state": "aro",
                "k": k,
                "boron": boron,
                "solves": 1,
                "burnups": b,
                "powers": p,
            }
            for efpd, p, b in zip((0.0, 10.0), powers, burnups)
        ],
    }


RATE_POWERS = [
    [[0.25, 0.25], [0.25, 0.25], [0.0, 0.0]],
    [[0.2, 0.2], [0.3, 0.3], [0.0, 0.0]],
]
BURNUPS = [
    [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
    [[0.4, 0.4], [0.6, 0.6], [0.0, 0.0]],
]


class Compare(unittest.TestCase):
    def test_differences_are_the_energy_minus_the_fission_rate_shares(self):
        rate = core(RATE_POWERS, BURNUPS)
        energy_powers = [
            [[0.25, 0.25], [0.25, 0.25], [0.0, 0.0]],
            [[0.19, 0.2], [0.31, 0.3], [0.0, 0.0]],
        ]
        energy_burnups = [
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            [[0.4, 0.38], [0.6, 0.62], [0.0, 0.0]],
        ]
        energy = core(energy_powers, energy_burnups, k=1.002, boron=610.0)
        c = compare_core_energy.compare(rate, energy)
        self.assertEqual(c["assemblies"], ["A", "B", "R"])
        first, later = c["points"]
        self.assertEqual(
            (later["efpd"], later["days"], later["state"]), (10.0, 10.0, "aro")
        )
        self.assertEqual(first["power"], 0.0)
        self.assertEqual(later["boron"], [600.0, 610.0])
        self.assertEqual(later["k"], [1.0, 1.002])
        # 反応度の差 (1/k_核分裂率 − 1/k_κΣf) × 1e5 pcm
        self.assertAlmostEqual(later["reactivity"], (1.0 - 1.0 / 1.002) * 1e5, places=9)
        # ノードの出力の割合の相対差の絶対値の最大（核分裂のないノードを除く）は A の下の 0.19 / 0.2 − 1
        self.assertAlmostEqual(later["power"], 0.05, places=12)
        self.assertAlmostEqual(later["burnup"], 0.05, places=12)
        # 集合体ごとの出力の割合（ノードの和）と、質量の重みの燃焼度
        self.assertEqual(len(later["assembly_powers"]["rate"]), 3)
        self.assertAlmostEqual(later["assembly_powers"]["rate"][0], 0.4, places=12)
        self.assertAlmostEqual(later["assembly_powers"]["energy"][0], 0.39, places=12)
        self.assertAlmostEqual(later["assembly_burnups"]["energy"][1], 0.61, places=12)
        self.assertEqual(later["assembly_powers"]["energy"][2], 0.0)

    def test_histories_must_share_the_points_and_the_core(self):
        rate = core(RATE_POWERS, BURNUPS)
        other = core(RATE_POWERS, BURNUPS)
        other["points"][1]["efpd"] = 11.0
        with self.assertRaises(SystemExit):
            compare_core_energy.compare(rate, other)
        other = core(RATE_POWERS, BURNUPS)
        other["heavy_metal"] = [[1.0, 1.0], [2.0, 1.0], [0.0, 0.0]]
        with self.assertRaises(SystemExit):
            compare_core_energy.compare(rate, other)

    def test_report_gives_the_largest_differences_per_case(self):
        rate = core(RATE_POWERS, BURNUPS)
        energy = core(
            [RATE_POWERS[0], [[0.19, 0.2], [0.31, 0.3], [0.0, 0.0]]],
            BURNUPS,
            boron=610.0,
        )
        text = "\n".join(
            compare_core_energy.report(
                {"cases": {"reload": compare_core_energy.compare(rate, energy)}}
            )
        )
        self.assertIn("reload: ホウ素の差の絶対値の最大 10.0000 ppm", text)
        self.assertIn("ノードの出力の割合の相対差の絶対値の最大 5.00e-02", text)
        self.assertIn("集合体の出力の割合の相対差の絶対値の最大 2.50e-02", text)


if __name__ == "__main__":
    unittest.main()
