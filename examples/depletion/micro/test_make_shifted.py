"""make_shifted.py（燃焼点 k の組成を、点 k - 1 の組成のうち追う核種だけを点 k の値にしたものにした履歴）。

実行: python3 -m unittest discover -s examples/depletion/micro -p 'test_make_shifted.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_shifted

NUCLIDES = ["O-16", "U-235", "Xe-135", "Cs-133"]


def history() -> dict:
    """燃焼点 3 つ × 領域 2 つ。数密度は (点, 領域, 核種) から決まる"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": NUCLIDES,
        "regions": [
            {"material": "fuel", "cell": 20, "regions": [3], "volume": 2.0},
            {"material": "fuel", "cell": 21, "regions": [4], "volume": 4.0},
        ],
        "transportable": NUCLIDES,
        "points": [
            {
                "burnup": burnup,
                "k_inf": 1.0 + p,
                "densities": [
                    [100 * p + 10 * r + n for n in range(len(NUCLIDES))]
                    for r in range(2)
                ],
            }
            for p, burnup in enumerate([0.0, 10.0, 20.0])
        ],
    }


class Shifted(unittest.TestCase):
    def test_tracked_nuclides_take_the_point_and_the_others_the_point_before(self):
        original = history()
        shifted = make_shifted.shifted(copy.deepcopy(original), ["U-235", "Xe-135"])
        # 燃焼度 0 は入力の組成と比べるので、そのまま
        self.assertEqual(shifted["points"][0], original["points"][0])
        for p in (1, 2):
            for r in range(2):
                self.assertEqual(
                    shifted["points"][p]["densities"][r],
                    [
                        100 * (p - 1) + 10 * r,
                        100 * p + 10 * r + 1,
                        100 * p + 10 * r + 2,
                        100 * (p - 1) + 10 * r + 3,
                    ],
                )
        # 燃焼度とほかの項目はそのまま
        self.assertEqual([p["burnup"] for p in shifted["points"]], [0.0, 10.0, 20.0])
        self.assertEqual(shifted["regions"], original["regions"])
        self.assertEqual(shifted["nuclides"], original["nuclides"])

    def test_tracked_nuclides_must_be_transported_nuclides_of_a_history(self):
        with self.assertRaises(SystemExit):
            make_shifted.shifted(history(), ["Pu-239"])
        not_transported = history()
        not_transported["transportable"] = ["O-16", "U-235"]
        with self.assertRaises(SystemExit):
            make_shifted.shifted(not_transported, ["Xe-135"])
        other = history()
        other["format"] = "xlc-table"
        with self.assertRaises(SystemExit):
            make_shifted.shifted(other, ["U-235"])

    def test_the_tracked_nuclides_are_the_micro_list_of_the_input(self):
        text = (
            '[history]\nfile = "h.msgpack"\nburnup = [0]\nmicro = ["U-235", "Xe-135"]\n'
        )
        self.assertEqual(make_shifted.tracked(text), ["U-235", "Xe-135"])
        with self.assertRaises(SystemExit):
            make_shifted.tracked('[history]\nfile = "h.msgpack"\nburnup = [0]\n')


if __name__ == "__main__":
    unittest.main()
