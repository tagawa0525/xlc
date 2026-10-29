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


def mean(history, p, i):
    """点 p の核種 i の、燃焼の領域の体積で重みを付けた数密度（負は 0）"""
    regions = history["regions"]
    densities = history["points"][p]["densities"]
    total = sum(d["volume"] for d in regions)
    return sum(d["volume"] * max(n[i], 0.0) for d, n in zip(regions, densities)) / total


class Spread(unittest.TestCase):
    def test_tracked_nuclides_keep_the_mean_of_the_point_and_the_spread_of_the_point_before(
        self,
    ):
        original = history()
        spread = make_shifted.spread(copy.deepcopy(original), ["U-235", "Xe-135"])
        self.assertEqual(spread["points"][0], original["points"][0])
        for p in (1, 2):
            before = original["points"][p - 1]["densities"]
            got = spread["points"][p]["densities"]
            for i in (1, 2):
                self.assertAlmostEqual(
                    mean(spread, p, i),
                    mean(original, p, i),
                    delta=1e-12 * mean(original, p, i),
                )
                # 領域の間の比は前の点のもの
                self.assertAlmostEqual(
                    got[1][i] / got[0][i], before[1][i] / before[0][i], delta=1e-12
                )
            # 追わない核種は前の点のまま
            for r in range(2):
                self.assertEqual(got[r][0], before[r][0])
                self.assertEqual(got[r][3], before[r][3])

    def test_tracked_nuclides_absent_at_the_point_before_spread_uniformly(self):
        # 前の点で 0 の核種（燃焼度 0 の Xe-135 など）は、燃焼の領域で一様（微視的断面積の重み 1。ADR-0216 の案 Z1）
        original = history()
        for r in range(2):
            original["points"][0]["densities"][r][2] = 0.0
        spread = make_shifted.spread(copy.deepcopy(original), ["Xe-135"])
        got = spread["points"][1]["densities"]
        self.assertEqual(got[0][2], got[1][2])
        self.assertAlmostEqual(got[0][2], mean(original, 1, 2), delta=1e-12 * got[0][2])

    def test_negative_densities_count_as_zero(self):
        # lattice は負の数密度を 0 として輸送し、均質化する
        original = history()
        original["points"][0]["densities"][0][1] = -1e-30
        spread = make_shifted.spread(copy.deepcopy(original), ["U-235"])
        got = spread["points"][1]["densities"]
        self.assertEqual(got[0][1], 0.0)
        self.assertAlmostEqual(mean(spread, 1, 1), mean(original, 1, 1), delta=1e-9)


class Previous(unittest.TestCase):
    def test_each_point_takes_the_densities_of_the_point_before(self):
        original = history()
        previous = make_shifted.previous(copy.deepcopy(original))
        self.assertEqual(previous["points"][0], original["points"][0])
        for p in (1, 2):
            self.assertEqual(
                previous["points"][p]["densities"],
                original["points"][p - 1]["densities"],
            )
            self.assertEqual(
                previous["points"][p]["burnup"], original["points"][p]["burnup"]
            )
        other = history()
        other["format"] = "xlc-table"
        with self.assertRaises(SystemExit):
            make_shifted.previous(other)


if __name__ == "__main__":
    unittest.main()
