"""compare_dilute.py（輸送しない核種を無限希釈の吸収体として数えた燃焼と、65 核種を輸送する核種にした燃焼の k∞ の差）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_dilute.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_dilute


def history(ks, transportable, dilute=(), shielding=1.0):
    h = {
        "format": "xlc-depletion-history",
        "transportable": list(transportable),
        "points": [{"burnup": b, "k_inf": k} for b, k in ks],
        "seconds": {
            "shielding": shielding,
            "moc": 2.0,
            "rates": 0.5,
            "matrix": 0.1,
            "cram": 0.2,
        },
        "shieldings": 3,
    }
    if dilute:
        h["dilute_absorbers"] = list(dilute)
    return h


# 輸送する核種 T、65 核種にあたる A（自己遮蔽する S を含む）、残り R
BASE = ["T"]
ADDED = ["A", "S"]
SHIELDED = ["S"]


def histories(k):
    """燃焼度 0 と 10 の点の k∞ を、名前ごとに `k[name]` とした 5 つの履歴"""
    return {
        "before": history([(0.0, 1.2), (10.0, k["before"])], BASE),
        "added": history([(0.0, 1.2), (10.0, k["added"])], BASE + ADDED),
        "dilute": history([(0.0, 1.2), (10.0, k["dilute"])], BASE, ADDED + ["R"], 1.1),
        "added_dilute": history(
            [(0.0, 1.2), (10.0, k["added_dilute"])], BASE + ADDED, ["R"]
        ),
        "shielded_dilute": history(
            [(0.0, 1.2), (10.0, k["shielded_dilute"])], BASE + SHIELDED, ["A", "R"], 1.2
        ),
    }


def rho(k_from, k_to):
    return (1.0 / k_from - 1.0 / k_to) * 1e5


class Compare(unittest.TestCase):
    def test_differences_split_the_dilute_error_from_the_rest_absorption(self):
        k = {
            "before": 1.0,
            "added": 0.99,
            "dilute": 0.985,
            "added_dilute": 0.986,
            "shielded_dilute": 0.9855,
        }
        summary = compare_dilute.compare(histories(k), ADDED, SHIELDED)
        row = summary["points"][1]
        self.assertEqual(row["burnup"], 10.0)
        self.assertAlmostEqual(row["dilute_pcm"], rho(0.99, 0.985))
        self.assertAlmostEqual(row["shielded_dilute_pcm"], rho(0.99, 0.9855))
        self.assertAlmostEqual(row["rest_pcm"], rho(0.99, 0.986))
        self.assertAlmostEqual(row["dilute_added_pcm"], rho(0.986, 0.985))
        self.assertAlmostEqual(row["dilute_unshielded_pcm"], rho(0.986, 0.9855))
        # 内訳が閉じる: 全体の差 = 残りの吸収 + 足した核種を無限希釈で数えた差
        self.assertAlmostEqual(
            row["dilute_pcm"], row["rest_pcm"] + row["dilute_added_pcm"]
        )
        self.assertEqual(summary["points"][0]["dilute_pcm"], 0.0)

    def test_shielding_seconds_are_compared_with_the_history_without_absorbers(self):
        k = {
            name: 1.0
            for name in ["before", "added", "dilute", "added_dilute", "shielded_dilute"]
        }
        summary = compare_dilute.compare(histories(k), ADDED, SHIELDED)
        self.assertAlmostEqual(summary["seconds"]["dilute"]["shielding_ratio"], 1.1)
        self.assertAlmostEqual(
            summary["seconds"]["shielded_dilute"]["shielding_ratio"], 1.2
        )
        self.assertEqual(summary["seconds"]["dilute"]["shieldings"], 3)
        self.assertEqual(summary["counts"]["dilute"], 3)

    def test_histories_with_other_nuclide_sets_are_errors(self):
        k = {
            name: 1.0
            for name in ["before", "added", "dilute", "added_dilute", "shielded_dilute"]
        }
        for name, key, value in [
            ("added", "transportable", BASE + ["A"]),
            ("dilute", "dilute_absorbers", ["A", "R"]),
            ("added_dilute", "dilute_absorbers", ["A", "R"]),
            ("shielded_dilute", "transportable", BASE),
            ("shielded_dilute", "dilute_absorbers", ["R"]),
            ("before", "dilute_absorbers", ["R"]),
        ]:
            h = histories(k)
            h[name][key] = value
            with self.assertRaises(SystemExit, msg=f"{name} {key}"):
                compare_dilute.compare(h, ADDED, SHIELDED)

    def test_histories_with_other_burnup_points_are_errors(self):
        k = {
            name: 1.0
            for name in ["before", "added", "dilute", "added_dilute", "shielded_dilute"]
        }
        h = histories(k)
        h["dilute"]["points"][1]["burnup"] = 20.0
        with self.assertRaises(SystemExit):
            compare_dilute.compare(h, ADDED, SHIELDED)


if __name__ == "__main__":
    unittest.main()
