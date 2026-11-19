"""compare_fissionable.py（核分裂の生成が効く核種も輸送する核種にして残りを数えた燃焼と、65 核種を輸送して残りを数えた燃焼の k∞ の差）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_fissionable.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_fissionable


def history(k, transportable, dilute=(), shielding=1.0):
    h = {
        "format": "xlc-depletion-history",
        "transportable": list(transportable),
        "points": [{"burnup": 0.0, "k_inf": 1.2}, {"burnup": 10.0, "k_inf": k}],
        "seconds": {"shielding": shielding},
        "shieldings": 3,
    }
    if dilute:
        h["dilute_absorbers"] = list(dilute)
    return h


# 輸送する核種 T、65 核種にあたる A（自己遮蔽する S と核分裂の F1 を含む）、A の外の核分裂の核種 F2、残り R
BASE = ["T"]
ADDED = ["A", "S", "F1"]
OUTSIDE = ["F2"]
SHIELDED = ["S"]
FISSIONABLE = ["F1", "F2"]


def histories(k):
    return {
        "before": history(k["before"], BASE),
        "added_dilute": history(k["added_dilute"], BASE + ADDED, ["F2", "R"], 2.0),
        "reference": history(k["reference"], BASE + ADDED + OUTSIDE, ["R"], 2.1),
        "form": history(k["form"], BASE + SHIELDED + FISSIONABLE, ["A", "R"], 1.1),
    }


def rho(k_from, k_to):
    return (1.0 / k_from - 1.0 / k_to) * 1e5


def compare(h):
    return compare_fissionable.compare(h, ADDED, OUTSIDE, SHIELDED, FISSIONABLE)


K = {"before": 1.0, "added_dilute": 0.98, "reference": 0.981, "form": 0.9805}


class Compare(unittest.TestCase):
    def test_the_form_is_compared_with_the_reference_that_also_transports_the_outside_fissions(
        self,
    ):
        summary = compare(histories(K))
        row = summary["points"][1]
        self.assertEqual(row["burnup"], 10.0)
        self.assertAlmostEqual(row["form_pcm"], rho(0.981, 0.9805))
        self.assertAlmostEqual(row["rest_fission_pcm"], rho(0.98, 0.981))
        self.assertAlmostEqual(row["form_from_65_pcm"], rho(0.98, 0.9805))
        self.assertAlmostEqual(
            row["form_from_65_pcm"], row["form_pcm"] + row["rest_fission_pcm"]
        )
        self.assertEqual(summary["points"][0]["form_pcm"], 0.0)

    def test_without_the_reference_only_the_difference_from_the_65_nuclides_is_given(
        self,
    ):
        h = histories(K)
        del h["reference"]
        summary = compare(h)
        row = summary["points"][1]
        self.assertAlmostEqual(row["form_from_65_pcm"], rho(0.98, 0.9805))
        self.assertNotIn("form_pcm", row)
        self.assertNotIn("reference", summary["seconds"])

    def test_shielding_seconds_are_compared_with_the_history_without_absorbers(self):
        summary = compare(histories(K))
        self.assertAlmostEqual(summary["seconds"]["form"]["shielding_ratio"], 1.1)
        self.assertAlmostEqual(summary["seconds"]["reference"]["shielding_ratio"], 2.1)
        self.assertEqual(
            summary["counts"],
            {"before": 0, "added_dilute": 2, "reference": 1, "form": 2},
        )

    def test_histories_with_other_nuclide_sets_are_errors(self):
        for name, key, value in [
            ("reference", "transportable", BASE + ADDED),
            ("reference", "dilute_absorbers", ["F2", "R"]),
            ("form", "transportable", BASE + SHIELDED + ["F1"]),
            ("form", "dilute_absorbers", ["R"]),
            ("added_dilute", "dilute_absorbers", ["R"]),
            ("before", "dilute_absorbers", ["R"]),
        ]:
            h = histories(K)
            h[name][key] = value
            with self.assertRaises(SystemExit, msg=f"{name} {key}"):
                compare(h)

    def test_the_outside_list_must_be_the_fissionable_nuclides_outside_the_added_ones(
        self,
    ):
        with self.assertRaises(SystemExit):
            compare_fissionable.compare(histories(K), ADDED, [], SHIELDED, FISSIONABLE)

    def test_histories_with_other_burnup_points_are_errors(self):
        h = histories(K)
        h["form"]["points"][1]["burnup"] = 20.0
        with self.assertRaises(SystemExit):
            compare(h)


if __name__ == "__main__":
    unittest.main()
