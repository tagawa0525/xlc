"""compare_depletion_versions.py の要約（版ごとのピンセルの燃焼の、基準の版との k∞、反応度、数密度の差）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_depletion_versions.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_depletion_versions as cdv

NUCLIDES = cdv.NUCLIDES


def density(name, burnup, scale):
    """核種ごとに違う数密度。燃焼度 0 は版によらず同じ（Xe-135 は 0）で、燃焼した点は `scale` 倍"""
    if burnup == 0.0:
        return 0.0 if name == "Xe-135" else 1e-3 * (ALL.index(name) + 1)
    return 1e-3 * (ALL.index(name) + 1) * scale


ALL = NUCLIDES + ["O-16"]


def history(k, burnup=(0.0, 17.0, 60.0), scale=1.0, nuclides=None, volume=0.5):
    """燃焼点ごとの k∞ と、燃焼した点の数密度が `scale` 倍の、1 領域の燃焼の履歴（lattice の `[depletion]` の出力の形）"""
    names = list(nuclides or ALL)
    return {
        "format": "xlc-depletion-history",
        "nuclides": names,
        "regions": [{"material": "fuel", "volume": volume}],
        "transportable": names,
        "points": [
            {
                "burnup": b,
                "k_inf": kk,
                "densities": [[density(n, b, scale) for n in names]],
            }
            for b, kk in zip(burnup, k)
        ],
    }


def four(**changes):
    """VII.0、VII.1、VIII.0、VIII.1 の履歴。`changes` は版ごとの history の引数"""
    base = {"k": (1.2, 1.1, 0.9)}
    return {
        v: history(**{**base, **changes.get(v.replace(".", ""), {})})
        for v in cdv.VERSIONS
    }


class Summarize(unittest.TestCase):
    def test_reactivity_differences_are_against_the_base_version(self):
        s = cdv.summarize(four(VIII0={"k": (1.2, 1.0, 0.8)}), "VII.1")
        diff = s["versions"]["VIII.0"]["reactivity_pcm"]
        expected = [
            (1 / k1 - 1 / k) * 1e5 for k1, k in zip((1.2, 1.1, 0.9), (1.2, 1.0, 0.8))
        ]
        for got, want in zip(diff, expected):
            self.assertAlmostEqual(got, want)
        self.assertLess(diff[1], 0.0)

    def test_the_base_version_has_no_differences(self):
        s = cdv.summarize(four(), "VII.1")
        self.assertEqual(s["base"], "VII.1")
        self.assertNotIn("reactivity_pcm", s["versions"]["VII.1"])
        self.assertEqual(s["versions"]["VII.1"]["k_inf"], [1.2, 1.1, 0.9])

    def test_k_differences_are_in_pcm_of_k(self):
        s = cdv.summarize(four(VII0={"k": (1.201, 1.1, 0.9)}), "VII.1")
        self.assertAlmostEqual(s["versions"]["VII.0"]["k_pcm"][0], 100.0)

    def test_densities_are_relative_to_the_base_version(self):
        s = cdv.summarize(four(VIII1={"scale": 1.02}), "VII.1")
        pin = s["versions"]["VIII.1"]["densities"]
        self.assertEqual(sorted(pin), sorted(NUCLIDES))
        self.assertAlmostEqual(pin["U-235"][2], 0.02)

    def test_a_density_that_is_zero_in_the_base_is_not_compared(self):
        s = cdv.summarize(four(), "VII.1")
        self.assertIsNone(s["versions"]["VII.0"]["densities"]["Xe-135"][0])
        self.assertAlmostEqual(s["versions"]["VII.0"]["densities"]["Xe-135"][1], 0.0)

    def test_densities_are_found_by_name_in_chains_of_other_orders(self):
        histories = four()
        histories["VII.0"] = history(
            (1.2, 1.1, 0.9), nuclides=list(reversed(ALL)), scale=1.03
        )
        s = cdv.summarize(histories, "VII.1")
        for n in NUCLIDES:
            self.assertAlmostEqual(s["versions"]["VII.0"]["densities"][n][2], 0.03)

    def test_a_nuclide_missing_from_a_chain_is_an_error(self):
        histories = four()
        histories["VII.0"] = history(
            (1.2, 1.1, 0.9), nuclides=[n for n in ALL if n != "Sm-149"]
        )
        with self.assertRaises(SystemExit):
            cdv.summarize(histories, "VII.1")

    def test_marks_pick_the_problem9_and_last_points(self):
        s = cdv.summarize(four(VIII0={"k": (1.2, 1.0, 0.8)}), "VII.1")
        marks = {m["burnup"]: m for m in s["marks"]}
        self.assertEqual(sorted(marks), [17.0, 60.0])
        self.assertAlmostEqual(
            marks[17.0]["reactivity_pcm"]["VIII.0"],
            s["versions"]["VIII.0"]["reactivity_pcm"][1],
        )

    def test_largest_differences_up_to_problem9_and_overall(self):
        s = cdv.summarize(four(VIII0={"k": (1.2, 1.0, 0.89)}), "VII.1")
        largest = s["versions"]["VIII.0"]["largest_reactivity_pcm"]
        self.assertEqual(largest["problem9"]["burnup"], 17.0)
        self.assertEqual(largest["all"]["burnup"], 17.0)
        self.assertAlmostEqual(
            largest["all"]["pcm"], s["versions"]["VIII.0"]["reactivity_pcm"][1]
        )

    def test_histories_must_have_the_same_burnup_points(self):
        h = four(VII0={"burnup": (0.0, 17.0, 59.0)})
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")

    def test_histories_must_reach_the_problem9_burnup(self):
        h = {
            v: history((1.2, 1.1, 0.9), burnup=(0.0, 10.0, 60.0)) for v in cdv.VERSIONS
        }
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")

    def test_histories_must_have_the_same_regions(self):
        h = four(VIII1={"volume": 0.6})
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")

    def test_histories_must_have_the_same_initial_composition(self):
        h = four()
        h["VII.0"] = copy.deepcopy(h["VII.0"])
        h["VII.0"]["points"][0]["densities"][0][0] *= 1.001
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")

    def test_all_four_versions_are_required(self):
        h = four()
        del h["VIII.1"]
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")

    def test_a_base_outside_the_versions_is_an_error(self):
        with self.assertRaises(SystemExit):
            cdv.summarize(four(), "VI.8")

    def test_another_format_is_an_error(self):
        h = four()
        h["VII.0"] = {**h["VII.0"], "format": "openmc-depletion"}
        with self.assertRaises(SystemExit):
            cdv.summarize(h, "VII.1")


class Arguments(unittest.TestCase):
    def test_labelled_paths_are_split_at_the_first_equals_sign(self):
        self.assertEqual(
            cdv.labelled(["VII.0=endfb70/a.msgpack", "VII.1=endfb71/b=c.msgpack"]),
            {"VII.0": Path("endfb70/a.msgpack"), "VII.1": Path("endfb71/b=c.msgpack")},
        )

    def test_a_path_without_a_label_is_an_error(self):
        with self.assertRaises(SystemExit):
            cdv.labelled(["endfb70/a.msgpack"])

    def test_a_label_given_twice_is_an_error(self):
        with self.assertRaises(SystemExit):
            cdv.labelled(["VII.0=a.msgpack", "VII.0=b.msgpack"])


class Report(unittest.TestCase):
    def test_report_shows_the_marks_and_the_base(self):
        s = {
            "format": "xlc-compare-depletion-versions",
            "inputs": [{"path": "endfb70/x.msgpack", "md5": "0" * 32}],
            **cdv.summarize(four(VIII0={"k": (1.2, 1.0, 0.8)}), "VII.1"),
        }
        text = "\n".join(cdv.report(s))
        self.assertIn("基準 VII.1", text)
        self.assertIn("17 MWd/kg", text)
        self.assertIn("60 MWd/kg", text)

    def test_report_needs_the_inputs(self):
        s = cdv.summarize(four(), "VII.1")
        with self.assertRaises(SystemExit):
            cdv.report(s)


if __name__ == "__main__":
    unittest.main()
