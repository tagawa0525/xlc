"""compare_history.py の要約（同じ入力を既定と厳密な経路で解いた 2 つの燃焼の履歴の、燃焼点ごとの反応度の差、物質の種類ごとの
核種の数密度の相対の差、ピン出力の相対の差）。

実行: python3 -m unittest test_compare_history（このディレクトリで）
"""

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_history as compare

NAMES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149", "B-10", "Gd-155", "Gd-157"]


def history(k: list[float], fuel: list[list[float]], ifba: list[float], pins):
    """2 点（燃焼度 0、10）の履歴。燃料 2 領域（体積 1 と 3）と IFBA 1 領域。`fuel` は領域ごとの 2 点の U-235、`ifba` は 2 点の B-10。
    Pu-239 は 0 の点で 0、10 の点で 0.1。"""
    regions = [
        {"material": "fuel", "cell": 6, "volume": 1.0},
        {"material": "fuel", "cell": 7, "volume": 3.0},
        {"material": "ifba", "cell": 6, "volume": 0.5},
    ]
    points = []
    for p, burnup in enumerate([0.0, 10.0]):
        densities = []
        for r, region in enumerate(regions):
            row = [0.0] * len(NAMES)
            if region["material"] == "fuel":
                row[NAMES.index("U-235")] = fuel[r][p]
                row[NAMES.index("U-238")] = 2.0
                row[NAMES.index("Pu-239")] = 0.1 * p
                row[NAMES.index("Sm-149")] = 0.01 * p
            else:
                row[NAMES.index("B-10")] = ifba[p]
            densities.append(row)
        points.append(
            {
                "burnup": burnup,
                "k_inf": k[p],
                "densities": densities,
                "pin_powers": pins[p],
                "few_group": None,
            }
        )
    return {
        "format": "xlc-depletion-history",
        "nuclides": NAMES,
        "regions": regions,
        "points": points,
    }


PINS = [[1.0, 1.0], [1.0, 0.0]]


def pair():
    """既定と厳密の履歴。10 の点で、k、燃料の 2 つ目の領域の U-235、B-10、ピン出力の 1 つが違う。"""
    exact = history([1.2, 1.1], [[2.0, 1.0], [2.0, 1.0]], [1.0, 0.5], [PINS, PINS])
    default = history(
        [1.2, 1.1 * (1 + 2e-6)],
        [[2.0, 1.0], [2.0, 1.0 + 4e-6]],
        [1.0, 0.5 * (1 - 1e-6)],
        [PINS, [[1.0, 1.0 + 3e-7], [1.0, 0.0]]],
    )
    return default, exact


class Summary(unittest.TestCase):
    def test_reactivity_difference_is_inverse_k_default_minus_exact_in_pcm(self):
        default, exact = pair()
        summary = compare.summarize(default, exact)
        self.assertEqual(summary["burnup"], [0.0, 10.0])
        self.assertEqual(summary["k"][0]["drho_pcm"], 0.0)
        expected = (1 / (1.1 * (1 + 2e-6)) - 1 / 1.1) * 1e5
        self.assertAlmostEqual(summary["k"][1]["drho_pcm"], expected, places=12)
        self.assertEqual(summary["max_drho_pcm"]["burnup"], 10.0)
        self.assertAlmostEqual(summary["max_drho_pcm"]["value"], expected, places=12)

    def test_densities_are_volume_weighted_per_kind_with_the_record_nuclides(self):
        default, exact = pair()
        summary = compare.summarize(default, exact)
        fuel = summary["kinds"]["fuel"]
        self.assertEqual(list(fuel), ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"])
        # 体積 1 と 3 の平均: (1 + 3 (1 + 4e-6)) / 4 = 1 + 3e-6
        self.assertAlmostEqual(fuel["U-235"]["series"][1], 3e-6, places=15)
        self.assertEqual(fuel["U-235"]["series"][0], 0.0)
        self.assertAlmostEqual(fuel["U-235"]["max"]["value"], 3e-6, places=15)
        self.assertEqual(fuel["U-235"]["max"]["burnup"], 10.0)
        # 領域の中の最大は、体積の平均より大きい 2 つ目の領域の差
        self.assertAlmostEqual(fuel["U-235"]["region_max"]["value"], 4e-6, places=15)
        # 両方 0 の点（0 の点の Pu-239、どの点の Xe-135）は比べない
        self.assertIsNone(fuel["Pu-239"]["series"][0])
        self.assertEqual(fuel["Pu-239"]["series"][1], 0.0)
        self.assertEqual(fuel["Xe-135"]["series"], [None, None])
        self.assertIsNone(fuel["Xe-135"]["max"])
        self.assertIsNone(fuel["Xe-135"]["region_max"])
        self.assertAlmostEqual(
            summary["kinds"]["ifba"]["B-10"]["max"]["value"], -1e-6, places=15
        )

    def test_the_largest_density_difference_names_its_kind_and_nuclide(self):
        default, exact = pair()
        largest = compare.summarize(default, exact)["max_density"]
        self.assertEqual(
            (largest["kind"], largest["nuclide"], largest["burnup"]),
            ("fuel", "U-235", 10.0),
        )
        self.assertAlmostEqual(largest["value"], 3e-6, places=15)

    def test_pin_powers_report_the_largest_relative_difference(self):
        default, exact = pair()
        summary = compare.summarize(default, exact)
        self.assertEqual(summary["pins"][0], 0.0)
        self.assertAlmostEqual(summary["pins"][1], 3e-7, places=15)
        self.assertAlmostEqual(summary["max_pins"]["value"], 3e-7, places=15)
        self.assertEqual(summary["max_pins"]["burnup"], 10.0)

    def test_mox_zones_use_the_fuel_nuclides(self):
        default, exact = pair()
        for h in (default, exact):
            h["regions"][0]["material"] = "mox-high"
            h["regions"][1]["material"] = "mox-low"
        summary = compare.summarize(default, exact)
        self.assertEqual(list(summary["kinds"]), ["mox-high", "mox-low", "ifba"])
        self.assertIn("Pu-239", summary["kinds"]["mox-low"])


class Mismatches(unittest.TestCase):
    def check_stops(self, change):
        default, exact = pair()
        change(default)
        with self.assertRaises(SystemExit):
            compare.summarize(default, exact)

    def test_other_burnup_points_stop(self):
        self.check_stops(lambda h: h["points"][1].update(burnup=11.0))

    def test_other_regions_stop(self):
        self.check_stops(lambda h: h["regions"][1].update(volume=2.0))
        self.check_stops(lambda h: h["regions"][2].update(material="pyrex"))

    def test_other_nuclides_stop(self):
        self.check_stops(lambda h: h.update(nuclides=list(reversed(NAMES))))

    def test_fission_in_a_cell_without_fission_on_one_side_stops(self):
        self.check_stops(lambda h: h["points"][1]["pin_powers"][1].__setitem__(1, 0.1))

    def test_a_density_zero_on_one_side_only_stops(self):
        self.check_stops(
            lambda h: h["points"][1]["densities"][0].__setitem__(
                NAMES.index("Xe-135"), 1e-9
            )
        )

    def test_few_group_constants_are_not_silently_skipped(self):
        self.check_stops(lambda h: h["points"][0].update(few_group={"total": [1.0]}))

    def test_other_formats_stop(self):
        self.check_stops(lambda h: h.update(format="xlc-other"))


class Files(unittest.TestCase):
    def test_json_lists_cases_with_input_md5_and_report_reads_it(self):
        default, exact = pair()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "d.msgpack").write_bytes(msgpack.packb(default))
            (tmp / "e.msgpack").write_bytes(msgpack.packb(copy.deepcopy(exact)))
            out = tmp / "out.json"
            compare.write_json(out, [("2A", tmp / "d.msgpack", tmp / "e.msgpack")])
            summary = json.loads(out.read_text())
            case = summary["2A"]
            self.assertEqual(case["inputs"]["default"]["path"], str(tmp / "d.msgpack"))
            self.assertEqual(len(case["inputs"]["exact"]["md5"]), 32)
            self.assertAlmostEqual(case["max_pins"]["value"], 3e-7, places=15)
            text = "\n".join(compare.report(summary))
            self.assertIn("2A", text)
            self.assertIn("pcm", text)


if __name__ == "__main__":
    unittest.main()
