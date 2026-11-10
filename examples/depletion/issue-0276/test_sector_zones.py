"""sector_zones.py（#0276 の燃料のピンを方位の扇形の組にも分けた燃焼の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0276 -p 'test_sector_zones.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sector_zones as sz


class PinPowerDifference(unittest.TestCase):
    def test_rms_and_largest_difference_cover_the_fuel_pins_only(self):
        base = [[1.0, 1.0], [2.0, 0.0]]
        other = [[1.01, 0.99], [1.96, 0.0]]
        d = sz.pin_power_difference(base, other)
        # 燃料の 3 本: +1%、−1%、−2%
        self.assertAlmostEqual(d["rms_percent"], (6 / 3) ** 0.5)
        self.assertAlmostEqual(d["max_percent"], -2.0)

    def test_power_in_a_pin_without_fuel_in_the_base_is_an_error(self):
        with self.assertRaises(SystemExit):
            sz.pin_power_difference([[1.0, 0.0]], [[1.0, 0.1]])

    def test_other_shapes_are_errors(self):
        with self.assertRaises(SystemExit):
            sz.pin_power_difference([[1.0, 1.0]], [[1.0], [1.0]])


class SectorSpread(unittest.TestCase):
    def test_spread_is_the_range_over_the_volume_weighted_mean(self):
        # 平均 (1 × 2 + 4 × 1) / 3 = 2
        self.assertAlmostEqual(sz.sector_spread([1.0, 4.0], [2.0, 1.0]), 150.0)

    def test_one_sector_orbit_has_no_spread(self):
        self.assertIsNone(sz.sector_spread([3.0], [1.0]))

    def test_none_in_the_mean_has_no_spread(self):
        self.assertIsNone(sz.sector_spread([0.0, 0.0], [1.0, 1.0]))


def region(material, cell, sector, volume, position):
    return {
        "material": material,
        "cell": cell,
        "cells": [position],
        "volume": volume,
        "sector": sector,
    }


class SectorRows(unittest.TestCase):
    def history(self):
        nuclides = ["U-235", "Pu-239"]
        regions = [
            region("fuel-0", 0, 0, 2.0, [0, 0]),
            region("fuel-0", 0, 1, 1.0, [0, 0]),
            region("fuel-1", 0, 0, 2.0, [0, 0]),
            region("fuel-1", 0, 1, 1.0, [0, 0]),
            region("fuel-0", 1, 0, 1.0, [1, 1]),
            region("fuel-1", 1, 0, 1.0, [1, 1]),
        ]
        densities = [
            [1.0, 1.0],
            [1.0, 4.0],
            [2.0, 2.0],
            [0.5, 2.0],
            [1.0, 1.0],
            [1.0, 1.0],
        ]
        return {
            "nuclides": nuclides,
            "regions": regions,
            "points": [{"burnup": 10.0, "densities": densities}],
        }

    def test_spreads_are_per_pin_orbit_and_ring_with_the_outer_ring_by_sector(self):
        (row,) = sz.sector_rows(self.history(), ["U-235", "Pu-239"])
        self.assertEqual(row["burnup"], 10.0)
        corner, centre = row["cells"]
        self.assertEqual(corner["position"], [0, 0])
        # 内のリングの Pu-239 は平均 2 に対し 1 と 4、外の U-235 は平均 1.5 に対し 2 と 0.5
        self.assertAlmostEqual(corner["ring_spread_percent"]["Pu-239"][0], 150.0)
        self.assertAlmostEqual(corner["ring_spread_percent"]["U-235"][1], 100.0)
        self.assertEqual(corner["ring_spread_percent"]["U-235"][0], 0.0)
        outer = corner["outer_by_sector_percent"]
        self.assertEqual(outer["sectors"], [0, 1])
        self.assertAlmostEqual(outer["U-235"][0], 100 * (2.0 / 1.5 - 1))
        self.assertAlmostEqual(outer["U-235"][1], 100 * (0.5 / 1.5 - 1))
        # 扇形の組が 1 つのピン（対称操作で扇形がすべて移り合う）は広がりを持たない
        self.assertEqual(centre["position"], [1, 1])
        self.assertIsNone(centre["ring_spread_percent"]["Pu-239"][1])

    def test_unsplit_histories_are_rejected(self):
        h = self.history()
        for r in h["regions"]:
            del r["sector"]
        with self.assertRaises(SystemExit):
            sz.sector_rows(h, ["U-235"])


class SpreadSummary(unittest.TestCase):
    def rows(self):
        def cell(position, spreads):
            return {"position": position, "ring_spread_percent": {"U-235": spreads}}

        return [
            {
                "burnup": 0.0,
                "cells": [cell([0, 0], [0.0, 0.0]), cell([1, 1], [None, None])],
            },
            {
                "burnup": 5.0,
                "cells": [cell([0, 0], [1.0, 3.0]), cell([1, 1], [None, None])],
            },
            {
                "burnup": 10.0,
                "cells": [cell([0, 0], [2.0, 2.5]), cell([1, 1], [None, None])],
            },
        ]

    def test_rows_are_the_selected_burnups_and_the_maximum_covers_every_point(self):
        s = sz.spread_summary("guide-rings5-sectors4", self.rows(), [0.0, 10.0])
        self.assertEqual([r["burnup"] for r in s["rows"]], [0.0, 10.0])
        self.assertEqual(
            s["max"]["U-235"],
            {"percent": 3.0, "burnup": 5.0, "position": [0, 0], "ring": 1},
        )


def point(burnup, k, pins=None):
    return {"burnup": burnup, "k": k, "pins": pins or [[1.0, 1.0]]}


class PairRows(unittest.TestCase):
    def test_rows_give_k_reactivity_and_pin_power_differences_against_the_base(self):
        base = {"path": "a", "points": [point(0, 1.0), point(1, 1.0)]}
        other = {
            "path": "b",
            "points": [point(0, 1.0), point(1, 1.001, [[1.002, 0.998]])],
        }
        rows = sz.pair_rows(base, other)
        self.assertEqual(rows[0]["dk_pcm"], 0.0)
        self.assertAlmostEqual(rows[1]["dk_pcm"], 100.0)
        # Δρ = 1/k_基準 − 1/k
        self.assertAlmostEqual(rows[1]["drho_pcm"], (1 - 1 / 1.001) * 1e5)
        self.assertAlmostEqual(rows[1]["pin_rms_percent"], 0.2)
        self.assertAlmostEqual(abs(rows[1]["pin_max_percent"]), 0.2)
        self.assertEqual(rows[1]["k_base"], 1.0)

    def test_histories_with_other_burnups_are_rejected(self):
        base = {"path": "a", "points": [point(0, 1.0)]}
        other = {"path": "b", "points": [point(0.1, 1.0)]}
        with self.assertRaises(SystemExit):
            sz.pair_rows(base, other)


class Report(unittest.TestCase):
    def test_extremes_cover_every_burnup_point(self):
        rows = [
            {
                "burnup": 0,
                "dk_pcm": 0.0,
                "drho_pcm": 0.0,
                "pin_rms_percent": 0.0,
                "pin_max_percent": 0.0,
            },
            {
                "burnup": 11,
                "dk_pcm": -9.94,
                "drho_pcm": -8.0,
                "pin_rms_percent": 0.05,
                "pin_max_percent": -0.08,
            },
            {
                "burnup": 40,
                "dk_pcm": 43.42,
                "drho_pcm": 52.87,
                "pin_rms_percent": 0.1,
                "pin_max_percent": 0.12,
            },
            {
                "burnup": 44,
                "dk_pcm": 42.3,
                "drho_pcm": 53.51,
                "pin_rms_percent": 0.09,
                "pin_max_percent": 0.11,
            },
        ]
        line = sz.extremes_line(rows)
        self.assertIn("最大 +43.42 pcm（40 MWd/kg）", line)
        self.assertIn("最小 −9.94 pcm（11 MWd/kg）", line)
        self.assertIn("Δρ +53.51 pcm（44 MWd/kg）", line)
        self.assertIn("rms の最大 0.100%（40 MWd/kg）", line)
        self.assertIn("絶対値が最大のピン出力の差 +0.120%（40 MWd/kg）", line)


if __name__ == "__main__":
    unittest.main()
