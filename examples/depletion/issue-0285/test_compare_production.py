"""compare_production.py のログの読み方、状態の種類の分け方、テーブルの差（#0285、ADR-0250）。

実行: cd examples/depletion/issue-0285 && python3 -m unittest test_compare_production
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_production import breakdown, differences, kinds, log

LOG = (
    "t  INFO shielding:dancoff: xlc_lattice::shielding: close time.busy=2.20s time.idle=1µs\n"
    "t  INFO shielding: xlc_lattice::shielding: close time.busy=30.0s time.idle=1µs\n"
    "t  INFO ratio-shielding: xlc_lattice::shielding: close time.busy=20.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc::symmetric: symmetry tracks=1\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=6.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3000ms time.idle=1µs\n"
    't  INFO lattice: derivative state=3 nuclide="Pu-240" iterations=14 seconds=4.0\n'
    "t  INFO lattice: lattice state=3 iterations=25 tracks=1 seconds=70.0 total_seconds=75\n"
    "t  INFO ratio-shielding: xlc_lattice::shielding: close time.busy=15.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3.0s time.idle=1µs\n"
    "t  INFO lattice: lattice state=10 iterations=13 tracks=1 seconds=20.0 total_seconds=95\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3.0s time.idle=1µs\n"
    "t  INFO lattice: lattice state=17 iterations=13 tracks=1 seconds=8.0 total_seconds=103\n"
    "t  INFO lattice: branches shieldings=1 ratios=2 derivative_shieldings=0 total_seconds=104.5\n"
)


class Log(unittest.TestCase):
    def test_spans_are_summed_per_state_and_kinds_follow_the_shielding(self):
        parsed = log(LOG)
        first, ratio, remix = parsed["states"]
        # 自己遮蔽の中の区間（dancoff）は数えず、MOC の区間は ms も秒にする
        self.assertEqual(first["shielding"], [30.0])
        self.assertEqual(first["moc"], [6.0, 3.0])
        self.assertEqual(first["derivative"], [4.0])
        # 状態 3 と 10 は燃料の温度の点（7 点の 3）が同じで、比の基準の静止した標的を解くのは最初の状態 3
        self.assertEqual(kinds(parsed["states"], 7), ["ratio-first", "ratio", "remix"])
        self.assertEqual((ratio["state"], remix["seconds"]), (10, 8.0))

    def test_the_first_ratio_state_is_per_fuel_temperature_point(self):
        # 状態の番号を 7 で割った余りが燃料の温度の点: 3 と 10 は点 3、98 は点 0
        states = [
            {
                "state": n,
                "shielding": [],
                "ratio-shielding": r,
                "moc": [],
                "derivative": [],
            }
            for n, r in [(3, [1.0]), (14, [30.0]), (10, [1.0]), (98, [1.0]), (17, [])]
        ]
        states[1] |= {"shielding": [30.0], "ratio-shielding": []}
        self.assertEqual(
            kinds(states, 7), ["ratio-first", "ufg", "ratio", "ratio-first", "remix"]
        )

    def test_the_breakdown_counts_and_sums_by_kind(self):
        b = breakdown(log(LOG), 7)
        self.assertEqual(b["total_seconds"], 104.5)
        self.assertEqual(b["state_seconds"], 98.0)
        self.assertEqual((b["shieldings"], b["ratios"]), (1, 2))
        first = b["kinds"]["ratio-first"]
        self.assertEqual((first["ufg"], first["ratio"]), ([1, 30.0], [1, 20.0]))
        # 状態の最初の MOC は状態の固有値計算、ほかは微分の固有値計算
        self.assertEqual((first["moc_state"], first["moc_other"]), ([1, 6.0], [1, 3.0]))
        self.assertEqual(b["kinds"]["ratio"]["ratio"], [1, 15.0])
        self.assertEqual(b["kinds"]["remix"]["moc_state"], [1, 3.0])
        self.assertEqual(
            [(r["state"], r["kind"], r["ratio"]) for r in b["states"]],
            [(3, "ratio-first", 20.0), (10, "ratio", 15.0), (17, "remix", 0.0)],
        )

    def test_an_unfinished_run_stops(self):
        with self.assertRaises(SystemExit):
            log(LOG.rsplit("t  INFO lattice: branches", 1)[0])


def state(k: float, absorption: list) -> dict:
    return {
        "k_inf": k,
        "absorption": absorption,
        "nu_fission": [0.01, 0.1],
        "kappa_fission": [1.0, 2.0],
        "diffusion": [1.4, 0.4],
        "scatter": [[0.5, 0.02], [0.0, 1.0]],
        "adf": [1.0, 1.0],
    }


def micro(total: float) -> dict:
    return {
        "total": [total, 60.0],
        "nu_fission": [1.5, 0.08],
        "fission": [0.6, 0.03],
        "scatter": [[12.0, 0.01], [0.03, 15.0]],
        "thermal_fission": 0.03,
    }


def tab(points: list, states: list[dict], totals: list[list[float]]) -> dict:
    return {
        "points": points,
        "states": states,
        "fuels": 7,
        "micro": {
            "nuclides": ["U-235", "U-238"],
            "states": [[micro(t) for t in row] for row in totals],
        },
    }


class Differences(unittest.TestCase):
    POINTS = ((0, 2, 2, 3), (0, 2, 6, 3))

    def test_the_worst_state_is_named_by_its_point(self):
        base = tab(
            list(self.POINTS),
            [state(1.1, [0.01, 0.1]), state(0.8, [0.01, 0.1])],
            [[30.0, 10.0], [30.0, 10.0]],
        )
        variant = tab(
            list(self.POINTS),
            [state(1.1, [0.01, 0.1]), state(0.8002, [0.01, 0.101])],
            [[30.0, 10.0], [30.0, 10.1]],
        )
        d = differences(variant, base)
        self.assertAlmostEqual(d["worst_k"]["value"], 20.0)
        self.assertEqual(d["worst_k"]["point"], self.POINTS[1])
        # 反応度の差は 1/k の差で、k が 1 より小さいと k∞ の差より大きい
        self.assertAlmostEqual(d["worst_rho"]["value"], (1 / 0.8 - 1 / 0.8002) * 1e5)
        self.assertAlmostEqual(d["items"]["absorption"]["value"], 0.001 / 0.101)
        self.assertEqual(d["items"]["nu_fission"]["value"], 0.0)
        self.assertEqual(
            (d["micro"]["total"]["nuclide"], d["micro"]["total"]["point"]),
            ("U-238", self.POINTS[1]),
        )
        self.assertAlmostEqual(d["micro"]["total"]["value"], 0.1 / 10.1)

    def test_tables_of_other_states_stop(self):
        base = tab(list(self.POINTS), [state(1.1, [0.01, 0.1])] * 2, [[30.0, 10.0]] * 2)
        other = base | {"points": [self.POINTS[0], (0, 2, 5, 3)]}
        with self.assertRaises(SystemExit):
            differences(other, base)

    def test_a_micro_value_missing_from_both_tables_is_no_difference(self):
        # 核分裂しない核種の thermal_fission は無い
        base = tab(list(self.POINTS), [state(1.1, [0.01, 0.1])] * 2, [[30.0, 10.0]] * 2)
        for row in base["micro"]["states"]:
            row[1] = row[1] | {"thermal_fission": None}
        d = differences(base, base)
        self.assertEqual(d["micro"]["thermal_fission"]["value"], 0.0)
        one = base | {
            "micro": base["micro"] | {"states": [[micro(30.0), micro(10.0)]] * 2}
        }
        with self.assertRaises(SystemExit):
            differences(one, base)

    def test_values_that_are_not_finite_stop(self):
        base = tab(list(self.POINTS), [state(1.1, [0.01, 0.1])] * 2, [[30.0, 10.0]] * 2)
        broken = base | {
            "states": [state(1.1, [0.01, 0.1]), state(1.1, [0.01, math.nan])]
        }
        with self.assertRaises(SystemExit):
            differences(broken, base)


if __name__ == "__main__":
    unittest.main()
