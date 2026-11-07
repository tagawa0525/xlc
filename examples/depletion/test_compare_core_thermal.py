"""compare_core_thermal.py の要約（熱水力のある core の燃焼と、同じサイクルを入口の温度の一様な状態で解いた燃焼の差、熱水力の効きの
見込み、固有値計算の数、時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_thermal.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_thermal
from test_compare_core_guess import core, log


def table():
    """燃焼度 0、10、ホウ素 0、1000、減速材 565、600 K、燃料 565、600 K の、2 軸の組の部分格子（制御棒の値 out と aic）のテーブル。
    k∞ の反応度は、燃焼度で −0.001 /(MWd/kg)、ホウ素で −1e-4 /ppm、減速材で −2e-4 /K、燃料で −1e-3 /√K"""

    def k_inf(u, b, m, f):
        rho = (
            0.1
            - 0.001 * u
            - 1e-4 * b
            - 2e-4 * (m - 565.0)
            - 1e-3 * (math.sqrt(f) - math.sqrt(565.0))
        )
        return 1.0 / (1.0 - rho)

    burnup, boron, moderator, fuel = (
        [0.0, 10.0],
        [0.0, 1000.0],
        [565.0, 600.0],
        [565.0, 600.0],
    )
    points, states = [], []
    for iu, u in enumerate(burnup):
        for ib, b in enumerate(boron):
            for im, m in enumerate(moderator):
                for jf, f in enumerate(fuel):
                    # 3 軸の同時の分岐は持たない
                    if ib == 0 and im == 1 and jf == 1:
                        continue
                    for rod in (0, 1):
                        points.append([iu, ib, im, jf, rod])
                        states.append({"k_inf": k_inf(u, b, m, f) - 0.05 * rod})
    return {
        "axes": [
            {"name": "burnup", "points": burnup},
            {"name": "boron", "points": boron, "nominal": 1000.0},
            {"name": "moderator", "points": moderator, "nominal": 565.0},
            {"name": "fuel", "points": fuel, "nominal": 565.0},
            {"name": "rod", "values": ["out", "aic"]},
        ],
        "points": points,
        "states": states,
    }


def thermal_core(boron=1000.0, moderator=567.0, fuel=580.0, solves=(2, 3)):
    """`core` の燃焼に、熱水力の状態（入口の温度 565 K）と点のノードの温度を足したもの"""
    c = core(boron=boron, solves=solves)
    c["states"] = {"hfp": {"inlet": 565.0}}
    for p in c["points"]:
        p["burnup"] = p["efpd"] * 0.1
        p["temperatures"] = [[[moderator, fuel]], [[moderator, moderator]]]
    return c


def isothermal_core(boron=1000.0, solves=(1, 2)):
    c = core(boron=boron, solves=solves)
    c["states"] = {"hfp": {"moderator": 565.0, "fuel": 565.0}}
    for p in c["points"]:
        p["burnup"] = p["efpd"] * 0.1
    return c


def thermal_log(iterations, midpoints, total, thermal):
    """`log` に、点と中点の熱水力の行 (繰り返しの回数, 出口の温度, 燃料の温度の最大) を足したもの"""
    lines = [
        f"2026-11-06T00:00:00Z  INFO xlc_core::depletion: thermal-hydraulics iterations={n} change=1e-6 "
        f"outlet={o} fuel_average=570.0 fuel_max={f} density_error=0.001"
        for n, o, f in thermal
    ]
    return log(iterations, midpoints=midpoints, total=total) + "\n".join(lines) + "\n"


class Compare(unittest.TestCase):
    def setUp(self):
        self.iso = isothermal_core()
        self.iso_log = log([100] * 4, midpoints=((1.0, 900.0, 1),), total=5.0)
        self.th = thermal_core(boron=1000.0 - 30.0)
        self.th_log = thermal_log(
            [50] * 7,
            midpoints=((1.0, 900.0 - 28.0, 2),),
            total=8.0,
            thermal=((2, 566.0, 581.0), (2, 566.5, 582.0), (3, 567.0, 583.0)),
        )

    def test_differences_are_the_thermal_minus_the_isothermal_depletion(self):
        c = compare_core_thermal.compare(
            self.iso, self.iso_log, self.th, self.th_log, table()
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["boron"], -30.0, delta=1e-12)
        self.assertEqual(later["solves"], [2, 3])
        self.assertAlmostEqual(c["midpoints"][0]["boron"], -28.0, delta=1e-12)
        self.assertEqual(c["runs"]["thermal"]["solves"], 7)
        self.assertEqual(c["runs"]["isothermal"]["solves"], 4)
        self.assertEqual(
            c["thermal"],
            {"solves": 3, "iterations": [2, 3], "outlet": 567.0, "fuel": 583.0},
        )

    def test_the_estimate_takes_the_table_coefficients_at_the_power_weighted_temperatures(
        self,
    ):
        c = compare_core_thermal.compare(
            self.iso, self.iso_log, self.th, self.th_log, table()
        )
        later = c["points"][1]
        # 核分裂のあるノードだけが重みを持つ（2 体目は出力 0）。減速材 +2 K、燃料 √580 − √565
        self.assertAlmostEqual(later["moderator_rise"], 2.0, delta=1e-12)
        self.assertAlmostEqual(later["fuel_rise"], 15.0, delta=1e-12)
        fuel = -1e-3 * (math.sqrt(580.0) - math.sqrt(565.0)) / 1e-4
        moderator = -2e-4 * 2.0 / 1e-4
        self.assertAlmostEqual(later["estimate"]["fuel"], fuel, delta=1e-6)
        self.assertAlmostEqual(later["estimate"]["moderator"], moderator, delta=1e-6)
        self.assertAlmostEqual(later["estimate"]["total"], fuel + moderator, delta=1e-6)
        # 探索した点の熱水力の効きは、ホウ素の差（k は同じ）
        self.assertAlmostEqual(later["effect"], -30.0, delta=1e-6)
        self.assertEqual((later["moderator_max"], later["fuel_max"]), (567.0, 580.0))

    def test_points_at_a_given_boron_count_the_k_difference_as_boron(self):
        # ホウ素を与えた点（同じ濃度）は、k の差の反応度を k∞ のホウ素の係数でホウ素に直す
        th = thermal_core()
        th["points"][0]["k"] = 1.0 / (1.0 + 1e-3)
        c = compare_core_thermal.compare(
            self.iso, self.iso_log, th, self.th_log, table()
        )
        first = c["points"][0]
        self.assertAlmostEqual(first["effect"], -1e-3 / 1e-4, delta=1e-6)

    def test_mismatched_states_and_points_are_errors(self):
        th = thermal_core()
        th["states"]["hfp"]["inlet"] = 560.0
        with self.assertRaises(SystemExit):
            compare_core_thermal.compare(
                self.iso, self.iso_log, th, self.th_log, table()
            )
        th = thermal_core()
        th["points"].pop()
        with self.assertRaises(SystemExit):
            compare_core_thermal.compare(
                self.iso, self.iso_log, th, self.th_log, table()
            )
        # 熱水力の燃焼の点はノードの温度を持つ
        th = thermal_core()
        del th["points"][1]["temperatures"]
        with self.assertRaises(SystemExit):
            compare_core_thermal.compare(
                self.iso, self.iso_log, th, self.th_log, table()
            )


if __name__ == "__main__":
    unittest.main()
