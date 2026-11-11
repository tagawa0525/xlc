"""fine_pin.py の要約（燃料の温度を細かく取ったピンセルの、当てはめの残差の跳びと、格子点から引いた足した点の誤差。#0274）。

実行: cd examples/depletion/issue-0274/fine-pin && python3 -m unittest test_fine_pin
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import fine_pin as fp
import test_fuel_models as tfm

FINE = [560.0] + [600.0 + 50.0 * n for n in range(16)]


def pin(g, fuel=FINE) -> dict:
    """燃料の温度の軸だけを持つピンセルのテーブル。各点の状態は g(燃料の温度) の (Σa1, νΣf2, Σs12)"""
    return {
        "format": "xlc-table",
        "axes": [
            {"name": "boron", "points": [650.0], "nominal": 650.0},
            {"name": "moderator", "points": [580.0], "nominal": 580.0},
            {"name": "fuel", "points": fuel, "nominal": 900.0},
        ],
        "states": [tfm.state(*g(t)) for t in fuel],
    }


def model(t):
    """物理のモデルの項: 高速群の吸収は √T、熱群の νΣf は T の 1 次"""
    return (0.01 + 2e-5 * math.sqrt(t), 0.12 - 1e-6 * t, 0.017)


def bias_at(points, size):
    """model に、points の温度だけ高速群の吸収を相対で size 動かしたもの"""

    def g(t):
        a, n, s = model(t)
        return (a * (1 + size) if t in points else a, n, s)

    return g


def curved(t):
    """model に、√T、T、T^1.5 で表せない滑らかな曲がり（ln T の 2 乗）を足したもの"""
    a, n, s = model(t)
    return (a * (1 + 1e-4 * math.log(t / 900.0) ** 2), n, s)


def fit(summary, basis):
    return next(f for f in summary["fits"] if f["basis"] == basis)


def prediction(summary, scheme):
    return next(p for p in summary["predictions"] if p["scheme"] == list(scheme))


class FitTest(unittest.TestCase):
    def test_a_model_inside_the_basis_leaves_no_residual_or_jump(self):
        f = fit(fp.summarize_pin(pin(model)), "√T、T")
        for values in [f["all"]["高速群の吸収"], f["all"]["k"], f["jumps"]["k"]]:
            for x in values:
                self.assertAlmostEqual(x, 0.0, places=6)

    def test_a_bias_on_the_grid_points_jumps_there_and_not_between(self):
        f = fit(fp.summarize_pin(pin(bias_at([1200.0], 5e-5))), "√T、T、T^1.5")
        jumps = dict(zip(FINE[1:-1], f["jumps"]["高速群の吸収"], strict=True))
        self.assertGreater(abs(jumps[1200.0]), 3.0)
        # 隣の点は片寄りの半分ほど逆に跳び、離れた点はほとんど跳ばない
        for t in [700.0, 800.0, 900.0]:
            self.assertLess(abs(jumps[t]), 0.3)

    def test_a_smooth_curve_outside_the_basis_jumps_much_less_than_its_residual(self):
        f = fit(fp.summarize_pin(pin(curved)), "√T")
        residual = max(abs(x) for x in f["all"]["高速群の吸収"])
        jump = max(abs(x) for x in f["jumps"]["高速群の吸収"])
        self.assertGreater(residual, 1.0)
        self.assertLess(jump, residual / 10)

    def test_the_grid_fit_uses_only_the_grid_points(self):
        # 足した点だけを動かしても、格子点の当てはめは変わらない
        added = [t for t in FINE if t not in fp.GRID]
        a = fit(fp.summarize_pin(pin(model)), "√T、T")["grid"]["k"]
        b = fit(fp.summarize_pin(pin(bias_at(added, 1e-3))), "√T、T")["grid"]["k"]
        for x, y in zip(a, b, strict=True):
            self.assertAlmostEqual(x, y, places=9)


class PredictionTest(unittest.TestCase):
    def test_schemes_predict_the_added_points_from_the_grid_points(self):
        s = fp.summarize_pin(pin(model))
        self.assertEqual(s["added"], [t for t in FINE if t not in fp.GRID])
        # √T の 3 次 Hermite は a + b√T + cT を厳密に再現する
        exact = prediction(s, ("interpolate", "√T", "hermite"))["pcm"]
        self.assertEqual(len(exact), 10)
        for x in exact:
            self.assertAlmostEqual(x, 0.0, places=6)
        # √T の線形は T の項を再現しない
        linear = prediction(s, ("interpolate", "√T", "linear"))["pcm"]
        self.assertGreater(max(abs(x) for x in linear), 0.01)

    def test_a_bias_on_an_added_point_shows_in_its_own_prediction_only(self):
        s = fp.summarize_pin(pin(bias_at([1100.0], 1e-4)))
        rows = prediction(s, ("interpolate", "√T", "hermite"))["pcm"]
        by_point = dict(zip(s["added"], rows, strict=True))
        self.assertGreater(abs(by_point[1100.0]), 0.1)
        for t, x in by_point.items():
            if t != 1100.0:
                self.assertAlmostEqual(x, 0.0, places=6)


class InputTest(unittest.TestCase):
    def test_a_pin_table_without_the_grid_points_is_an_error(self):
        with self.assertRaises(SystemExit):
            fp.summarize_pin(pin(model, fuel=[560.0, 650.0, 750.0, 900.0, 1050.0]))

    def test_a_pin_table_with_two_boron_points_is_an_error(self):
        t = pin(model)
        t["axes"][0]["points"] = [0.0, 650.0]
        with self.assertRaises(SystemExit):
            fp.summarize_pin(t)

    def test_a_pin_table_with_only_the_grid_points_is_an_error(self):
        with self.assertRaises(SystemExit):
            fp.summarize_pin(pin(model, fuel=fp.GRID))

    def test_a_pin_table_with_a_point_outside_the_grid_range_is_an_error(self):
        # 格子点の範囲の外の点は、方式で引くと外挿になる
        with self.assertRaises(SystemExit):
            fp.summarize_pin(pin(model, fuel=[500.0, *FINE]))


class ProductTest(unittest.TestCase):
    def test_the_product_fit_takes_the_nominal_state_of_the_burnup_point(self):
        def g(u, b, m, f):
            a, n, s = model(f)
            # 公称の状態（ホウ素 325、減速材 570）の 10 MWd/kg だけ 1200 K が片寄る
            hit = (u, b, m, f) == (10.0, 325.0, 570.0, 1200.0)
            return (a * (1 + 5e-5) if hit else a, n, s)

        product = tfm.table(g)
        at = {
            u: fit(
                fp.summarize_product(
                    fp.load_product(product), u, {"boron": 325.0, "moderator": 570.0}
                ),
                "√T、T",
            )["grid"]["高速群の吸収"]
            for u in [0.0, 10.0]
        }
        self.assertTrue(all(abs(x) < 1e-6 for x in at[0.0]))
        self.assertGreater(abs(at[10.0][fp.GRID.index(1200.0)]), 1.0)

    def test_a_pin_state_off_the_product_nominal_is_an_error(self):
        product = fp.load_product(tfm.table(lambda u, b, m, f: model(f)))
        with self.assertRaises(SystemExit):
            fp.summarize_product(product, 0.0, {"boron": 650.0, "moderator": 570.0})


class ReferenceTest(unittest.TestCase):
    def test_the_same_states_on_the_grid_differ_by_nothing(self):
        self.assertEqual(
            fp.reference_difference(pin(model), pin(model, fuel=fp.GRID)), 0.0
        )

    def test_a_changed_grid_state_shows_its_relative_difference(self):
        diff = fp.reference_difference(pin(model), pin(bias_at([900.0], 1e-6), fp.GRID))
        self.assertAlmostEqual(diff, 1e-6, delta=1e-9)


class LibraryTest(unittest.TestCase):
    def test_only_the_shared_nuclide_temperatures_are_compared(self):
        a = {
            "nuclides": [
                {"name": "U-238", "temperatures": [{"temperature": 900.0, "x": [1]}]},
                {"name": "O-16", "temperatures": [{"temperature": 950.0, "x": [2]}]},
            ]
        }
        b = {
            "nuclides": [
                {"name": "U-238", "temperatures": [{"temperature": 900.0, "x": [3]}]},
                {"name": "O-16", "temperatures": [{"temperature": 900.0, "x": [2]}]},
            ]
        }
        self.assertEqual(
            fp.library_difference(a, b), {"shared": 1, "differ": [["U-238", 900.0]]}
        )


class RoundedTest(unittest.TestCase):
    def test_floats_keep_six_significant_digits_and_other_values_stay(self):
        got = fp.rounded(
            {"pcm": [1.23456789e-3, -98765.4321, 0.0], "n": 17, "md5": "8efb40eb"}
        )
        self.assertEqual(
            got, {"pcm": [1.23457e-3, -98765.4, 0.0], "n": 17, "md5": "8efb40eb"}
        )

    def test_k_inf_keeps_every_digit(self):
        # k∞ の 6 桁は 1 より大きいと小数 5 桁で、反応度の差の pcm の小数 1 桁がずれる
        k = [1.2531162345678, 0.8956091234567]
        self.assertEqual(
            fp.rounded({"pins": [{"k_inf": k, "x": 1.2531162345678}]}),
            {"pins": [{"k_inf": k, "x": 1.25312}]},
        )


if __name__ == "__main__":
    unittest.main()
