"""compare_burnup_interpolation.py の要約（燃焼度の軸の内挿の方式ごとの、間引いた点での k∞ の誤差とオーバーシュート。ADR-0192）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_burnup_interpolation.py'
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_burnup_interpolation as cbi


def state(nu_fission2: float) -> dict:
    """2 群の状態。上方散乱なし、k∞ = (νΣf1 + νΣf2 Σs12 / Σr2) / Σr1 で、νΣf2 に比例する部分を持つ"""
    return {
        "k_inf": 0.0,
        "total": [0.5, 1.0],
        "absorption": [0.01, 0.1],
        "nu_fission": [0.005, nu_fission2],
        "fission": [0.002, 0.04],
        "chi": [1.0, 0.0],
        "diffusion": [1.4, 0.4],
        "scatter": [[0.47, 0.02], [0.0, 0.9]],
        "adf": [[1.0, 1.0]] * 6,
        "vacuum_ratio": [None] * 6,
        "normalized_adf": [False] * 6,
        "height": None,
        "pins": None,
    }


def k_of(nu_fission2: float) -> float:
    # Σr1 = 0.5 − 0.47 = 0.03、Σr2 = 1.0 − 0.9 = 0.1、Σs12 = 0.02
    return (0.005 + nu_fission2 * 0.02 / 0.1) / 0.03


def table(points: list[float], g) -> dict:
    """燃焼度の軸と、点が 1 つのほかの軸。各点の熱群の νΣf が g(燃焼度)"""
    return {
        "format": "xlc-table",
        "split": 1,
        "axes": [
            {"name": "burnup", "points": points},
            {"name": "boron", "points": [1300.0], "b10": 0.199},
            {"name": "moderator", "points": [565.0], "density": [0.743]},
            {"name": "fuel", "points": [565.0]},
        ],
        "states": [{**state(g(b)), "k_inf": k_of(g(b))} for b in points],
    }


def quadratic(b: float) -> float:
    """燃焼度 0 だけ曲線の外（Xe と Sm のない燃料）の 2 次"""
    return 0.2 if b == 0.0 else 0.1 + 0.004 * b - 1e-4 * b * b


class Interpolation(unittest.TestCase):
    def test_hermite_reproduces_quadratics_past_the_first_point(self):
        x = np.array([0.0, 0.5, 2.0, 5.0, 10.0, 18.0])
        y = np.array([[quadratic(b)] for b in x])
        q = np.array([0.5, 0.7, 3.1, 7.5, 17.9])
        got = cbi.interpolate("hermite", x, y, q)[:, 0]
        np.testing.assert_allclose(got, [quadratic(b) for b in q], rtol=1e-13)

    def test_the_first_span_from_zero_is_linear_and_zero_does_not_shape_the_rest(self):
        x = np.array([0.0, 0.5, 2.0, 5.0])
        y = np.array([[quadratic(b)] for b in x])
        got = cbi.interpolate("hermite", x, y, np.array([0.2]))[0, 0]
        self.assertAlmostEqual(got, 0.6 * 0.2 + 0.4 * quadratic(0.5), places=15)
        other = y.copy()
        other[0, 0] = 1.0
        q = np.array([0.6, 1.0, 4.0])
        np.testing.assert_array_equal(
            cbi.interpolate("hermite", x, y, q), cbi.interpolate("hermite", x, other, q)
        )
        # 0 も傾きに使う方式は、0 の値で先の区間が変わる
        self.assertFalse(
            np.array_equal(
                cbi.interpolate("hermite-with-zero", x, y, q),
                cbi.interpolate("hermite-with-zero", x, other, q),
            )
        )

    def test_every_method_returns_the_grid_points(self):
        x = np.array([0.0, 0.5, 2.0, 5.0, 10.0, 18.0])
        y = np.array([[quadratic(b), 2.0 * quadratic(b)] for b in x])
        for method in cbi.METHODS:
            np.testing.assert_allclose(
                cbi.interpolate(method, x, y, x), y, rtol=1e-14, err_msg=method
            )

    def test_unknown_methods_are_errors(self):
        x = np.array([0.0, 1.0, 2.0])
        with self.assertRaises(SystemExit):
            cbi.interpolate("spline", x, x[:, None], x)


class KInf(unittest.TestCase):
    def test_two_group_k_infinity_without_upscatter(self):
        self.assertAlmostEqual(cbi.k_inf(state(0.1)), k_of(0.1), places=13)


class Summary(unittest.TestCase):
    def test_removed_points_take_the_full_table_as_the_truth(self):
        points = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        s = cbi.summarize_table(
            table(points, quadratic), [{"name": "half", "burnup": [0.0, 1.0, 3.0, 5.0]}]
        )
        (sel,) = s["selections"]
        self.assertEqual(sel["removed"], [2.0, 4.0])
        hermite = sel["results"]["hermite"]
        # k は νΣf2 について 1 次なので、2 次の νΣf2 は Hermite でちょうど再現する
        self.assertLess(hermite["max_abs_pcm"], 1e-9)
        linear = sel["results"]["linear"]
        # 線形の誤差: νΣf2 が h²/8 × 2e-4 = 1e-4 小さい（h = 2）
        want = [
            (1 / k_of(quadratic(b)) - 1 / k_of(quadratic(b) - 1e-4)) * 1e5
            for b in (2.0, 4.0)
        ]
        np.testing.assert_allclose(linear["rho_pcm"], want, rtol=1e-6)
        self.assertAlmostEqual(linear["max_abs_pcm"], max(abs(w) for w in want))
        self.assertAlmostEqual(
            linear["rms_pcm"], float(np.sqrt(np.mean(np.square(want))))
        )
        self.assertLess(s["table_k_max_abs_pcm"], 1e-9)

    def test_overshoot_counts_only_items_whose_truth_stays_in_the_bracket(self):
        # 熱群の νΣf2 が 10 で折れ、傾きが 0 から 0.005 に変わる（Gd の燃え尽きを模す）。間隔 2
        points = [2.0 * i for i in range(1, 11)]
        g = lambda b: 0.1 + 0.005 * max(b - 10.0, 0.0)
        selection = points[::2] + [points[-1]]
        s = cbi.summarize_table(
            table(points, g), [{"name": "every other", "burnup": selection}]
        )
        (sel,) = s["selections"]
        # 間引いた軸は間隔 4。折れ曲がりの手前の区間 [6, 10] で、Hermite は 8 の点でデータの範囲 0.1 を下回る
        # （10 の傾きは 0.0025、8 は区間の t = 1/2 で、H11 × h × 傾き = −1/8 × 4 × 0.0025）
        self.assertAlmostEqual(
            sel["results"]["hermite"]["overshoot"], 1.25e-3 / 0.1, places=12
        )
        self.assertEqual(sel["results"]["linear"]["overshoot"], 0.0)
        self.assertEqual(sel["results"]["pchip"]["overshoot"], 0.0)

    def test_tables_with_more_than_one_state_are_errors(self):
        t = table([0.0, 1.0, 2.0], quadratic)
        t["axes"][1]["points"] = [1000.0, 1300.0]
        with self.assertRaises(SystemExit):
            cbi.summarize_table(t, [{"name": "all", "burnup": [0.0, 1.0, 2.0]}])

    def test_selections_must_be_table_points_from_the_first_to_the_last(self):
        t = table([0.0, 1.0, 2.0, 3.0], quadratic)
        for burnup in ([0.0, 1.5, 3.0], [0.0, 1.0, 2.0], [1.0, 2.0, 3.0]):
            with self.assertRaises(SystemExit, msg=str(burnup)):
                cbi.summarize_table(t, [{"name": "bad", "burnup": burnup}])


class Report(unittest.TestCase):
    def test_report_names_the_worst_point_of_each_method(self):
        points = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        summary = {
            "format": "xlc-compare-burnup-interpolation",
            "inputs": [{"path": "x.msgpack", "md5": "0" * 32}],
            "tables": [
                {
                    "name": "T",
                    **cbi.summarize_table(
                        table(points, quadratic),
                        [{"name": "half", "burnup": [0.0, 1.0, 3.0, 5.0]}],
                    ),
                }
            ],
        }
        text = "\n".join(cbi.report(summary))
        self.assertIn("T half", text)
        self.assertIn("linear", text)


if __name__ == "__main__":
    unittest.main()
