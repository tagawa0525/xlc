"""fuel_models.py の要約（燃料の温度の軸の、物理のモデルの式との比較と、当てはめと内挿の候補の誤差。#0274）。

実行: cd examples/depletion/issue-0274 && python3 -m unittest test_fuel_models
"""

import itertools
import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import compare_grid_interpolation as cgi
import fuel_models as fm

BURNUP = [0.0, 10.0]
BORON = [0.0, 325.0, 650.0, 975.0]
MODERATOR = [560.0, 570.0, 580.0, 590.0]
FUEL = [560.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0]


def state(absorption1: float, nu_fission2: float, scatter12: float) -> dict:
    """2 群の状態（上方散乱なし。全断面積は吸収と散乱の和）"""
    total1 = absorption1 + scatter12 + 0.47
    return {
        "k_inf": 1.0,
        "total": [total1, 1.0],
        "absorption": [absorption1, 0.08],
        "nu_fission": [0.005, nu_fission2],
        "fission": [0.002, nu_fission2 / 2.4],
        "chi": [1.0, 0.0],
        "diffusion": [1.4, 0.4],
        "scatter": [[0.47, scatter12], [0.0, 0.92]],
        "adf": [[1.0, 1.0]] * 6,
        "vacuum_ratio": [None] * 6,
        "normalized_adf": [False] * 6,
        "height": None,
        "pins": None,
    }


def table(g, fuel=FUEL, boron=BORON) -> dict:
    """直積の順（最初の軸がいちばん遅く変わる）で、各状態が g(燃焼度, ホウ素, 減速材, 燃料) の (Σa1, νΣf2, Σs12) の状態"""
    axes = [
        {"name": "burnup", "points": BURNUP},
        {
            "name": "boron",
            "points": boron,
            "b10": 0.199,
            "nominal": 325.0,
            "transform": "linear",
        },
        {
            "name": "moderator",
            "points": MODERATOR,
            "density": [0.75, 0.73, 0.71, 0.69],
            "nominal": 570.0,
            "transform": "linear",
        },
        {"name": "fuel", "points": fuel, "nominal": 900.0, "transform": "sqrt"},
    ]
    sizes = [len(a["points"]) for a in axes]
    return {
        "format": "xlc-table",
        "split": 1,
        "axes": axes,
        "points": [list(p) for p in itertools.product(*[range(n) for n in sizes])],
        "states": [
            state(*g(u, b, m, f))
            for u, b, m, f in itertools.product(BURNUP, boron, MODERATOR, fuel)
        ],
    }


def doppler(u, b, m, f):
    """高速群の吸収は √T、熱群の νΣf は T の 1 次（物理のモデルの 2 つの項）"""
    return (
        0.0093 * (1 + 0.004 * (1 + 1e-4 * b) * math.sqrt(f)),
        0.139 * (1 - 2e-6 * f) * (1 + 1e-4 * (m - 560)),
        0.0185 * (1 + 1e-3 * u),
    )


def jitter(u, b, m, f):
    """doppler に、状態ごとに符号の変わる √T の傾きを足したもの（係数が状態で跳ぶ）"""
    a, n, s = doppler(u, b, m, f)
    sign = (-1) ** (BORON.index(b) + MODERATOR.index(m))
    return a * (1 + 0.002 * sign * math.sqrt(f)), n, s


def design(summary, name):
    return next(d for d in summary["designs"] if d["design"] == name)


def scheme_max(summary, design_name, scheme):
    d = design(summary, design_name)
    row = next(s for s in d["schemes"] if s["scheme"] == list(scheme))
    return max(abs(p["worst_pcm"]) for c in row["centers"] for p in c["points"])


def exponent(summary, item):
    rows = next(e for e in summary["exponents"] if e["item"] == item)["burnups"]
    return [(r["min"], r["max"]) for r in rows]


def fit_item(summary, basis, item):
    f = next(x for x in summary["fits"] if x["basis"] == basis)
    return next(i for i in f["items"] if i["item"] == item)


class FitWeightsTest(unittest.TestCase):
    def test_the_square_root_and_linear_basis_reproduces_its_functions_between_points(
        self,
    ):
        q = 820.0
        weights = fm.fit_weights(fm.BASES["√T、T"], FUEL, q)
        for f in (
            lambda t: 1.0,
            math.sqrt,
            lambda t: t,
            lambda t: 3 - 2 * math.sqrt(t) + 1e-3 * t,
        ):
            value = sum(c * f(FUEL[k]) for k, c in weights.items())
            self.assertAlmostEqual(value, f(q), places=9)

    def test_a_fit_needs_as_many_points_as_terms(self):
        with self.assertRaises(ValueError):
            fm.fit_weights(fm.BASES["√T、T、T^1.5"], FUEL[:3], 700.0)


class ExponentTest(unittest.TestCase):
    def test_the_fitted_exponents_are_one_half_for_the_square_root_and_one_for_the_linear_term(
        self,
    ):
        summary = fm.summarize(table(doppler))
        for low, high in exponent(summary, "高速群の吸収"):
            self.assertEqual((low, high), (0.5, 0.5))
        for low, high in exponent(summary, "熱群の νΣf"):
            self.assertEqual((low, high), (1.0, 1.0))

    def test_items_that_do_not_change_with_the_temperature_have_no_exponent(self):
        summary = fm.summarize(table(doppler))
        rows = next(e for e in summary["exponents"] if e["item"] == "1→2 の散乱")[
            "burnups"
        ]
        for r in rows:
            self.assertEqual((r["min"], r["median"], r["max"]), (None, None, None))

    def test_exponents_at_the_end_of_the_grid_are_errors(self):
        ts = np.array([560.0, 750.0, 1000.0, 1350.0])
        with self.assertRaises(SystemExit):
            fm.best_exponents(ts, np.array([[1 + 1e-9 * t**3 for t in ts]]))

    def test_the_exponent_zero_is_the_logarithm(self):
        ts = [560.0, 750.0, 1000.0, 1350.0]
        p = fm.best_exponents(
            np.array(ts), np.array([[2 + 0.3 * math.log(t) for t in ts]])
        )
        self.assertEqual(p.tolist(), [0.0])


class FitTest(unittest.TestCase):
    def test_a_basis_holding_the_model_leaves_no_residual(self):
        summary = fm.summarize(table(doppler))
        for item in ("高速群の吸収", "熱群の νΣf"):
            i = fit_item(summary, "√T、T", item)
            self.assertLess(max(map(abs, i["residual_min"] + i["residual_max"])), 1e-6)
        f = next(x for x in summary["fits"] if x["basis"] == "√T、T")
        self.assertLess(
            max(abs(p["worst_pcm"]) for r in f["k"] for p in r["points"]), 1e-6
        )

    def test_the_square_root_alone_leaves_the_linear_term(self):
        summary = fm.summarize(table(doppler))
        i = fit_item(summary, "√T", "熱群の νΣf")
        self.assertGreater(max(map(abs, i["residual_min"] + i["residual_max"])), 1.0)

    def test_coefficients_that_follow_the_state_smoothly_are_not_rough(self):
        summary = fm.summarize(table(doppler))
        for c in fit_item(summary, "√T", "高速群の吸収")["coefficients"]:
            self.assertLess(c["coefficients"][0]["roughness"], 1e-2)

    def test_coefficients_that_jump_between_states_are_rough(self):
        summary = fm.summarize(table(jitter))
        for c in fit_item(summary, "√T", "高速群の吸収")["coefficients"]:
            self.assertGreater(c["coefficients"][0]["roughness"], 0.1)


class DesignTest(unittest.TestCase):
    def test_leaving_one_out_matches_the_grid_comparison(self):
        t = table(jitter)
        summary = fm.summarize(t)
        grid = cgi.summarize(
            t,
            [
                {"axis": "fuel", "variable": v, "method": m}
                for v, m in (("sqrt", "linear"), ("sqrt", "hermite"), ("log", "spline"))
            ],
        )
        for (variable, method), ours in zip(
            (("sqrt", "linear"), ("sqrt", "hermite"), ("log", "spline")),
            (
                ("interpolate", "√T", "linear"),
                ("interpolate", "√T", "hermite"),
                ("interpolate", "ln T", "spline"),
            ),
            strict=True,
        ):
            theirs = next(
                s
                for s in grid["schemes"]
                if (s["variable"], s["method"]) == (variable, method)
            )
            d = design(summary, "loo")
            row = next(s for s in d["schemes"] if s["scheme"] == list(ours))
            # cgi はホウ素と減速材と燃料の全部の組を見るので、最大どうしを比べる
            for a, b in zip(theirs["centers"], row["centers"], strict=True):
                self.assertEqual(a["point"], b["point"])
                for pa, pb in zip(a["points"], b["points"], strict=True):
                    self.assertAlmostEqual(pa["worst_pcm"], pb["worst_pcm"], places=9)

    def test_three_points_hold_the_model_through_the_square_root_parabola_and_fit(self):
        summary = fm.summarize(table(doppler))
        for scheme in (("interpolate", "√T", "hermite"), ("fit", "√T、T", "最小 2 乗")):
            self.assertLess(scheme_max(summary, "3 点", scheme), 1e-6)
        self.assertGreater(
            scheme_max(summary, "3 点", ("interpolate", "√T", "linear")), 1.0
        )

    def test_designs_skip_schemes_that_need_more_points(self):
        summary = fm.summarize(table(doppler))
        skipped = design(summary, "3 点")["skipped"]
        self.assertIn(["interpolate", "√T", "lagrange"], skipped)
        self.assertIn(["fit", "√T、T、T^1.5", "最小 2 乗"], skipped)
        self.assertNotIn(["fit", "√T、T", "最小 2 乗"], skipped)
        self.assertEqual(design(summary, "loo")["skipped"], [])


class SurrogateTest(unittest.TestCase):
    def test_the_square_root_hermite_holds_the_model_curve_between_points(self):
        summary = fm.summarize(table(doppler))
        for g in summary["surrogates"]:
            rows = {tuple(s["scheme"]): s for s in g["schemes"]}
            worst = {
                k: max(abs(p["worst_pcm"]) for c in s["centers"] for p in c["points"])
                for k, s in rows.items()
            }
            self.assertLess(worst[("interpolate", "√T", "hermite")], 1e-6)
            self.assertGreater(worst[("interpolate", "√T", "linear")], 0.1)
            self.assertEqual(
                len(rows[("interpolate", "√T", "linear")]["centers"]),
                len(g["points"]) - 1,
            )

    def test_fits_are_not_scored_against_their_own_curve(self):
        summary = fm.summarize(table(doppler))
        for g in summary["surrogates"]:
            own = ["fit", g["basis"], "最小 2 乗"]
            self.assertIn(own, g["skipped"])
            self.assertNotIn(own, [s["scheme"] for s in g["schemes"]])

    def test_no_interpolation_reproduces_a_surrogate_curve_exactly(self):
        # 4 点のラグランジュは変数の 3 次式を、3 次 Hermite は 2 次式を厳密に再現するので、真値の曲線がそれだと形の誤差が 0 になる
        def curved(u, b, m, f):
            a, n, s = doppler(u, b, m, f)
            return a * (1 + 0.05 / math.sqrt(f) + 2e-3 * math.log(f)), n, s

        summary = fm.summarize(table(curved))
        for g in summary["surrogates"]:
            for s in g["schemes"]:
                if s["scheme"][0] != "interpolate":
                    continue
                worst = max(
                    abs(p["worst_pcm"]) for c in s["centers"] for p in c["points"]
                )
                self.assertGreater(worst, 1e-6, (g["basis"], g["design"], s["scheme"]))

    def test_an_offset_at_one_point_enters_leaving_one_out_more_than_the_surrogate(
        self,
    ):
        def offset(u, b, m, f):
            a, n, s = doppler(u, b, m, f)
            return (a * (1 + 1e-4) if f == 1200.0 else a), n, s

        summary = fm.summarize(table(offset))
        hermite = ("interpolate", "√T", "hermite")
        loo = scheme_max(summary, "loo", hermite)
        for g in summary["surrogates"]:
            if g["design"] != "loo":
                continue
            row = next(s for s in g["schemes"] if tuple(s["scheme"]) == hermite)
            between = max(
                abs(p["worst_pcm"]) for c in row["centers"] for p in c["points"]
            )
            self.assertLess(between, loo / 10)


class InputTest(unittest.TestCase):
    def test_states_out_of_product_order_are_errors(self):
        t = table(doppler)
        t["points"][0], t["points"][1] = t["points"][1], t["points"][0]
        with self.assertRaises(SystemExit):
            fm.summarize(t)

    def test_tables_that_split_the_assembly_are_errors(self):
        t = table(doppler)
        t["split"] = 2
        with self.assertRaises(SystemExit):
            fm.summarize(t)

    def test_the_thinned_designs_need_the_seven_fuel_points(self):
        with self.assertRaises(SystemExit):
            fm.summarize(table(doppler, fuel=[560.0, 750.0, 900.0, 1050.0, 1350.0]))

    def test_the_roughness_needs_four_points_on_the_boron_and_moderator_axes(self):
        with self.assertRaises(SystemExit):
            fm.summarize(table(doppler, boron=[0.0, 325.0, 650.0]))

    def test_tables_without_the_product_axes_are_errors(self):
        t = table(doppler)
        t["axes"] = t["axes"][1:]
        with self.assertRaises(SystemExit):
            fm.summarize(t)

    def test_the_report_needs_the_inputs(self):
        with self.assertRaises(SystemExit):
            fm.report(fm.summarize(table(doppler)))


if __name__ == "__main__":
    unittest.main()
