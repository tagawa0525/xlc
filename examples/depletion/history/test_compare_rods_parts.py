"""compare_rods_parts.py（制御棒を途中まで入れて抜く燃焼の core − 参照を、数密度の効きと定数の側に分ける）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_rods_parts.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_forgetting as cf
import compare_rods_parts as c
from compare_residual import added, coefficients, corrected, k_inf, rho
from compare_rods_switch import production
from test_compare_rods_switch import few, residual

NUCLIDES = ["U-235", "U-238", "Pu-240"]
BURNUPS = [0.0, 10.0, 20.0]
NOMINAL = [[7.0e-4, 2.2e-2, 0.0], [6.0e-4, 2.19e-2, 1.0e-5], [5.0e-4, 2.18e-2, 2.0e-5]]
# 制御棒を入れた履歴の分岐の数密度（制御棒の値 out の R の点）
BRANCH = [NOMINAL[0], [5.8e-4, 2.19e-2, 1.4e-5], [4.7e-4, 2.18e-2, 2.9e-5]]
# 制御棒を抜いた燃焼度
SWITCH = 10.0


def sigma(s):
    return {
        "total": [11.0 + s, 130.0 + s],
        "nu_fission": [1.1, 170.0 + s],
        "fission": [0.4, 70.0],
        "scatter": [[9.0, 0.1], [0.01, 8.0]],
    }


def table():
    """燃焼度 0、10、20、制御棒の軸（out、aic）、追う核種 3 つの公称のテーブルと、2 つの変数の補正"""
    axes = [
        {"name": "burnup", "points": BURNUPS},
        {"name": "boron", "points": [1300.0], "nominal": 1300.0},
        {"name": "moderator", "points": [565.0, 600.0], "nominal": 565.0},
        {"name": "fuel", "points": [565.0], "nominal": 565.0},
        {"name": "rod", "values": ["out", "aic"]},
    ]
    points = [[b, 0, m, 0, r] for b in range(3) for m in range(2) for r in range(2)]

    def rod(value, scale):
        return {
            "rod": value,
            "burnup": BURNUPS,
            "untracked": [residual(0.0), residual(0.3 * scale), residual(0.5 * scale)],
            "linearization": [residual(0.0), residual(scale), residual(2.0 * scale)],
            "densities": BRANCH,
        }

    return {
        "axes": axes,
        "points": points,
        "states": [few(float(i)) for i in range(len(points))],
        "micro": {
            "nuclides": NUCLIDES,
            "densities": NOMINAL,
            "states": [
                [sigma(float(i)), sigma(0.5 * i), sigma(2.0 * i)]
                for i in range(len(points))
            ],
        },
        "histories": [
            {
                "variable": "moderator-density",
                "nominal": 0.743,
                "branch": 0.661,
                "rods": [rod("out", 0.7)],
            },
            {
                "variable": "rods",
                "nominal": 0.0,
                "branch": 1.0,
                "rods": [rod("aic", -1.0), rod("out", 1.0)],
            },
        ],
    }


def expected_node(t, b, densities, w):
    """燃焼度の点 `b`（格子点）のノードの定数（core の `node_constants` の式を書き下したもの）"""
    p = t["points"].index([b, 0, 0, 0, 0])
    factors = coefficients(t["micro"]["states"][p])
    nominal = t["micro"]["densities"][b]
    deltas = [x - y for x, y in zip(densities, nominal)]
    f = corrected(t["states"][p], factors, deltas)
    out = t["histories"][1]["rods"][1]
    branch = [x - y for x, y in zip(out["densities"][b], nominal)]
    f = added(f, c.scaled(out["untracked"][b], w))
    # 履歴の分岐の数密度が公称と同じ点（燃焼度 0）は L を足さない（core の `production_ratio`）
    if all(d == 0.0 for d in branch):
        return f, factors, None
    a = production(factors, deltas) / production(factors, branch)
    return added(f, c.scaled(out["linearization"][b], a)), factors, a


# core と参照の数密度（燃焼度の点 0、10、20）
CORE_DENSITIES = [NOMINAL[0], [5.9e-4, 2.19e-2, 1.3e-5], [4.9e-4, 2.18e-2, 2.6e-5]]
REFERENCE_DENSITIES = [NOMINAL[0], [5.9e-4, 2.19e-2, 1.3e-5], [4.8e-4, 2.18e-2, 2.7e-5]]
# 参照の k∞ の燃焼度 0 と 10 の値（20 は `reference_k`）
REFERENCE_K = [1.2, 1.1]
# 燃焼度 20 の参照に合う U(1) の重み（core の w は 0.5）
REFERENCE_WEIGHT = 0.3


def reference_k(t):
    """燃焼度 20 の参照の k∞: core のノードの定数の数密度を参照のものにし、U(1) の重みを `REFERENCE_WEIGHT` にした定数の 2 群の k"""
    f, factors, _ = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
    at_reference = corrected(
        f, factors, [y - x for x, y in zip(CORE_DENSITIES[2], REFERENCE_DENSITIES[2])]
    )
    u = t["histories"][1]["rods"][1]["untracked"][2]
    return k_inf(added(at_reference, c.scaled(u, REFERENCE_WEIGHT - 0.5)))


def histories():
    """表 `table` で解いたことにする core の燃焼の履歴（1 ノード）と、格子計算の参照の燃焼の履歴"""
    t = table()
    h = [1.0, 1.0, SWITCH / 20.0]
    points = []
    for b, u in enumerate(BURNUPS):
        f, _, _ = expected_node(t, b, CORE_DENSITIES[b], h[b])
        points.append(
            {
                "burnups": [[u]],
                "state": "rods" if u < SWITCH else "out",
                "k": k_inf(f),
                "densities": [[CORE_DENSITIES[b]]],
                "histories": [[[0.743, h[b]]]],
            }
        )
    core = {
        "format": "xlc-core-history",
        "assemblies": ["A01"],
        "axial_mesh": [0.0, 21.5],
        "nuclides": NUCLIDES,
        "history_variables": ["moderator-density", "rods"],
        "points": points,
    }
    reference = {
        "nuclides": NUCLIDES,
        "regions": [{"volume": 1.0}],
        "points": [
            {"burnup": u, "k_inf": k, "densities": [d]}
            for u, k, d in zip(
                BURNUPS, [*REFERENCE_K, reference_k(t)], REFERENCE_DENSITIES
            )
        ],
        "continued": {"burnup": SWITCH},
    }
    return t, reference, core


class BurnupTerms(unittest.TestCase):
    def test_grid_points_take_the_point_itself(self):
        self.assertEqual(c.burnup_terms([0.0, 1.0, 2.0, 4.0], 2.0), [(2, 1.0)])

    def test_the_first_interval_from_zero_is_linear(self):
        # 燃焼度 0 は曲線に使わない（core の burnup_terms）
        terms = c.burnup_terms([0.0, 1.0, 2.0, 4.0], 0.25)
        self.assertEqual(terms, [(0, 0.75), (1, 0.25)])

    def test_the_curve_reproduces_a_parabola(self):
        # 3 点の放物線の傾きの 3 次 Hermite は 2 次式を厳密に通す
        points = [0.0, 1.0, 2.0, 4.0, 7.0]
        for x in (1.5, 3.0, 5.5):
            value = sum(w * points[i] ** 2 for i, w in c.burnup_terms(points, x))
            self.assertAlmostEqual(value, x * x, delta=1e-12)

    def test_burnups_outside_the_points_are_errors(self):
        with self.assertRaises(SystemExit):
            c.burnup_terms([0.0, 1.0, 2.0], 2.5)


class Compare(unittest.TestCase):
    def compare(self):
        return c.compare(*histories())

    def test_rows_are_the_points_after_the_switch(self):
        result = self.compare()
        self.assertEqual(result["switch"], SWITCH)
        self.assertEqual([r["burnup"] for r in result["rows"]], [20.0])

    def test_the_node_constants_are_those_of_the_core(self):
        t, reference, core = histories()
        (row,) = c.compare(t, reference, core)["rows"]
        f, _, a = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
        self.assertAlmostEqual(row["weight"], 0.5, delta=1e-15)
        self.assertAlmostEqual(row["index_ratio"], a, delta=1e-12)
        self.assertAlmostEqual(row["reconstruct_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(
            row["difference_pcm"], rho(k_inf(f)) - rho(reference_k(t)), delta=1e-9
        )

    def test_a_core_k_other_than_its_constants_is_an_error(self):
        t, reference, core = histories()
        core["points"][2]["k"] *= 1.001
        with self.assertRaises(SystemExit):
            c.compare(t, reference, core)

    def test_the_composition_effect_removes_the_density_differences(self):
        t, reference, core = histories()
        (row,) = c.compare(t, reference, core)["rows"]
        f, factors, _ = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
        at_reference = corrected(
            f,
            factors,
            [y - x for x, y in zip(CORE_DENSITIES[2], REFERENCE_DENSITIES[2])],
        )
        self.assertAlmostEqual(
            row["composition_pcm"],
            rho(k_inf(f)) - rho(k_inf(at_reference)),
            delta=1e-9,
        )
        self.assertAlmostEqual(
            row["constant_pcm"],
            row["difference_pcm"] - row["composition_pcm"],
            delta=1e-12,
        )

    def test_equal_densities_have_no_composition_effect(self):
        t, reference, core = histories()
        reference["points"][2]["densities"] = [CORE_DENSITIES[2]]
        (row,) = c.compare(t, reference, core)["rows"]
        self.assertEqual(row["composition_pcm"], 0.0)
        for n in row["nuclides"]:
            self.assertEqual(n["percent"], 0.0)
            self.assertEqual(n["pcm"], 0.0)

    def test_nuclide_rows_hold_every_tracked_nuclide(self):
        (row,) = self.compare()["rows"]
        self.assertEqual([n["nuclide"] for n in row["nuclides"]], NUCLIDES)
        u235 = row["nuclides"][0]
        self.assertAlmostEqual(
            u235["percent"],
            100.0 * (CORE_DENSITIES[2][0] / REFERENCE_DENSITIES[2][0] - 1.0),
            delta=1e-12,
        )

    def test_reference_weights_bring_the_node_to_the_reference(self):
        t, reference, core = histories()
        (row,) = c.compare(t, reference, core)["rows"]
        f, factors, _ = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
        at_reference = corrected(
            f,
            factors,
            [y - x for x, y in zip(CORE_DENSITIES[2], REFERENCE_DENSITIES[2])],
        )
        self.assertAlmostEqual(row["untracked_weight"], REFERENCE_WEIGHT, delta=1e-9)
        out = t["histories"][1]["rods"][1]
        u, l = out["untracked"][2], out["linearization"][2]
        for key, part, now in (
            ("untracked_weight", u, row["weight"]),
            ("linearization_weight", l, row["index_ratio"]),
        ):
            moved = added(at_reference, c.scaled(part, row[key] - now))
            self.assertAlmostEqual(rho(k_inf(moved)), rho(reference_k(t)), delta=1e-6)

    def test_parts_are_the_effects_of_removing_them(self):
        t, reference, core = histories()
        (row,) = c.compare(t, reference, core)["rows"]
        self.assertAlmostEqual(
            row["without_untracked_pcm"] - row["constant_pcm"],
            -row["untracked_pcm"],
            delta=1e-9,
        )
        self.assertAlmostEqual(
            row["without_linearization_pcm"] - row["constant_pcm"],
            -row["linearization_pcm"],
            delta=1e-9,
        )

    def test_a_moderator_density_history_off_nominal_is_an_error(self):
        t, reference, core = histories()
        core["points"][2]["histories"] = [[[0.7, 0.5]]]
        with self.assertRaises(SystemExit):
            c.compare(t, reference, core)

    def test_table_items_outside_the_constants_are_not_interpolated(self):
        # 格子計算のテーブルの状態は、k の組み立てに使わない項目（値の無い vacuum_ratio、真偽の normalized_adf、None の pins）と、
        # 微視的断面積の反応の一覧を持つ
        t, reference, core = histories()
        for state in t["states"]:
            state |= {"vacuum_ratio": [None], "normalized_adf": [False], "pins": None}
        for nuclides in t["micro"]["states"]:
            for n in nuclides:
                n["reactions"] = [{"mt": 16, "channels": [[0.0, None]]}]
        (row,) = c.compare(t, reference, core)["rows"]
        self.assertAlmostEqual(row["reconstruct_pcm"], 0.0, delta=1e-9)


class Forgetting(unittest.TestCase):
    def test_the_infinite_time_constant_is_the_core(self):
        result = c.compare(*histories(), taus=[None, 5.0])
        (row,) = result["rows"]
        infinite = result["forgetting"][0]
        self.assertIsNone(infinite["tau"])
        (p,) = infinite["points"]
        self.assertAlmostEqual(p["weight"], row["weight"], delta=1e-15)
        self.assertAlmostEqual(p["constant_pcm"], row["constant_pcm"], delta=1e-9)

    def test_the_weight_is_one_less_the_exponential_weight_of_the_rest(self):
        # 制御棒を入れた区間が燃焼度 0〜B₀ なので、w_τ = 1 − (B₀ から後の重み)
        result = c.compare(*histories(), taus=[5.0, -5.0])
        for r in result["forgetting"]:
            (p,) = r["points"]
            self.assertAlmostEqual(
                p["weight"], 1.0 - cf.weight(20.0, SWITCH, r["tau"]), delta=1e-15
            )

    def test_the_best_time_constant_has_the_smallest_largest_constant_side(self):
        result = c.compare(*histories(), taus=[None, 5.0, -5.0, 50.0])
        best = min(result["forgetting"], key=lambda r: r["largest_pcm"])
        self.assertEqual(result["best"], best["tau"])


class Report(unittest.TestCase):
    def test_report_names_the_parts_and_the_best_time_constant(self):
        text = c.report(c.compare(*histories(), taus=[None, 5.0]))
        for word in ("数密度の効き", "定数の側", "τ"):
            self.assertIn(word, text)
        self.assertNotRegex(text, r"(?<![A-Za-z])-[0-9]")


if __name__ == "__main__":
    unittest.main()
