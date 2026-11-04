"""compare_rods_split.py（制御棒を途中まで入れて抜く燃焼の定数の側を、今の core の組成で U の外れと L の外れに分ける）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_rods_split.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_forgetting as cf
import compare_rods_split as c
from compare_residual import added, coefficients, corrected, k_inf, rho
from test_compare_rods_parts import (
    CORE_DENSITIES,
    NOMINAL,
    REFERENCE_DENSITIES,
    SWITCH,
    expected_node,
    histories,
)


def scaled(r, w):
    return c.scaled(r, w)


def tables(untracked_weight=0.5, linearization_weight=None):
    """公称のテーブル、参照、core と、参照の組成の分岐とハイブリッドのテーブル。切り替えの点（燃焼度 10）の参照の U と L はテーブルの
    値 out の U(1) と L(1)、燃焼度 20 の U_switch はその U(1) の `untracked_weight` 倍、L_switch は L(1) の `linearization_weight` 倍
    （None なら core の P / P_j 倍）"""
    t, reference, core = histories()
    moderator = t["axes"][2]
    moderator["density"] = [0.743, 0.661]
    out = t["histories"][1]["rods"][1]
    _, _, a = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
    lw = a if linearization_weight is None else linearization_weight
    parts = [
        (1, out["densities"][1], out["untracked"][1], out["linearization"][1]),
        (
            2,
            REFERENCE_DENSITIES[2],
            scaled(out["untracked"][2], untracked_weight),
            scaled(out["linearization"][2], lw),
        ),
    ]
    branch_states, hybrid_states, densities = [], [], []
    for b, d, u, lin in parts:
        p = t["points"].index([b, 0, 0, 0, 0])
        linear = corrected(
            t["states"][p],
            coefficients(t["micro"]["states"][p]),
            [x - y for x, y in zip(d, NOMINAL[b])],
        )
        hybrid = added(linear, lin)
        branch_states.append(added(hybrid, u))
        hybrid_states.append(hybrid)
        densities.append(d)

    def single(states):
        return {
            "axes": [
                {"name": "burnup", "points": [SWITCH, 20.0]},
                {"name": "boron", "points": [1300.0], "nominal": 1300.0},
                {
                    "name": "moderator",
                    "points": [565.0],
                    "density": [0.743],
                    "nominal": 565.0,
                },
                {"name": "fuel", "points": [565.0], "nominal": 565.0},
            ],
            "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
            "states": states,
            "micro": {
                "nuclides": t["micro"]["nuclides"],
                "densities": [list(x) for x in densities],
                "states": [t["micro"]["states"][0], t["micro"]["states"][0]],
            },
        }

    return t, reference, single(branch_states), single(hybrid_states), core


class Compare(unittest.TestCase):
    def test_rows_are_the_points_after_the_switch(self):
        result = c.compare(*tables())
        self.assertEqual(result["switch"], SWITCH)
        self.assertEqual([r["burnup"] for r in result["rows"]], [20.0])

    def test_the_switch_point_compares_the_branch_with_the_rods_out_residual(self):
        s = c.compare(*tables())["switch_point"]
        self.assertEqual(s["burnup"], SWITCH)
        for key in ("untracked", "linearization", "densities"):
            self.assertLess(s[key], 1e-15, key)

    def test_reference_parts_that_the_core_adds_leave_no_gap(self):
        # U_switch = w U(1)、L_switch = (P / P_j) L(1) なら外れは 0 で、定数の側は数密度の効きのほかの残りだけ
        (row,) = c.compare(*tables())["rows"]
        self.assertAlmostEqual(row["untracked_gap_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(row["linearization_gap_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(row["gap_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(row["remainder_pcm"], row["constant_pcm"], delta=1e-9)

    def test_gaps_are_the_effects_on_the_reference_constants(self):
        t, reference, branch, hybrid, core = tables(0.2, 0.1)
        (row,) = c.compare(t, reference, branch, hybrid, core)["rows"]
        out = t["histories"][1]["rods"][1]
        _, _, a = expected_node(t, 2, CORE_DENSITIES[2], 0.5)
        f = branch["states"][1]
        du = scaled(out["untracked"][2], 0.5 - 0.2)
        dl = scaled(out["linearization"][2], a - 0.1)
        for key, gap in (
            ("untracked_gap_pcm", du),
            ("linearization_gap_pcm", dl),
            ("gap_pcm", added(du, dl)),
        ):
            self.assertAlmostEqual(
                row[key], rho(k_inf(added(f, gap))) - rho(k_inf(f)), delta=1e-9
            )

    def test_reference_weights_are_ratios_of_the_effects(self):
        # 参照の U と L が U(1) と L(1) そのものなら、重みは 1
        (row,) = c.compare(*tables(1.0, 1.0))["rows"]
        self.assertAlmostEqual(row["untracked_weight"], 1.0, delta=1e-12)
        self.assertAlmostEqual(row["linearization_weight"], 1.0, delta=1e-12)

    def test_the_composition_effect_moves_the_reference_to_the_core_densities(self):
        t, reference, branch, hybrid, core = tables()
        (row,) = c.compare(t, reference, branch, hybrid, core)["rows"]
        f = branch["states"][1]
        p = t["points"].index([2, 0, 0, 0, 0])
        dense = corrected(
            f,
            coefficients(t["micro"]["states"][p]),
            [x - y for x, y in zip(CORE_DENSITIES[2], REFERENCE_DENSITIES[2])],
        )
        self.assertAlmostEqual(
            row["composition_pcm"], rho(k_inf(dense)) - rho(k_inf(f)), delta=1e-9
        )

    def test_a_hybrid_with_other_densities_is_an_error(self):
        t, reference, branch, hybrid, core = tables()
        hybrid["micro"]["densities"][1][0] *= 1.01
        with self.assertRaises(SystemExit):
            c.compare(t, reference, branch, hybrid, core)

    def test_a_branch_off_the_nominal_state_is_an_error(self):
        t, reference, branch, hybrid, core = tables()
        for table in (branch, hybrid):
            table["axes"][2]["density"] = [0.661]
        with self.assertRaises(SystemExit):
            c.compare(t, reference, branch, hybrid, core)

    def test_branch_densities_other_than_the_reference_are_errors(self):
        t, reference, branch, hybrid, core = tables()
        reference["points"][2]["densities"] = [CORE_DENSITIES[2]]
        with self.assertRaises(SystemExit):
            c.compare(t, reference, branch, hybrid, core)

    def test_a_core_k_other_than_its_constants_is_an_error(self):
        t, reference, branch, hybrid, core = tables()
        core["points"][2]["k"] *= 1.001
        with self.assertRaises(SystemExit):
            c.compare(t, reference, branch, hybrid, core)


class Forgetting(unittest.TestCase):
    def test_the_infinite_time_constant_is_the_core(self):
        result = c.compare(*tables(0.2, 0.1), taus=[None, 5.0])
        (row,) = result["rows"]
        (p,) = result["forgetting"][0]["points"]
        self.assertAlmostEqual(p["weight"], row["weight"], delta=1e-15)
        self.assertAlmostEqual(
            p["untracked_gap_pcm"], row["untracked_gap_pcm"], delta=1e-9
        )
        self.assertAlmostEqual(p["gap_pcm"], row["gap_pcm"], delta=1e-9)

    def test_the_weight_is_one_less_the_exponential_weight_of_the_rest(self):
        result = c.compare(*tables(), taus=[5.0, -5.0])
        for r in result["forgetting"]:
            (p,) = r["points"]
            self.assertAlmostEqual(
                p["weight"], 1.0 - cf.weight(20.0, SWITCH, r["tau"]), delta=1e-15
            )

    def test_the_best_time_constant_has_the_smallest_largest_gap(self):
        # 参照の U が U(1) の 0.1 倍なら、w_τ の小さい正の τ が良い
        result = c.compare(*tables(0.1), taus=[None, 1.0, 50.0, -50.0])
        best = min(result["forgetting"], key=lambda r: r["largest_pcm"])
        self.assertEqual(result["best"], best["tau"])
        self.assertEqual(result["best"], 1.0)
        self.assertTrue(result["best_at_edge"])


class Report(unittest.TestCase):
    def test_report_names_the_gaps_and_the_best_time_constant(self):
        text = c.report(c.compare(*tables(0.2, 0.1), taus=[None, 5.0]))
        for word in ("U の外れ", "L の外れ", "定数の側", "τ"):
            self.assertIn(word, text)
        self.assertNotRegex(text, r"(?<![A-Za-z])-[0-9]")


if __name__ == "__main__":
    unittest.main()
