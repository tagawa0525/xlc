"""compare_weight.py（履歴の重み ω の core の差と、密度だけの経路の R(0.5) と R(1)/2 の効き）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_weight.py'
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_residual as r
import compare_weight as c
from test_compare_core_history import branch_table, core, sigma
from test_compare_core_history import table as residual_table
from test_compare_residual import NUCLIDES, axes, few, residual
from test_compare_switch import reference

REACTIONS = [("U-235", "fission"), ("Pu-240", 102)]
# 公称のテーブルの燃焼度 0 と 20 の N̄ᵢ⁰（test_compare_core_history の table）
N0 = [[7.0e-4, 0.0], [5.0e-4, 2.0e-5]]
# 密度だけの経路の R(1)（燃焼度 0 と 20）
R_ONE = [residual(few(1.2), few(1.0)), residual(few(4.9), few(4.5))]


def condition_tables(density, densities, residuals, sigma_shift):
    """565 K、`density` g/cc の 1 状態の、公称の組成のテーブルと、R が `residuals` になる分岐の組成のテーブル（燃焼度 0 と 20）"""
    nominal_sigma = [
        [sigma(3, False), sigma(4, True)],
        [sigma(5, False), sigma(6, True)],
    ]
    nominal = {
        "format": "xlc-table",
        "axes": axes([0.0, 20.0], [565.0], [density], 565.0),
        "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
        "states": [few(1.0), few(4.5)],
        "micro": {"nuclides": NUCLIDES, "densities": N0, "states": nominal_sigma},
    }
    states = []
    for k, rk in enumerate(residuals):
        deltas = [a - b for a, b in zip(densities[k], N0[k])]
        linear = r.corrected(
            nominal["states"][k], r.coefficients(nominal_sigma[k]), deltas
        )
        states.append(r.added(linear, rk))
    branch = {
        "format": "xlc-table",
        "axes": axes([0.0, 20.0], [565.0], [density], 565.0),
        "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
        "states": states,
        "micro": {
            "nuclides": NUCLIDES,
            "densities": densities,
            "states": [
                [sigma(3 + sigma_shift, False), sigma(4 + sigma_shift, False)],
                [sigma(5 + sigma_shift, False), sigma(6 + sigma_shift, False)],
            ],
        },
    }
    return nominal, branch


def residual_summary(scale):
    """R(0.5) = `scale` R(1) の、密度だけの経路の比べた要約"""
    half = condition_tables(
        0.702,
        [[7.0e-4, 0.0], [5.1e-4, 2.3e-5]],
        [c.scaled(x, scale) for x in R_ONE],
        0.5,
    )
    one = condition_tables(0.661, [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]], R_ONE, 1.0)
    with mock.patch.object(c, "REACTIONS", REACTIONS):
        return c.compare_residuals(
            residual_table(), residual_table(), branch_table(), half, one
        )


def switch(estimates):
    """途中で条件を変えた燃焼の要約（10 MWd/kg から続け、参照は 10 と 20 MWd/kg）"""
    return {
        "continued_burnup": 10.0,
        "rows": [
            {
                "burnup": 10.0,
                "weight": 0.0,
                "needed_weight": None,
                "difference_pcm": -1.0,
            },
            {
                "burnup": 20.0,
                "weight": 0.5,
                "needed_weight": 0.6,
                "difference_pcm": 30.0,
            },
        ],
        "omegas": [{"omega": 2.0, "weights": [0.0, 0.6], "estimate_pcm": estimates}],
        "inputs": {
            "table-residual-nominal-histories.msgpack": "aaa",
            "table-switch-density-branch.msgpack": "bbb",
            "core-2A-residual-switch.msgpack": "ddd",
        },
    }


class Residuals(unittest.TestCase):
    def test_proportional_residuals_give_no_difference_from_half_of_r_one(self):
        rows = residual_summary(0.5)["rows"]
        for x in rows:
            self.assertAlmostEqual(x["difference_pcm"], 0.0, delta=1e-9)
            self.assertAlmostEqual(x["half_pcm"], x["one_half_pcm"], delta=1e-9)
        # Δσ も比例なら比は 0.5（燃焼度 0 の Pu-240 のように Δσ(1) が 0 の群は None）
        for name, values in rows[1]["sigma_ratios"].items():
            for v in values:
                self.assertAlmostEqual(v, 0.5, places=12, msg=name)

    def test_residuals_larger_than_proportional_add_their_excess_effect(self):
        rows = residual_summary(0.6)["rows"]
        last = rows[-1]
        # R(0.5) = 0.6 R(1) の効きは、R(1)/2 の効きより同じ符号に大きい（効きは R に 1 次でほぼ比例）
        self.assertNotAlmostEqual(last["one_half_pcm"], 0.0, places=6)
        self.assertGreater(last["difference_pcm"] * last["one_half_pcm"], 0.0)
        self.assertAlmostEqual(last["half_pcm"] / last["one_half_pcm"], 1.2, delta=0.01)

    def test_residual_effects_are_taken_on_the_nominal_state_with_the_half_densities(
        self,
    ):
        rows = residual_summary(0.6)["rows"]
        t = residual_table()
        # 燃焼度 20 の公称の状態は点 3、N̄ᵢ⁰ は燃焼度の点 2
        f0 = t["states"][3]
        deltas = [a - b for a, b in zip([5.1e-4, 2.3e-5], t["micro"]["densities"][2])]
        base = r.corrected(f0, r.coefficients(t["micro"]["states"][3]), deltas)
        expected = r.rho(r.k_inf(r.added(base, c.scaled(R_ONE[1], 0.6)))) - r.rho(
            r.k_inf(base)
        )
        self.assertAlmostEqual(rows[1]["half_pcm"], expected, delta=1e-9)
        # 同じ R(1) を w = 0.5 の定数に足した効き（R(1)/2 の効きとの違いは R の大きさだけ、R(1) の効きとの違いは足す定数だけ）
        whole = r.rho(r.k_inf(r.added(base, R_ONE[1]))) - r.rho(r.k_inf(base))
        self.assertAlmostEqual(rows[1]["one_on_half_pcm"], whole, delta=1e-9)

    def test_the_temperature_path_residual_is_rebuilt_from_the_tables(self):
        rows = residual_summary(0.5)["rows"]
        for x in rows:
            self.assertLess(x["temperature_path_rebuilt"], 1e-15)

    def test_conditions_that_change_more_than_the_density_are_errors(self):
        half = condition_tables(0.702, [[7.0e-4, 0.0], [5.1e-4, 2.3e-5]], R_ONE, 0.5)
        one = condition_tables(0.661, [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]], R_ONE, 1.0)
        for table in one:
            table["axes"][2]["points"] = [600.0]
        with (
            self.assertRaises(SystemExit),
            mock.patch.object(c, "REACTIONS", REACTIONS),
        ):
            c.compare_residuals(
                residual_table(), residual_table(), branch_table(), half, one
            )

    def test_nominal_compositions_must_match_the_nominal_table(self):
        half = condition_tables(0.702, [[7.0e-4, 0.0], [5.1e-4, 2.3e-5]], R_ONE, 0.5)
        one = condition_tables(0.661, [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]], R_ONE, 1.0)
        half[0]["micro"]["densities"] = [[7.0e-4, 0.0], [5.0e-4, 2.1e-5]]
        with (
            self.assertRaises(SystemExit),
            mock.patch.object(c, "REACTIONS", REACTIONS),
        ):
            c.compare_residuals(
                residual_table(), residual_table(), branch_table(), half, one
            )


# 途中で条件を変えた燃焼の core の、R の点の w と w R の効き（core-history.json の switch_rows）
CORE_SWITCH = [
    {"burnup": 10.0, "weight": 0.0, "residual_pcm": 0.0},
    {"burnup": 20.0 + 1e-12, "weight": 0.5, "residual_pcm": -120.0},
]


class Quadratic(unittest.TestCase):
    def test_the_quadratic_changes_the_difference_by_the_excess_of_r_half(self):
        rows = residual_summary(0.6)["rows"]
        (q,) = c.quadratic(switch([None, 5.0]), rows, CORE_SWITCH)
        x = rows[1]
        ratio = x["half_pcm"] / (2.0 * x["one_half_pcm"])
        # core が足した w R(1) の効きを、w R(1) + 4 w (1 − w) (R(0.5) − R(1)/2) に置き換えた変化:
        # 4 (1 − w) (r − 1/2) × (w R(1) の効き)
        change = 4.0 * 0.5 * (ratio - 0.5) * -120.0
        self.assertAlmostEqual(q["ratio"], ratio, delta=1e-15)
        self.assertAlmostEqual(q["change_pcm"], change, delta=1e-12)
        self.assertAlmostEqual(q["estimate_pcm"], 30.0 + change, delta=1e-12)
        self.assertEqual(q["needed_weight"], 0.6)

    def test_proportional_residuals_leave_the_difference(self):
        (q,) = c.quadratic(
            switch([None, 5.0]), residual_summary(0.5)["rows"], CORE_SWITCH
        )
        self.assertAlmostEqual(q["change_pcm"], 0.0, delta=1e-9)

    def test_core_runs_of_other_weights_are_errors(self):
        other = [dict(x) for x in CORE_SWITCH]
        other[1]["weight"] = 0.4
        with self.assertRaises(SystemExit):
            c.quadratic(switch([None, 5.0]), residual_summary(0.6)["rows"], other)


class Omegas(unittest.TestCase):
    def cores(self):
        h20 = 0.743 + 0.6 * (0.661 - 0.743)
        ref = reference()
        k20 = ref["states"][1]["k_inf"] * 1.0002
        nodes = core(
            [
                (0.0, 0.0, "nominal", 1.2, [7.0e-4, 0.0], 0.743),
                (
                    5.0,
                    10.0,
                    "branch",
                    ref["states"][0]["k_inf"],
                    [6.0e-4, 1.0e-5],
                    0.743,
                ),
                (10.0, 20.0, "branch", k20, [5.25e-4, 2.5e-5], h20),
            ]
        )
        return ref, k20, {2.0: nodes}

    def test_measured_differences_are_set_beside_the_estimates(self):
        ref, k20, cores = self.cores()
        (o,) = c.compare_omegas(residual_table(), ref, switch([None, 5.0]), cores)
        self.assertEqual(o["omega"], 2.0)
        first, last = o["rows"]
        self.assertAlmostEqual(first["difference_pcm"], 0.0, delta=1e-9)
        self.assertIsNone(first["estimate_pcm"])
        self.assertAlmostEqual(last["weight"], 0.6, places=12)
        self.assertAlmostEqual(
            last["difference_pcm"],
            r.rho(k20) - r.rho(ref["states"][1]["k_inf"]),
            delta=1e-9,
        )
        self.assertEqual(last["estimate_pcm"], 5.0)

    def test_weights_without_an_estimate_are_errors(self):
        ref, _, cores = self.cores()
        with self.assertRaises(SystemExit):
            c.compare_omegas(
                residual_table(), ref, switch([None, 5.0]), {1.5: cores[2.0]}
            )

    def test_switch_summaries_of_other_tables_are_errors(self):
        md5s = {
            "table-residual-nominal-histories.msgpack": "aaa",
            "table-switch-density-branch.msgpack": "ccc",
            "core-2A-residual-switch.msgpack": "ddd",
        }
        with self.assertRaises(SystemExit):
            c.check_switch(switch([None, 5.0]), md5s)
        md5s["table-switch-density-branch.msgpack"] = "bbb"
        c.check_switch(switch([None, 5.0]), md5s)
        # core-history.json も同じ core の途中で条件を変えた燃焼から作ったもの
        with self.assertRaises(SystemExit):
            c.check_switch(
                switch([None, 5.0]),
                md5s,
                {"inputs": {"core-2A-residual-switch.msgpack": "eee"}},
            )
        c.check_switch(
            switch([None, 5.0]),
            md5s,
            {"inputs": {"core-2A-residual-switch.msgpack": "ddd"}},
        )


if __name__ == "__main__":
    unittest.main()
