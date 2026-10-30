"""compare_split.py（履歴の残り R を追わない核種の分と線形化の誤差に分ける）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_split.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_residual as r
import compare_split as c
import compare_weight as w
from test_compare_core_history import branch_table
from test_compare_residual import few
from test_compare_switch import reference
from test_compare_switch import table as nominal_table
from test_compare_weight import R_ONE, condition_tables


def hybrid(nominal, branch, linearizations):
    """`branch` と同じ条件、燃焼度、追う核種の数密度で、線形化の誤差が `linearizations`（燃焼度の点ごと）になる定数の 1 状態の
    テーブル。`nominal` は同じ条件の状態を持つ公称の組成のテーブル"""
    out = copy.deepcopy(branch)
    condition, burnups, states = w.condition_table(branch, "branch")
    at = c.solved(nominal, condition)
    for b, (u, p) in enumerate(zip(burnups, states)):
        bn, pn = at[u]
        deltas = [
            x - y
            for x, y in zip(
                branch["micro"]["densities"][b], nominal["micro"]["densities"][bn]
            )
        ]
        linear = r.corrected(
            nominal["states"][pn],
            r.coefficients(nominal["micro"]["states"][pn]),
            deltas,
        )
        state = r.added(linear, linearizations[b])
        state["k_inf"] = r.k_inf(state)
        out["states"][p] = state
    return out


def part(scale):
    """燃焼度の点に依らない、`scale` に比例する定数の差（効きが差に 1 次でほぼ比例する小ささ）"""
    return w.scaled(r.difference(few(1.03), few(1.0)), scale)


def switch_reference():
    """途中で条件を変えた燃焼の参照（燃焼度 0 と 20 MWd/kg を条件 j で）。R(1) と同じ燃焼度の点"""
    ref = reference()
    ref["axes"][0]["points"] = [0.0, 20.0]
    ref["micro"]["densities"][0] = [7.0e-4, 0.0]
    return ref


def summary(weights):
    """switch.json の行（参照の燃焼度 0 と 20 MWd/kg）"""
    return {
        "continued_burnup": 0.0,
        "rows": [
            {
                "burnup": u,
                "weight": x,
                "difference_pcm": d,
                "composition_pcm": m,
            }
            for u, x, d, m in zip([0.0, 20.0], weights, [-1.0, 30.0], [-1.0, -5.0])
        ],
    }


def density_tables(half_scale, half_linearization):
    """密度だけの経路の（公称の組成、分岐、ハイブリッド）を w = 0.5 と 1 で。R(0.5) = `half_scale` R(1)、L(1) = part(1)、
    L(0.5) = `half_linearization` L(1)"""
    half = condition_tables(
        0.702,
        [[7.0e-4, 0.0], [5.1e-4, 2.3e-5]],
        [w.scaled(x, half_scale) for x in R_ONE],
        0.5,
    )
    one = condition_tables(0.661, [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]], R_ONE, 1.0)
    return (
        (*half, hybrid(*half, [part(0.0), part(half_linearization)])),
        (*one, hybrid(*one, [part(0.0), part(1.0)])),
    )


def tables(half_scale=0.5, half_linearization=0.5, switch_linearization=None):
    nominal = nominal_table()
    one = branch_table()
    ref = switch_reference()
    half, full = density_tables(half_scale, half_linearization)
    # R(1) の線形化の誤差は燃焼度 0 で 0（ハイブリッドが分岐と同じ定数）、20 MWd/kg で part(1)
    one_linearization = [part(0.0), part(1.0)]
    if switch_linearization is None:
        switch_linearization = [part(0.2), part(0.4)]
    return {
        "nominal": nominal,
        "one": one,
        "one_hybrid": hybrid(nominal, one, one_linearization),
        "switch": ref,
        "switch_hybrid": hybrid(nominal, ref, switch_linearization),
        "nominal_half": half[0],
        "half": half[1],
        "half_hybrid": half[2],
        "nominal_full": full[0],
        "full": full[1],
        "full_hybrid": full[2],
    }


def rho(state):
    return r.rho(r.k_inf(state))


class Split(unittest.TestCase):
    def test_untracked_and_linearization_parts_add_up_to_the_residual(self):
        result = c.compare(tables(), summary([0.0, 0.5]))
        for x in result["one"]:
            self.assertAlmostEqual(
                x["untracked_condition_pcm"] + x["linearization_condition_pcm"],
                x["residual_condition_pcm"],
                delta=1e-9,
            )
            self.assertLess(x["split_check"], 1e-15)

    def test_the_condition_residual_is_the_branch_less_the_linear_correction(self):
        t = tables()
        last = c.compare(t, summary([0.0, 0.5]))["one"][-1]
        # 燃焼度 20 の条件 j の状態は点 4、N̄ᵢ⁰ は燃焼度の点 2
        nominal = t["nominal"]
        deltas = [
            a - b
            for a, b in zip(
                t["one"]["micro"]["densities"][1], nominal["micro"]["densities"][2]
            )
        ]
        linear = r.corrected(
            nominal["states"][4], r.coefficients(nominal["micro"]["states"][4]), deltas
        )
        self.assertAlmostEqual(
            last["residual_condition_pcm"],
            rho(t["one"]["states"][1]) - rho(linear),
            delta=1e-9,
        )
        self.assertAlmostEqual(
            last["linearization_condition_pcm"],
            rho(r.added(linear, part(1.0))) - rho(linear),
            delta=1e-9,
        )

    def test_a_hybrid_at_the_linear_correction_leaves_the_residual_untracked(
        self,
    ):
        first = c.compare(tables(), summary([0.0, 0.5]))["one"][0]
        self.assertAlmostEqual(first["linearization_pcm"], 0.0, delta=1e-12)
        self.assertAlmostEqual(
            first["untracked_pcm"], first["residual_pcm"], delta=1e-9
        )

    def test_effects_are_taken_on_the_nominal_state_with_the_history_densities(self):
        t = tables()
        last = c.compare(t, summary([0.0, 0.5]))["one"][-1]
        nominal = t["nominal"]
        # 燃焼度 20 の公称の状態は点 3
        deltas = [
            a - b
            for a, b in zip(
                t["one"]["micro"]["densities"][1], nominal["micro"]["densities"][2]
            )
        ]
        base = r.corrected(
            nominal["states"][3], r.coefficients(nominal["micro"]["states"][3]), deltas
        )
        self.assertAlmostEqual(
            last["linearization_pcm"],
            rho(r.added(base, part(1.0))) - rho(base),
            delta=1e-9,
        )

    def test_hybrids_with_other_tracked_densities_are_errors(self):
        t = tables()
        t["one_hybrid"]["micro"]["densities"][1][0] *= 1.001
        with self.assertRaises(SystemExit):
            c.compare(t, summary([0.0, 0.5]))

    def test_hybrids_of_other_burnups_are_errors(self):
        t = tables()
        t["switch_hybrid"]["axes"][0]["points"] = [10.0, 21.0]
        with self.assertRaises(SystemExit):
            c.compare(t, summary([0.0, 0.5]))

    def test_density_paths_need_the_nominal_composition_of_the_nominal_table(self):
        t = tables()
        # 数密度の列は test_compare_weight の N0 と共有なので、コピーしてから変える
        t["nominal_half"] = copy.deepcopy(t["nominal_half"])
        t["nominal_half"]["micro"]["densities"][1][1] *= 1.001
        with self.assertRaises(SystemExit):
            c.compare(t, summary([0.0, 0.5]))


class Switch(unittest.TestCase):
    def test_parts_proportional_to_w_explain_none_of_the_core_difference(self):
        # R(1) の 20 MWd/kg の線形化の誤差は part(1)。switch の L = w part(1)（w = 0.5）
        t = tables(switch_linearization=[part(0.0), part(0.5)])
        last = c.compare(t, summary([0.0, 0.5]))["switch"][-1]
        self.assertAlmostEqual(last["linearization_gap_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(
            last["linearization_pcm"], last["scaled_linearization_pcm"], delta=1e-9
        )

    def test_gaps_are_the_reference_with_the_part_replaced_by_its_share_of_r_one(
        self,
    ):
        t = tables()
        result = c.compare(t, summary([0.0, 0.5]))
        last = result["switch"][-1]
        # switch の L は part(0.4)、R(1) の L は part(1)、w = 0.5: 外れは part(0.4) − 0.5 part(1) = part(−0.1)
        ref = t["switch"]["states"][1]
        self.assertAlmostEqual(
            last["linearization_gap_pcm"],
            rho(r.added(ref, part(0.1))) - rho(ref),
            delta=1e-9,
        )
        self.assertAlmostEqual(
            last["remainder_pcm"],
            last["difference_pcm"] - last["composition_pcm"] - last["gap_pcm"],
            delta=1e-12,
        )

    def test_weights_of_the_parts_are_the_ratios_of_their_effects(self):
        last = c.compare(tables(), summary([0.0, 0.5]))["switch"][-1]
        # L の効きは part に 1 次でほぼ比例: switch の part(0.4) と R(1) の part(1) の比
        self.assertAlmostEqual(last["linearization_weight"], 0.4, delta=0.002)
        self.assertAlmostEqual(
            last["linearization_weight"],
            last["linearization_pcm"] / last["one_linearization_pcm"],
            delta=1e-12,
        )

    def test_summaries_of_other_burnups_are_errors(self):
        s = summary([0.0, 0.5])
        s["rows"][1]["burnup"] = 22.0
        with self.assertRaises(SystemExit):
            c.compare(tables(), s)


class Density(unittest.TestCase):
    def test_parts_proportional_to_w_give_no_difference_from_half_of_r_one(self):
        rows = c.compare(tables(0.5, 0.5), summary([0.0, 0.5]))["density"]
        last = rows[-1]
        for key in ("untracked", "linearization", "residual"):
            self.assertAlmostEqual(last[f"{key}_difference_pcm"], 0.0, delta=1e-9)
            self.assertAlmostEqual(last[f"{key}_ratio"], 0.5, delta=1e-6)
        # 燃焼度 0 は比を出さない
        self.assertIsNone(rows[0]["linearization_ratio"])

    def test_the_nonlinear_part_shows_its_own_ratio(self):
        last = c.compare(tables(0.5, 0.6), summary([0.0, 0.5]))["density"][-1]
        # L(0.5) = 0.6 L(1)、R は比例: L の比は 0.6、U はその分だけ逆に外れる
        self.assertAlmostEqual(last["linearization_ratio"], 0.6, delta=0.002)
        self.assertAlmostEqual(last["residual_difference_pcm"], 0.0, delta=1e-9)
        self.assertGreater(
            -last["untracked_difference_pcm"] * last["linearization_difference_pcm"],
            0.0,
        )


class Report(unittest.TestCase):
    def test_report_names_the_parts_at_the_last_points(self):
        result = c.compare(tables(), summary([0.0, 0.5]))
        text = c.report(result)
        self.assertIn("追わない核種の分", text)
        self.assertIn("線形化の誤差", text)


if __name__ == "__main__":
    unittest.main()
