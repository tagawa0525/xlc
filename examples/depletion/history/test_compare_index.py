"""compare_index.py（線形化の誤差と追わない核種の分を、組成の指標で表せるか）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_index.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_index as c
import compare_residual as r
import compare_split as s
from test_compare_residual import NUCLIDES, few
from test_compare_split import hybrid, summary, tables

WEIGHTS = [0.0, 0.5]


def rho(state):
    return r.rho(r.k_inf(state))


def split_of(t, weights=WEIGHTS):
    """compare_split.py が同じテーブルから書く split.json"""
    return s.compare(t, summary(weights))


def compared(t=None, weights=WEIGHTS):
    t = tables() if t is None else t
    return c.compare(t, summary(weights), split_of(t, weights), NUCLIDES)


def index_of(t, table, burnup):
    """`table` の燃焼度 `burnup`（テーブルの燃焼度の点の番号）の数密度の、公称の状態での線形の補正の効き"""
    nominal = t["nominal"]
    densities = t[table]["micro"]["densities"][burnup]
    u = t[table]["axes"][0]["points"][burnup]
    _, p = s.nominal_state(nominal, u)
    return rho(s.base(nominal, u, densities)) - rho(nominal["states"][p])


# switch の 20 MWd/kg の数密度（公称と R(1) の間。指標の比が 1 でない）
SWITCH20 = [6.0e-4, 1.5e-5]


def thermal(scale):
    """熱群の吸収だけを `scale` に比例して変える定数の差（`part` は定数を一様に拡大するので、拡大した定数の k を変えない）"""
    out = r.difference(few(1.0), few(1.0))
    out["total"][1] = 0.01 * scale
    out["absorption"][1] = 0.01 * scale
    return out


def fixture(switch_linearization=0.4):
    """R(1) の L が thermal(1)、switch の L が thermal(`switch_linearization`) のテーブル（20 MWd/kg）"""
    t = tables()
    t["switch"]["micro"]["densities"][1] = list(SWITCH20)
    t["one_hybrid"] = hybrid(t["nominal"], t["one"], [thermal(0.0), thermal(1.0)])
    t["switch_hybrid"] = hybrid(
        t["nominal"], t["switch"], [thermal(0.0), thermal(switch_linearization)]
    )
    return t


def proportional_tables():
    """switch の L が R(1) の L の、指標の比 E_switch / E(1) 倍のテーブル（20 MWd/kg）"""
    ratio = index_of(fixture(), "switch", 1) / index_of(fixture(), "one", 1)
    return fixture(ratio), ratio


class Index(unittest.TestCase):
    def test_the_index_is_the_effect_of_the_linear_correction_on_the_nominal_state(
        self,
    ):
        t = tables()
        result = compared(t)
        for name, table in (("one", "one"), ("switch", "switch")):
            last = result["histories"][name][-1]
            self.assertAlmostEqual(last["index_pcm"], index_of(t, table, 1), delta=1e-9)

    def test_the_parts_are_the_effects_that_the_split_gives(self):
        t = tables()
        result = compared(t)
        split = split_of(t)
        for x, y in zip(result["histories"]["one"], split["one"]):
            self.assertEqual(x["linearization_pcm"], y["linearization_pcm"])
            self.assertEqual(x["untracked_pcm"], y["untracked_pcm"])
        for x, y in zip(result["histories"]["switch"], split["switch"]):
            self.assertEqual(x["linearization_pcm"], y["linearization_pcm"])
        # 密度だけの経路の w = 0.5 は、その履歴の数密度の線形の補正を足した定数（split の half_*）
        for x, y in zip(result["histories"]["half"], split["density"]):
            self.assertAlmostEqual(
                x["linearization_pcm"], y["half_linearization_pcm"], delta=1e-12
            )

    def test_ratios_are_of_the_part_to_the_index_and_its_square(self):
        last = compared()["histories"]["one"][-1]
        e = last["index_pcm"]
        self.assertAlmostEqual(
            last["linearization_first"], last["linearization_pcm"] / e, delta=1e-15
        )
        self.assertAlmostEqual(
            last["linearization_second"],
            last["linearization_pcm"] / e**2,
            delta=1e-18,
        )
        self.assertAlmostEqual(
            last["untracked_first"], last["untracked_pcm"] / e, delta=1e-15
        )
        # R(1) は w = 1
        self.assertEqual(last["weight"], 1.0)
        self.assertAlmostEqual(
            last["untracked_weight"], last["untracked_pcm"], delta=1e-15
        )

    def test_ratios_to_r_one_are_one_for_r_one(self):
        last = compared()["histories"]["one"][-1]
        for key in (
            "linearization_first",
            "linearization_second",
            "untracked_first",
            "untracked_weight",
        ):
            self.assertAlmostEqual(last[f"{key}_to_one"], 1.0, delta=1e-12)

    def test_points_at_the_nominal_composition_have_no_ratios(self):
        # R(1) の燃焼度 0 は公称の組成: 指標は 0 で、比を出さない
        first = compared()["histories"]["one"][0]
        self.assertEqual(first["index_pcm"], 0.0)
        self.assertIsNone(first["linearization_first"])
        self.assertIsNone(first["linearization_second"])
        self.assertIsNone(first["linearization_first_to_one"])
        # switch の w = 0 の点は U/w を出さない
        self.assertIsNone(compared()["histories"]["switch"][0]["untracked_weight"])

    def test_density_paths_carry_their_weights(self):
        histories = compared()["histories"]
        self.assertEqual({x["weight"] for x in histories["half"]}, {0.5})
        self.assertEqual({x["weight"] for x in histories["full"]}, {1.0})
        self.assertEqual([x["weight"] for x in histories["switch"]], WEIGHTS)

    def test_relative_densities_are_against_the_nominal_composition(self):
        t = tables()
        last = compared(t)["histories"]["one"][-1]
        # 燃焼度 20 の公称の組成は燃焼度の点 2
        n0 = t["nominal"]["micro"]["densities"][2]
        nx = t["one"]["micro"]["densities"][1]
        for i, name in enumerate(NUCLIDES):
            self.assertAlmostEqual(
                last["densities"][name], (nx[i] - n0[i]) / n0[i], delta=1e-12
            )

    def test_nuclides_absent_from_the_nominal_composition_have_no_relative_density(
        self,
    ):
        t = tables()
        # 燃焼度 0 の公称の組成に Pu-240 は無い
        self.assertEqual(t["nominal"]["micro"]["densities"][0][1], 0.0)
        first = compared(t)["histories"]["one"][0]
        self.assertIsNone(first["densities"][NUCLIDES[1]])

    def test_untracked_nuclides_are_errors(self):
        t = tables()
        with self.assertRaises(SystemExit):
            c.compare(t, summary(WEIGHTS), split_of(t), ["Xe-135"])


class Switch(unittest.TestCase):
    def test_a_part_proportional_to_the_index_leaves_no_first_order_gap(self):
        t, ratio = proportional_tables()
        last = compared(t)["switch"][-1]
        self.assertAlmostEqual(last["index_ratio"], ratio, delta=1e-15)
        self.assertAlmostEqual(last["linearization_gap_index_pcm"], 0.0, delta=1e-9)
        # w（0.5）で内挿した外れは残る
        self.assertGreater(abs(last["linearization_gap_weight_pcm"]), 1e-3)

    def test_gaps_are_the_reference_with_the_part_replaced_by_its_share_of_r_one(
        self,
    ):
        t = fixture()
        last = compared(t)["switch"][-1]
        ratio = last["index_ratio"]
        ref = t["switch"]["states"][1]
        # switch の L は thermal(0.4)、R(1) の L は thermal(1)
        for key, factor in (
            ("index", ratio),
            ("index_squared", ratio**2),
            ("weight", 0.5),
        ):
            expected = rho(r.added(ref, thermal(factor - 0.4))) - rho(ref)
            self.assertGreater(abs(expected), 1.0)
            self.assertAlmostEqual(
                last[f"linearization_gap_{key}_pcm"], expected, delta=1e-9
            )

    def test_gaps_at_w_are_those_of_the_split(self):
        t = fixture()
        result = compared(t)
        split = split_of(t)
        for x, y in zip(result["switch"], split["switch"]):
            for key in ("untracked", "linearization"):
                self.assertAlmostEqual(
                    x[f"{key}_gap_weight_pcm"], y[f"{key}_gap_pcm"], delta=1e-12
                )

    def test_splits_of_other_tables_are_errors(self):
        t = fixture()
        split = split_of(t)
        split["switch"][-1]["linearization_gap_pcm"] += 1e-6
        with self.assertRaises(SystemExit):
            c.compare(t, summary(WEIGHTS), split, NUCLIDES)

    def test_both_gaps_take_the_untracked_part_at_w_and_the_linearization_at_the_index(
        self,
    ):
        t, _ = proportional_tables()
        last = compared(t)["switch"][-1]
        # L は指標に比例するので、両方の外れは U の w の外れだけ
        self.assertAlmostEqual(
            last["gap_mixed_pcm"], last["untracked_gap_weight_pcm"], delta=1e-9
        )

    def test_the_condition_index_is_the_linear_correction_at_the_condition(self):
        t = fixture()
        last = compared(t)["switch"][-1]
        nominal = t["nominal"]
        condition, _ = s.split(nominal, t["one"], t["one_hybrid"], "R(1)")
        _, p = s.solved(nominal, condition)[20.0]
        own = s.split(nominal, t["switch"], t["switch_hybrid"], "switch")[1][-1]
        one = s.split(nominal, t["one"], t["one_hybrid"], "R(1)")[1][-1]
        state = nominal["states"][p]
        self.assertAlmostEqual(
            last["condition_index_ratio"],
            (rho(own["linear"]) - rho(state)) / (rho(one["linear"]) - rho(state)),
            delta=1e-12,
        )

    def test_density_ratios_are_of_the_switch_to_r_one(self):
        t = fixture()
        last = compared(t)["switch"][-1]
        n0 = t["nominal"]["micro"]["densities"][2]
        ns = t["switch"]["micro"]["densities"][1]
        n1 = t["one"]["micro"]["densities"][1]
        for i, name in enumerate(NUCLIDES):
            expected = None if n1[i] == n0[i] else (ns[i] - n0[i]) / (n1[i] - n0[i])
            if expected is None:
                self.assertIsNone(last["density_ratios"][name])
            else:
                self.assertAlmostEqual(
                    last["density_ratios"][name], expected, delta=1e-12
                )

    def test_weights_of_the_parts_are_the_split_ratios(self):
        t = fixture()
        last = compared(t)["switch"][-1]
        split = split_of(t)["switch"][-1]
        self.assertEqual(last["linearization_weight"], split["linearization_weight"])
        self.assertEqual(last["untracked_weight"], split["untracked_weight"])


def group_index(t, table, burnup, key):
    """`table` の燃焼度 `burnup`（燃焼度の点の番号）の数密度の、公称の状態の係数での熱群（最後の群）の線形の補正（`key` は
    `nu_fission` か `absorption`）"""
    nominal = t["nominal"]
    densities = t[table]["micro"]["densities"][burnup]
    u = t[table]["axes"][0]["points"][burnup]
    b, p = s.nominal_state(nominal, u)
    deltas = [x - y for x, y in zip(densities, nominal["micro"]["densities"][b])]
    coefficients = r.coefficients(nominal["micro"]["states"][p])
    if key == "absorption":
        return sum(
            (x["total"][-1] - sum(x["scatter"][-1])) * d
            for x, d in zip(coefficients, deltas)
        )
    return sum(x[key][-1] * d for x, d in zip(coefficients, deltas))


def production_tables():
    """switch の L が R(1) の L の、熱群の νΣf の補正の比 P_switch / P(1) 倍のテーブル（20 MWd/kg）"""
    t = fixture()
    ratio = group_index(t, "switch", 1, "nu_fission") / group_index(
        t, "one", 1, "nu_fission"
    )
    return fixture(ratio), ratio


class GroupIndices(unittest.TestCase):
    def test_group_indices_are_the_thermal_linear_corrections_of_the_nominal_state(
        self,
    ):
        t = fixture()
        result = compared(t)
        for name in ("one", "switch"):
            last = result["histories"][name][-1]
            self.assertAlmostEqual(
                last["production"], group_index(t, name, 1, "nu_fission"), delta=1e-15
            )
            self.assertAlmostEqual(
                last["thermal_absorption"],
                group_index(t, name, 1, "absorption"),
                delta=1e-15,
            )
        # 公称の組成の点は 0
        self.assertEqual(result["histories"]["one"][0]["production"], 0.0)

    def test_a_part_proportional_to_the_production_index_leaves_no_gap(self):
        t, ratio = production_tables()
        last = compared(t)["switch"][-1]
        self.assertAlmostEqual(last["production_ratio"], ratio, delta=1e-15)
        self.assertAlmostEqual(
            last["linearization_gap_production_pcm"], 0.0, delta=1e-9
        )
        self.assertGreater(abs(last["linearization_gap_index_pcm"]), 1e-3)
        # 両方の外れ（U は w、L は熱群の νΣf の補正の比）は U の w の外れだけ
        self.assertAlmostEqual(
            last["gap_production_pcm"], last["untracked_gap_weight_pcm"], delta=1e-9
        )

    def test_the_absorption_index_gap_takes_the_ratio_of_the_thermal_absorption(self):
        t = fixture()
        last = compared(t)["switch"][-1]
        ratio = group_index(t, "switch", 1, "absorption") / group_index(
            t, "one", 1, "absorption"
        )
        self.assertAlmostEqual(last["absorption_ratio"], ratio, delta=1e-15)
        ref = t["switch"]["states"][1]
        expected = rho(r.added(ref, thermal(ratio - 0.4))) - rho(ref)
        self.assertAlmostEqual(
            last["linearization_gap_absorption_pcm"], expected, delta=1e-9
        )

    def test_density_pairs_set_the_index_ratios_against_the_ratio_of_the_parts(self):
        result = compared()
        pairs = result["density_pairs"]
        half, full = result["histories"]["half"], result["histories"]["full"]
        self.assertEqual([x["burnup"] for x in pairs], [x["burnup"] for x in half])
        last, h, f = pairs[-1], half[-1], full[-1]
        for key, name in (
            ("index_pcm", "index_ratio"),
            ("production", "production_ratio"),
            ("thermal_absorption", "absorption_ratio"),
            ("linearization_pcm", "linearization_ratio"),
        ):
            self.assertAlmostEqual(last[name], h[key] / f[key], delta=1e-15)
        # 公称の組成の点は比を出さない
        self.assertIsNone(pairs[0]["production_ratio"])


class Report(unittest.TestCase):
    def test_report_names_the_group_indices(self):
        text = c.report(compared(fixture()))
        self.assertIn("熱群の νΣf の補正の比", text)
        self.assertIn("熱群の吸収の補正の比", text)

    def test_report_names_the_index_and_the_gaps(self):
        text = c.report(compared())
        self.assertIn("L/E", text)
        self.assertIn("指標の比", text)
        self.assertNotRegex(text, r"(?<![A-Za-z])-[0-9]")
        self.assertRegex(text, r"−[0-9]")

    def test_report_sets_the_parts_of_the_density_paths_against_their_index(self):
        result = compared()
        text = c.report(result)
        half = result["histories"]["half"][-1]
        full = result["histories"]["full"][-1]
        ratio = half["linearization_pcm"] / full["linearization_pcm"]
        self.assertIn("密度だけの経路", text)
        self.assertIn(f"{ratio:.4f}".replace("-", "−"), text)

    def test_report_keeps_the_hyphens_of_the_nuclide_names(self):
        text = c.report(compared(fixture()))
        self.assertIn("U-235", text)
        self.assertNotIn("U−235", text)

    def test_summaries_of_other_tables_are_errors(self):
        md5s = {name: "aaa" for name in c.INPUTS}
        made = {"inputs": dict(md5s)}
        c.check_summaries({"split.json": made, "switch.json": made}, md5s)
        changed = copy.deepcopy(made)
        changed["inputs"][s.FILES["switch"]] = "bbb"
        with self.assertRaises(SystemExit):
            c.check_summaries({"split.json": changed, "switch.json": made}, md5s)


if __name__ == "__main__":
    unittest.main()
