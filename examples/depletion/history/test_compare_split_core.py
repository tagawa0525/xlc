"""compare_split_core.py（途中で条件を変えた燃焼の定数の側を、今の core の組成で U の外れと L の外れに分ける）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_split_core.py'
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_residual as r
import compare_split_core as c
import compare_switch as s
import compare_weight as w
from test_compare_switch import (
    N20,
    REFERENCE20,
    SUMMARY,
    cores,
    history,
    hybrid,
    reference,
    table,
)

# 履歴の補正の core の 20 MWd/kg の点（w = 1/2）
WEIGHT = 0.5


def rho(state):
    return r.rho(r.k_inf(state))


def tables(core=None):
    return {
        "nominal": table(),
        "switch": reference(),
        "switch_hybrid": hybrid(),
        "core": cores()[0] if core is None else core,
    }


def switch_summary(core=None):
    """compare_switch.py の要約（switch.json）。fixture の追う核種は U-235 と Pu-240 だけ"""
    corrected, plain = cores()
    with mock.patch.object(s, "GROUPS", ["U-235", "Pu-240"]):
        return s.compare(
            table(),
            reference(),
            history(),
            corrected if core is None else core,
            plain,
            SUMMARY,
            hybrid(),
        )


def j_state(t, b):
    """公称のテーブルの燃焼度の点 `b` の条件 j の状態の番号"""
    return t["points"].index([b, 0, 1, 0])


def parts_at_20():
    """20 MWd/kg の、参照の組成の U_switch と L_switch、core が足す w U(1) と (P / P_j) L(1)、参照の定数、core の数密度の線形の補正の差
    （モジュールの説明の式を、ツールを通さずに組み立てる）"""
    t, ref, hyb = table(), reference(), hybrid()
    p = j_state(t, 2)
    nominal, coeff = t["states"][p], r.coefficients(t["micro"]["states"][p])
    n0 = t["micro"]["densities"][2]
    linear = r.corrected(nominal, coeff, [a - b for a, b in zip(REFERENCE20, n0)])
    untracked = r.difference(ref["states"][1], hyb["states"][1])
    linearization = r.difference(hyb["states"][1], linear)
    rod = t["histories"][0]["rods"][0]
    production = lambda d: sum(x["nu_fission"][-1] * y for x, y in zip(coeff, d))
    index = production([a - b for a, b in zip(N20, n0)]) / production(
        [a - b for a, b in zip(rod["densities"][1], n0)]
    )
    return {
        "reference": ref["states"][1],
        "untracked": untracked,
        "linearization": linearization,
        "core_untracked": w.scaled(rod["untracked"][1], WEIGHT),
        "core_linearization": w.scaled(rod["linearization"][1], index),
        "index": index,
        "composition": r.difference(
            r.corrected(
                ref["states"][1], coeff, [a - b for a, b in zip(N20, REFERENCE20)]
            ),
            ref["states"][1],
        ),
    }


def with_core_k(k):
    """20 MWd/kg の core の k を `k` にした、履歴の補正の core"""
    core = cores()[0]
    core["points"][2]["k"] = k
    return core


class Gaps(unittest.TestCase):
    def test_gaps_are_what_the_core_adds_less_the_parts_of_the_switch(self):
        x = parts_at_20()
        row = c.compare(tables(), switch_summary())["rows"][1]
        self.assertAlmostEqual(row["weight"], WEIGHT, places=12)
        self.assertAlmostEqual(row["index_ratio"], x["index"], places=12)
        ref = x["reference"]
        untracked = r.difference(x["core_untracked"], x["untracked"])
        linearization = r.difference(x["core_linearization"], x["linearization"])
        self.assertAlmostEqual(
            row["untracked_gap_pcm"],
            rho(r.added(ref, untracked)) - rho(ref),
            places=9,
        )
        self.assertAlmostEqual(
            row["linearization_gap_pcm"],
            rho(r.added(ref, linearization)) - rho(ref),
            places=9,
        )
        self.assertAlmostEqual(
            row["gap_pcm"],
            rho(r.added(r.added(ref, untracked), linearization)) - rho(ref),
            places=9,
        )
        # core のノードの定数から外れを引いた効き
        node = r.added(
            r.added(r.added(ref, x["composition"]), x["core_untracked"]),
            x["core_linearization"],
        )
        node = r.added(
            r.added(node, w.scaled(x["untracked"], -1.0)),
            w.scaled(x["linearization"], -1.0),
        )
        self.assertAlmostEqual(
            row["node_untracked_gap_pcm"],
            rho(node) - rho(r.added(node, w.scaled(untracked, -1.0))),
            places=9,
        )
        self.assertAlmostEqual(
            row["node_linearization_gap_pcm"],
            rho(node) - rho(r.added(node, w.scaled(linearization, -1.0))),
            places=9,
        )

    def test_the_continued_point_adds_nothing_and_has_no_index(self):
        # w = 0 の点は core が U も L も足さない。続ける点の参照の組成は公称と同じでも、ハイブリッドの定数は別の解
        row = c.compare(tables(), switch_summary())["rows"][0]
        self.assertEqual(row["weight"], 0.0)
        self.assertIsNone(row["index_ratio"])
        self.assertIsNone(row["untracked_weight"])
        self.assertIsNone(row["linearization_weight"])
        ref, hyb = reference(), hybrid()
        untracked = r.difference(ref["states"][0], hyb["states"][0])
        self.assertAlmostEqual(
            row["untracked_gap_pcm"],
            rho(r.added(ref["states"][0], w.scaled(untracked, -1.0)))
            - rho(ref["states"][0]),
            places=9,
        )

    def test_the_parts_are_set_on_the_nominal_state_with_the_switch_densities(self):
        # 成分の大きさは compare_split.py の比べ方（公称の状態の定数に参照の数密度の線形の補正を足した定数）
        x = parts_at_20()
        t = table()
        p = t["points"].index([2, 0, 0, 0])
        base = r.corrected(
            t["states"][p],
            r.coefficients(t["micro"]["states"][p]),
            [a - b for a, b in zip(REFERENCE20, t["micro"]["densities"][2])],
        )
        row = c.compare(tables(), switch_summary())["rows"][1]
        for key, value in (
            ("untracked_pcm", x["untracked"]),
            ("scaled_untracked_pcm", x["core_untracked"]),
            ("linearization_pcm", x["linearization"]),
            ("scaled_linearization_pcm", x["core_linearization"]),
        ):
            self.assertAlmostEqual(row[key], w.effect(base, value), places=9)
        self.assertAlmostEqual(
            row["untracked_weight"],
            row["untracked_pcm"] / row["one_untracked_pcm"],
            places=12,
        )
        self.assertAlmostEqual(
            row["linearization_weight"],
            row["linearization_pcm"] / row["one_linearization_pcm"],
            places=12,
        )

    def test_branch_densities_at_the_nominal_ones_add_no_linearization(self):
        # core の `production_ratio` は、履歴の分岐の数密度が公称と同じ点で L を足さない（U は w で足す）
        nominal = table()
        rod = nominal["histories"][0]["rods"][0]
        rod["densities"][1] = list(nominal["micro"]["densities"][2])
        x = parts_at_20()
        row = c.compare(tables() | {"nominal": nominal}, switch_summary())["rows"][1]
        self.assertIsNone(row["index_ratio"])
        self.assertIsNone(row["linearization_weight"])
        self.assertEqual(row["scaled_linearization_pcm"], 0.0)
        ref = x["reference"]
        self.assertAlmostEqual(
            row["linearization_gap_pcm"],
            rho(r.added(ref, w.scaled(x["linearization"], -1.0))) - rho(ref),
            places=9,
        )
        self.assertAlmostEqual(
            row["untracked_gap_pcm"],
            rho(r.added(ref, r.difference(x["core_untracked"], x["untracked"])))
            - rho(ref),
            places=9,
        )


class Closure(unittest.TestCase):
    def test_the_remainder_splits_into_reconstruction_identity_additivity_and_basis(
        self,
    ):
        for row in c.compare(tables(), switch_summary())["rows"]:
            self.assertAlmostEqual(
                row["constant_pcm"],
                row["difference_pcm"] - row["composition_pcm"],
                places=12,
            )
            self.assertAlmostEqual(
                row["remainder_pcm"], row["constant_pcm"] - row["gap_pcm"], places=12
            )
            self.assertAlmostEqual(
                row["remainder_pcm"],
                row["reconstruct_pcm"]
                + row["identity_pcm"]
                + row["additivity_pcm"]
                + row["basis_pcm"],
                places=9,
            )
            # core の定数は参照の定数 + 数密度の線形の補正 + 両方の外れ（線形の補正は数密度の差に 1 次）
            self.assertLess(abs(row["identity_pcm"]), 1e-9)

    def test_a_core_at_the_k_of_its_constants_reconstructs_to_zero(self):
        x = parts_at_20()
        node = r.added(
            r.added(r.added(x["reference"], x["composition"]), x["core_untracked"]),
            x["core_linearization"],
        )
        node = r.added(
            r.added(node, w.scaled(x["untracked"], -1.0)),
            w.scaled(x["linearization"], -1.0),
        )
        core = with_core_k(r.k_inf(node))
        row = c.compare(tables(core), switch_summary(core))["rows"][1]
        self.assertLess(abs(row["reconstruct_pcm"]), 1e-9)
        self.assertAlmostEqual(
            row["basis_pcm"],
            rho(x["reference"]) - r.rho(x["reference"]["k_inf"]),
            places=9,
        )

    def test_estimates_drop_one_gap_and_keep_the_densities(self):
        x = parts_at_20()
        ref = x["reference"]
        untracked = r.difference(x["core_untracked"], x["untracked"])
        linearization = r.difference(x["core_linearization"], x["linearization"])
        row = c.compare(tables(), switch_summary())["rows"][1]
        dense = r.added(ref, x["composition"])
        self.assertAlmostEqual(
            row["without_untracked_pcm"],
            rho(r.added(dense, linearization)) - r.rho(ref["k_inf"]),
            places=9,
        )
        self.assertAlmostEqual(
            row["without_linearization_pcm"],
            rho(r.added(dense, untracked)) - r.rho(ref["k_inf"]),
            places=9,
        )


class Inputs(unittest.TestCase):
    def test_a_switch_summary_of_other_values_is_an_error(self):
        for key in ("weight", "difference_pcm", "composition_pcm", "index_ratio"):
            summary = switch_summary()
            summary["rows"][1][key] += 1e-6
            with self.assertRaises(SystemExit, msg=key):
                c.compare(tables(), summary)

    def test_a_switch_summary_of_other_inputs_is_an_error(self):
        md5s = {f: f"md5-{f}" for f in c.FILES.values()}
        summary = {"inputs": {f: md5s[f] for f in s.FILES.values() if f in md5s}}
        c.check_summary(summary, md5s)
        summary["inputs"][c.FILES["core"]] = "other"
        with self.assertRaises(SystemExit):
            c.check_summary(summary, md5s)


class Report(unittest.TestCase):
    def test_the_report_names_both_gaps_and_the_check(self):
        text = c.report(c.compare(tables(), switch_summary()))
        self.assertIn("U の外れ", text)
        self.assertIn("L の外れ", text)
        self.assertIn("定数の側", text)
        self.assertNotIn("-", text)

    def test_later_points_without_a_correction_are_reported_with_dashes(self):
        summary = c.compare(tables(), switch_summary())
        summary["continued_burnup"] = 0.0
        text = c.report(summary)
        self.assertIn("P / P_j –、", text)


if __name__ == "__main__":
    unittest.main()
