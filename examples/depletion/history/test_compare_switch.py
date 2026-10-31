"""compare_switch.py（途中で条件を変えた燃焼: 格子計算の参照の k∞ と、履歴の補正の core と補正を切った core の k の差）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_switch.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_residual as r
import compare_switch as c
from test_compare_core_history import core, sigma
from test_compare_core_history import table as residual_table
from test_compare_residual import NUCLIDES, axes, few

# 10 MWd/kg まで公称、その後を条件 j。参照は 10 と 20 MWd/kg
N10 = [6.0e-4, 1.0e-5]
N20 = [5.25e-4, 2.5e-5]
PLAIN20 = [5.3e-4, 2.4e-5]
INDEX20 = [5.22e-4, 2.55e-5]
REFERENCE20 = [5.2e-4, 2.6e-5]


def table():
    """R を持つ公称のテーブルに、10 MWd/kg の条件 j の状態を足したもの"""
    t = residual_table()
    t["points"].append([1, 0, 1, 0])
    t["states"].append(few(2.6))
    t["micro"]["states"].append([sigma(7, False), sigma(8, True)])
    return t


def reference():
    """参照（途中で条件を変えた燃焼の履歴の 10、20 MWd/kg を条件 j で解いた分岐。600 K、0.661 g/cc）"""
    return {
        "format": "xlc-table",
        "axes": axes([10.0, 20.0], [600.0], [0.661], 600.0),
        "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
        "states": [few(2.5), few(4.4)],
        "micro": {
            "nuclides": NUCLIDES,
            "densities": [N10, REFERENCE20],
            "states": [
                [sigma(2, False), sigma(3, True)],
                [sigma(4, False), sigma(5, True)],
            ],
        },
    }


def history():
    """途中で条件を変えた格子計算の燃焼の履歴（10 MWd/kg の点から続けた）"""
    return {
        "format": "xlc-depletion-history",
        "points": [
            {"burnup": 0.0, "k_inf": 1.2},
            {"burnup": 10.0, "k_inf": few(2.5)["k_inf"] * (1.0 + 2e-7)},
            {"burnup": 20.0, "k_inf": few(4.4)["k_inf"] * (1.0 - 1e-7)},
        ],
        "continued": {"file": "../p2/lattice-2A.msgpack", "burnup": 10.0},
    }


H20 = (0.743 * 10.0 + 0.661 * 10.0) / 20.0


def cores():
    """履歴の補正の core（k は 20 MWd/kg で 0.95）、補正を切った core（0.96）と、反応率の Δσ を P / P_j で足した core（0.952）"""
    n0 = [7.0e-4, 0.0]
    corrected = core(
        [
            (0.0, 0.0, "nominal", 1.2, n0, 0.743),
            (5.0, 10.0, "branch", 1.01, N10, 0.743),
            (10.0, 20.0, "branch", 0.95, N20, H20),
        ]
    )
    plain = core(
        [
            (0.0, 0.0, "nominal", 1.2, n0, 0.743),
            (5.0, 10.0, "branch", 1.01, N10, 0.743),
            (10.0, 20.0, "branch", 0.96, PLAIN20, H20),
        ]
    )
    index = core(
        [
            (0.0, 0.0, "nominal", 1.2, n0, 0.743),
            (5.0, 10.0, "branch", 1.01, N10, 0.743),
            (10.0, 20.0, "branch", 0.952, INDEX20, H20),
        ]
    )
    return corrected, plain, index


SUMMARY = {
    "format": "xlc-compare-core-history",
    "branch_rows": [
        {"burnup": 0.0, "difference_pcm": 0.0, "nominal_difference_pcm": 0.0},
        {"burnup": 10.0, "difference_pcm": -1.0, "nominal_difference_pcm": -0.5},
        {"burnup": 20.0, "difference_pcm": -2.0, "nominal_difference_pcm": -1.5},
    ],
}


def compared():
    corrected, plain, index = cores()
    return c.compare(table(), reference(), history(), corrected, plain, SUMMARY, index)


class Compare(unittest.TestCase):
    def test_the_cores_are_compared_with_the_reference_at_its_burnups(self):
        s = compared()
        self.assertEqual(s["format"], "xlc-compare-history-switch")
        self.assertEqual(s["continued_burnup"], 10.0)
        rows = s["rows"]
        self.assertEqual([x["burnup"] for x in rows], [10.0, 20.0])
        last = rows[1]
        k_reference = few(4.4)["k_inf"]
        self.assertEqual(last["reference_k_inf"], k_reference)
        self.assertAlmostEqual(
            last["difference_pcm"], r.rho(0.95) - r.rho(k_reference), places=9
        )
        self.assertAlmostEqual(
            last["plain_difference_pcm"], r.rho(0.96) - r.rho(k_reference), places=9
        )
        self.assertEqual((last["k"], last["plain_k"]), (0.95, 0.96))
        # w = (H − H⁰) / (Hᵇ − H⁰) は 20 MWd/kg で 1/2、続ける点で 0
        self.assertAlmostEqual(last["weight"], 0.5, places=12)
        self.assertEqual(rows[0]["weight"], 0.0)
        # 参照の分岐は格子計算の履歴の点を再現する
        self.assertAlmostEqual(
            last["reference_history_pcm"],
            r.rho(k_reference) - r.rho(history()["points"][2]["k_inf"]),
            places=9,
        )
        # 一生を条件 j で燃焼させた core と公称の燃焼の core の差を同じ点に並べる
        self.assertEqual(last["branch_difference_pcm"], -2.0)
        self.assertEqual(last["nominal_difference_pcm"], -1.5)

    def test_the_continued_point_matches_the_nominal_table_at_condition_j(self):
        # 続ける点の参照は、公称の組成を条件 j で解いたもので、公称のテーブルの条件 j の状態と同じ
        s = compared()
        t = table()
        k_table = t["states"][t["points"].index([1, 0, 1, 0])]["k_inf"]
        self.assertAlmostEqual(
            s["check_pcm"], r.rho(few(2.5)["k_inf"]) - r.rho(k_table), places=9
        )

    def test_the_difference_is_split_by_the_tracked_densities(self):
        # 参照の定数に、core の数密度と参照の数密度の差の線形の補正を足した効き（参照の係数）
        s = compared()
        last = s["rows"][1]
        ref = reference()
        # 係数は core と同じ、R を持つ公称のテーブルの同じ燃焼度の条件 j の状態のもの（微分を持つ核種は微分）
        t = table()
        coeff = r.coefficients(t["micro"]["states"][t["points"].index([2, 0, 1, 0])])
        base = ref["states"][1]
        for key, n in (("composition_pcm", N20), ("plain_composition_pcm", PLAIN20)):
            whole = r.corrected(base, coeff, [a - b for a, b in zip(n, REFERENCE20)])
            self.assertAlmostEqual(
                last[key], r.rho(r.k_inf(whole)) - r.rho(r.k_inf(base)), places=9
            )
        # 相対差（%）の絶対値が最大の核種と、主な核種の相対差
        self.assertEqual(last["density_nuclide"], "Pu-240")
        self.assertAlmostEqual(
            last["density_percent"], (N20[1] / REFERENCE20[1] - 1.0) * 100.0, places=9
        )
        self.assertAlmostEqual(
            last["densities_percent"]["U-235"],
            (N20[0] / REFERENCE20[0] - 1.0) * 100.0,
            places=9,
        )
        self.assertAlmostEqual(
            last["plain_densities_percent"]["Pu-240"],
            (PLAIN20[1] / REFERENCE20[1] - 1.0) * 100.0,
            places=9,
        )
        self.assertEqual(set(last["densities_percent"]), set(NUCLIDES))
        # 続ける点は同じ組成
        self.assertAlmostEqual(s["rows"][0]["composition_pcm"], 0.0, places=12)
        text = c.report(s)
        self.assertIn("参照", text)
        self.assertIn("補正を切った", text)

    def test_the_composition_is_split_by_nuclide(self):
        # 数密度の効きを核種ごと（その核種だけの差）に分け、効きの大きい順に並べる
        s = compared()
        last = s["rows"][1]
        ref = reference()
        t = table()
        coeff = r.coefficients(t["micro"]["states"][t["points"].index([2, 0, 1, 0])])
        base = ref["states"][1]
        expected = {}
        for i, name in enumerate(NUCLIDES):
            only = [
                a - b if j == i else 0.0
                for j, (a, b) in enumerate(zip(N20, REFERENCE20))
            ]
            expected[name] = r.rho(r.k_inf(r.corrected(base, coeff, only))) - r.rho(
                r.k_inf(base)
            )
        got = last["composition_nuclides"]
        self.assertEqual(
            [x["nuclide"] for x in got],
            sorted(NUCLIDES, key=lambda n: -abs(expected[n])),
        )
        for x in got:
            self.assertAlmostEqual(x["pcm"], expected[x["nuclide"]], places=9)
        text = c.report(s)
        # 数密度の効きを除いた残り（定数の側）と、最後の点の核種ごとの効き
        rest = last["difference_pcm"] - last["composition_pcm"]
        self.assertIn(f"{rest:+.3f}".replace("-", "−"), text.replace("-", "−"))
        self.assertIn(got[0]["nuclide"], text)

    def test_weights_of_other_omegas_estimate_the_difference(self):
        # core の H（ω = 1）から各ステップの x を戻し、H(B₀ + ΔB) = (H(B₀) B₀ + ω x ΔB) / (B₀ + ωΔB) で進め直した w と、補正の量が w に
        # 比例すると置いた差の見込み: 補正を切った差 − w(ω) / w(1) × (補正を切った差 − 補正した差)
        s = compared()
        # 0〜10 MWd/kg は公称（x = 0.743、w の x は 0）、10〜20 は条件 j（w の x は 1）
        omega = next(x for x in s["omegas"] if x["omega"] == 2.0)
        w2 = (0.0 * 10.0 + 2.0 * 1.0 * 10.0) / (10.0 + 2.0 * 10.0)
        self.assertEqual(omega["weights"][0], 0.0)
        self.assertAlmostEqual(omega["weights"][1], w2, places=12)
        last = s["rows"][1]
        fill = last["plain_difference_pcm"] - last["difference_pcm"]
        self.assertAlmostEqual(
            omega["estimate_pcm"][1],
            last["plain_difference_pcm"] - w2 / 0.5 * fill,
            places=9,
        )
        self.assertIsNone(omega["estimate_pcm"][0])
        # 参照に合う w（補正の量が w に比例すると置いて）
        self.assertAlmostEqual(
            last["needed_weight"], 0.5 * last["plain_difference_pcm"] / fill, places=12
        )
        self.assertIsNone(s["rows"][0]["needed_weight"])
        self.assertEqual([x["omega"] for x in s["omegas"]], [1.25, 1.5, 2.0, 2.5])
        self.assertIn("ω = 2", c.report(s))

    def test_reference_densities_follow_the_burnup_of_each_state(self):
        # N̄ᵢ は燃焼度の点の順、状態は points の順。状態を燃焼度の逆に並べても同じ要約になる
        corrected, plain, index = cores()
        swapped = reference()
        swapped["points"] = [[1, 0, 0, 0], [0, 0, 0, 0]]
        swapped["states"] = swapped["states"][::-1]
        swapped["micro"]["states"] = swapped["micro"]["states"][::-1]
        s = c.compare(table(), swapped, history(), corrected, plain, SUMMARY, index)
        self.assertEqual(s["rows"], compared()["rows"])

    def test_the_core_history_summary_is_made_from_the_same_table_and_core(self):
        md5s = {
            c.FILES["residual_table"]: "a",
            c.FILES["core"]: "b",
            c.FILES["plain_table"]: "c",
        }
        c.check_core_history({"inputs": dict(md5s)}, md5s)
        for key in md5s:
            with self.assertRaises(SystemExit):
                c.check_core_history({"inputs": dict(md5s, **{key: "x"})}, md5s)

    def test_points_where_the_correction_does_nothing_have_no_needed_weight(self):
        corrected, plain, index = cores()
        same = dict(
            corrected,
            points=[
                dict(p, k=q["k"]) for p, q in zip(corrected["points"], plain["points"])
            ],
        )
        s = c.compare(table(), reference(), history(), same, plain, SUMMARY, index)
        self.assertIsNone(s["rows"][1]["needed_weight"])
        # 補正の量が 0 なら、見込みは補正を切った差のまま
        self.assertEqual(
            s["omegas"][0]["estimate_pcm"][1], s["rows"][1]["plain_difference_pcm"]
        )
        c.report(s)

    def test_runs_that_do_not_match_the_reference_are_errors(self):
        corrected, plain, index = cores()
        # 続けた点が参照の最初の点でない
        other = dict(history(), continued={"file": "x", "burnup": 5.0})
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), other, corrected, plain, SUMMARY, index)
        # 続けない履歴
        plain_history = {k: v for k, v in history().items() if k != "continued"}
        with self.assertRaises(SystemExit):
            c.compare(
                table(), reference(), plain_history, corrected, plain, SUMMARY, index
            )
        # 参照の燃焼度の点を core が通らない
        short = dict(plain, points=plain["points"][:2])
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, short, SUMMARY, index)
        # 参照は燃焼度ごとに 1 状態
        two = reference()
        two["axes"][2]["points"] = [565.0, 600.0]
        with self.assertRaises(SystemExit):
            c.compare(table(), two, history(), corrected, plain, SUMMARY, index)
        # 一生を条件 j で燃焼させた比較の要約が参照の点を持たない
        missing = dict(SUMMARY, branch_rows=SUMMARY["branch_rows"][:2])
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, plain, missing, index)

    def test_the_core_with_the_rate_residuals_at_the_production_index_is_compared_too(
        self,
    ):
        # 反応率の Δσ を P / P_j で足した core（`rate-residual-weight = "production-index"`）も、参照からの差と数密度の効きに分ける
        s = compared()
        last = s["rows"][1]
        k_reference = few(4.4)["k_inf"]
        self.assertEqual(last["index_k"], 0.952)
        self.assertAlmostEqual(
            last["index_difference_pcm"], r.rho(0.952) - r.rho(k_reference), places=9
        )
        t = table()
        coeff = r.coefficients(t["micro"]["states"][t["points"].index([2, 0, 1, 0])])
        base = reference()["states"][1]
        whole = r.corrected(base, coeff, [a - b for a, b in zip(INDEX20, REFERENCE20)])
        self.assertAlmostEqual(
            last["index_composition_pcm"],
            r.rho(r.k_inf(whole)) - r.rho(r.k_inf(base)),
            places=9,
        )
        self.assertEqual(
            {x["nuclide"] for x in last["index_composition_nuclides"]}, set(NUCLIDES)
        )
        self.assertEqual(last["index_density_nuclide"], "Pu-240")
        self.assertAlmostEqual(
            last["index_densities_percent"]["U-235"],
            (INDEX20[0] / REFERENCE20[0] - 1.0) * 100.0,
            places=9,
        )
        self.assertIn("P / P_j", c.report(s))

    def test_a_production_index_core_of_another_history_is_an_error(self):
        # 2 つの core は同じ燃焼（同じ H）で、Δσ の重みだけが違う
        corrected, plain, index = cores()
        other = dict(
            index,
            points=[
                dict(p, histories=[[[0.7]]]) if i == 2 else p
                for i, p in enumerate(index["points"])
            ],
        )
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, plain, SUMMARY, other)


if __name__ == "__main__":
    unittest.main()
