"""compare_switch.py（途中で条件を変えた燃焼: 格子計算の参照の k∞ と、履歴の補正の core と補正を切った core の k の差）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_switch.py'
"""

import sys
import unittest
from pathlib import Path
from unittest import mock

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
REFERENCE20 = [5.2e-4, 2.6e-5]


def rate(a, b):
    """反応率の微視的断面積の差 a − b（lattice histories の形。収率の境界より下は割合 × 群の和の σf の差）"""
    below = lambda x: 0.9 * sum(x["fission"])
    return {
        "fission": [p - q for p, q in zip(a["fission"], b["fission"])],
        "reactions": [
            {
                "mt": 102,
                "channels": [
                    [
                        p - q
                        for p, q in zip(
                            a["reactions"][0]["channels"][0],
                            b["reactions"][0]["channels"][0],
                        )
                    ]
                ],
            }
        ],
        "thermal_fission": below(a) - below(b),
    }


# 履歴の分岐（一生を条件 j）の R の点 0、20 の、追う核種ごとの反応率の微視的断面積の差 Δσ_U（分岐 − ハイブリッド）と Δσ_L（ハイブリッド −
# 公称の条件 j）
RATE_UNTRACKED = [
    [
        rate(sigma(1.5, False), sigma(1.0, False)),
        rate(sigma(2.5, True), sigma(2.0, True)),
    ],
    [
        rate(sigma(4.5, False), sigma(4.25, False)),
        rate(sigma(8.5, True), sigma(8.25, True)),
    ],
]
RATE_LINEARIZATION = [
    [
        rate(sigma(1.0, False), sigma(1.0, False)),
        rate(sigma(2.0, True), sigma(2.0, True)),
    ],
    [
        rate(sigma(4.25, False), sigma(4.0, False)),
        rate(sigma(8.25, True), sigma(8.0, True)),
    ],
]


def table():
    """R を持つ公称のテーブルに、10 MWd/kg の条件 j の状態と、反応率の微視的断面積の差を足したもの"""
    t = residual_table()
    t["points"].append([1, 0, 1, 0])
    t["states"].append(few(2.6))
    t["micro"]["states"].append([sigma(7, False), sigma(8, True)])
    t["histories"][0]["rods"][0]["micro"] = {
        "untracked": RATE_UNTRACKED,
        "linearization": RATE_LINEARIZATION,
    }
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


def hybrid():
    """参照の追う核種だけを公称の履歴に入れた組成（ハイブリッド）を、参照と同じ点と条件 j で解いた分岐"""
    out = reference()
    out["states"] = [few(2.45), few(4.3)]
    out["micro"]["states"] = [
        [sigma(2.1, False), sigma(3.2, True)],
        [sigma(4.1, False), sigma(4.6, True)],
    ]
    return out


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
    """履歴の補正の core（k は 20 MWd/kg で 0.95）と、補正を切った core（0.96）"""
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
    return corrected, plain


SUMMARY = {
    "format": "xlc-compare-core-history",
    "branch_rows": [
        {"burnup": 0.0, "difference_pcm": 0.0, "nominal_difference_pcm": 0.0},
        {"burnup": 10.0, "difference_pcm": -1.0, "nominal_difference_pcm": -0.5},
        {"burnup": 20.0, "difference_pcm": -2.0, "nominal_difference_pcm": -1.5},
    ],
}


def compared():
    corrected, plain = cores()
    return c.compare(
        table(), reference(), history(), corrected, plain, SUMMARY, hybrid()
    )


def absorption(x):
    """核分裂と捕獲（MT 102）の和の、群ごとの微視的断面積"""
    return [f + ch for f, ch in zip(x["fission"], x["reactions"][0]["channels"][0])]


def rate_absorption(d):
    """反応率の微視的断面積の差の、核分裂と捕獲の和"""
    return [f + ch for f, ch in zip(d["fission"], d["reactions"][0]["channels"][0])]


class Compare(unittest.TestCase):
    def setUp(self):
        # fixture の追う核種は U-235 と Pu-240 だけ
        patcher = mock.patch.object(c, "GROUPS", ["U-235", "Pu-240"])
        patcher.start()
        self.addCleanup(patcher.stop)

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
        corrected, plain = cores()
        swapped = reference()
        swapped["points"] = [[1, 0, 0, 0], [0, 0, 0, 0]]
        swapped["states"] = swapped["states"][::-1]
        swapped["micro"]["states"] = swapped["micro"]["states"][::-1]
        swapped_hybrid = hybrid()
        swapped_hybrid["points"] = swapped["points"]
        swapped_hybrid["states"] = swapped_hybrid["states"][::-1]
        swapped_hybrid["micro"]["states"] = swapped_hybrid["micro"]["states"][::-1]
        s = c.compare(
            table(), swapped, history(), corrected, plain, SUMMARY, swapped_hybrid
        )
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
        corrected, plain = cores()
        same = dict(
            corrected,
            points=[
                dict(p, k=q["k"]) for p, q in zip(corrected["points"], plain["points"])
            ],
        )
        s = c.compare(table(), reference(), history(), same, plain, SUMMARY, hybrid())
        self.assertIsNone(s["rows"][1]["needed_weight"])
        # 補正の量が 0 なら、見込みは補正を切った差のまま
        self.assertEqual(
            s["omegas"][0]["estimate_pcm"][1], s["rows"][1]["plain_difference_pcm"]
        )
        c.report(s)

    def test_runs_that_do_not_match_the_reference_are_errors(self):
        corrected, plain = cores()
        h = hybrid()
        # 続けた点が参照の最初の点でない
        other = dict(history(), continued={"file": "x", "burnup": 5.0})
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), other, corrected, plain, SUMMARY, h)
        # 続けない履歴
        plain_history = {k: v for k, v in history().items() if k != "continued"}
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), plain_history, corrected, plain, SUMMARY, h)
        # 参照の燃焼度の点を core が通らない
        short = dict(plain, points=plain["points"][:2])
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, short, SUMMARY, h)
        # 参照は燃焼度ごとに 1 状態
        two = reference()
        two["axes"][2]["points"] = [565.0, 600.0]
        with self.assertRaises(SystemExit):
            c.compare(table(), two, history(), corrected, plain, SUMMARY, h)
        # 一生を条件 j で燃焼させた比較の要約が参照の点を持たない
        missing = dict(SUMMARY, branch_rows=SUMMARY["branch_rows"][:2])
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, plain, missing, h)

    def test_the_composition_is_split_by_the_named_nuclides_and_the_rest(self):
        # 数密度の効きを、名前を挙げた核種ごと（その核種だけの差）と、残りの核種をまとめた差に分ける
        with mock.patch.object(c, "GROUPS", ["Pu-240"]):
            s = compared()
        last = s["rows"][1]
        t = table()
        coeff = r.coefficients(t["micro"]["states"][t["points"].index([2, 0, 1, 0])])
        base = reference()["states"][1]
        deltas = [a - b for a, b in zip(N20, REFERENCE20)]

        def effect(mask):
            only = [d if m else 0.0 for d, m in zip(deltas, mask)]
            return r.rho(r.k_inf(r.corrected(base, coeff, only))) - r.rho(r.k_inf(base))

        groups = last["composition_groups"]
        self.assertEqual([x["group"] for x in groups], ["Pu-240", "rest"])
        self.assertAlmostEqual(groups[0]["pcm"], effect([False, True]), places=9)
        self.assertAlmostEqual(groups[1]["pcm"], effect([True, False]), places=9)
        # 追わない核種を挙げるのはエラー（黙って 0 にしない）
        with mock.patch.object(c, "GROUPS", ["Pu-239"]), self.assertRaises(SystemExit):
            compared()

    def test_the_rate_residual_outliers_are_split_into_the_untracked_part_and_the_linearization(
        self,
    ):
        # core の反応率は σ(公称、条件 j) + w Δσ_U + (P / P_j) Δσ_L、参照は σ(参照)。外れを、ハイブリッドで分けた参照の
        # Δσ_U,switch = σ(参照) − σ(ハイブリッド) と Δσ_L,switch = σ(ハイブリッド) − σ(公称、条件 j) からの差に分ける（核分裂と
        # 捕獲の和、群ごと、参照の σ に対する %）
        with mock.patch.object(c, "GROUPS", ["U-235", "Pu-240"]):
            s = compared()
        last = s["rows"][1]
        t = table()
        point = t["points"].index([2, 0, 1, 0])
        coeff = r.coefficients(t["micro"]["states"][point])
        nominal = t["micro"]["densities"][2]
        branch = t["histories"][0]["rods"][0]["densities"][1]
        production = lambda d: sum(x["nu_fission"][-1] * y for x, y in zip(coeff, d))
        index = production([a - b for a, b in zip(N20, nominal)]) / production(
            [a - b for a, b in zip(branch, nominal)]
        )
        self.assertAlmostEqual(last["index_ratio"], index, places=12)
        w = last["weight"]
        self.assertGreater(abs(index - w), 0.05)
        outliers = last["rate_outliers"]
        self.assertEqual([x["nuclide"] for x in outliers], ["U-235", "Pu-240"])
        for i, x in enumerate(outliers):
            ref = absorption(reference()["micro"]["states"][1][i])
            hyb = absorption(hybrid()["micro"]["states"][1][i])
            nom = absorption(t["micro"]["states"][point][i])
            du = rate_absorption(RATE_UNTRACKED[1][i])
            dl = rate_absorption(RATE_LINEARIZATION[1][i])
            for g in range(2):
                self.assertAlmostEqual(
                    x["untracked_percent"][g],
                    (w * du[g] - (ref[g] - hyb[g])) / ref[g] * 100.0,
                    places=9,
                )
                self.assertAlmostEqual(
                    x["linearization_percent"][g],
                    (index * dl[g] - (hyb[g] - nom[g])) / ref[g] * 100.0,
                    places=9,
                )
            # 参照の Δσ_L,switch の Δσ_L に対する実効の比（群ごと）と、参照の数密度の差の履歴の分岐の差に対する比（この核種だけの指標）
            for g in range(2):
                self.assertAlmostEqual(
                    x["linearization_ratio"][g], (hyb[g] - nom[g]) / dl[g], places=12
                )
            self.assertAlmostEqual(
                x["density_ratio"],
                (REFERENCE20[i] - nominal[i]) / (branch[i] - nominal[i]),
                places=12,
            )
        # 続ける点（w = 0）は外れを持たない（core は補正を足さない）
        self.assertIsNone(s["rows"][0]["rate_outliers"])
        text = c.report(s)
        self.assertIn("Δσ_U", text)
        self.assertIn("Δσ_L", text)

    def test_a_hybrid_of_other_points_or_densities_is_an_error(self):
        corrected, plain = cores()
        other = hybrid()
        other["micro"]["densities"] = [N10, N20]
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, plain, SUMMARY, other)
        other = hybrid()
        other["axes"] = axes([10.0, 20.0], [565.0], [0.743], 565.0)
        with self.assertRaises(SystemExit):
            c.compare(table(), reference(), history(), corrected, plain, SUMMARY, other)


if __name__ == "__main__":
    unittest.main()
