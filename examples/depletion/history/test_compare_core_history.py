"""compare_core_history.py（履歴の補正の core の確認: 公称の出力の一致、条件 j の k、途中で条件を変えた燃焼の w と R の効き）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_core_history.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_history as c
import compare_residual as r
from test_compare_residual import NUCLIDES, axes, few, history_rod, residual
from test_compare_residual import sigma as plain_sigma


def sigma(s: float, derivative: bool) -> dict:
    """2 群の微視的断面積に、捕獲（MT 102）の経路 1 つを足したもの"""
    out = plain_sigma(s, derivative)
    out["reactions"] = [
        {"mt": 102, "channels": [[0.3 * (1.0 + s / 5.0), 20.0 * (1.0 + s / 5.0)]]}
    ]
    return out


STATES = {
    "nominal": {
        "boron": 1300.0,
        "b10": 0.199,
        "moderator": 565.0,
        "fuel": 565.0,
        "banks": {},
    },
    "branch": {
        "boron": 1300.0,
        "b10": 0.199,
        "moderator": 600.0,
        "fuel": 565.0,
        "banks": {},
    },
}


def table():
    """公称のテーブル（燃焼度 0、10、20。減速材の公称 565 K と条件 600 K。条件は 0 と 20 で解く）と、R（燃焼度 0 と 20）"""
    points = [[0, 0, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [2, 0, 0, 0], [2, 0, 1, 0]]
    states = [few(p) for p in range(5)]
    micro = {
        "nuclides": NUCLIDES,
        "densities": [[7.0e-4, 0.0], [6.0e-4, 1.0e-5], [5.0e-4, 2.0e-5]],
        "states": [[sigma(p, False), sigma(2 * p, True)] for p in range(5)],
    }
    rs = []
    for p, u, nj, fb in (
        (1, 0, [7.0e-4, 0.0], few(1.0)),
        (4, 2, [5.2e-4, 2.6e-5], few(4.5)),
    ):
        deltas = [a - b for a, b in zip(nj, micro["densities"][u])]
        rs.append(
            residual(
                fb, r.corrected(states[p], r.coefficients(micro["states"][p]), deltas)
            )
        )
    return {
        "format": "xlc-table",
        "axes": axes([0.0, 10.0, 20.0], [565.0, 600.0], [0.743, 0.661], 565.0),
        "points": points,
        "states": states,
        "micro": micro,
        "histories": [
            {
                "variable": "moderator-density",
                "nominal": 0.743,
                "branch": 0.661,
                "rods": [
                    history_rod([0.0, 20.0], rs, [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]])
                ],
            }
        ],
    }


def branch_table():
    """第 1 項のテーブル（600 K、0.661 g/cc。燃焼度 0 と 20）"""
    return {
        "format": "xlc-table",
        "axes": axes([0.0, 20.0], [600.0], [0.661], 600.0),
        "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
        "states": [few(1.0), few(4.5)],
        "micro": {
            "nuclides": NUCLIDES,
            "densities": [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]],
            "states": [
                [sigma(9, False), sigma(9, False)],
                [sigma(10, False), sigma(11, False)],
            ],
        },
    }


def core(points):
    """集合体 1 体・ノード 1 つの core の履歴。`points` は (EFPD, 燃焼度, 状態, k, 数密度, H)"""
    return {
        "format": "xlc-core-history",
        "assemblies": ["A01"],
        "axial_mesh": [0.0, 21.5],
        "states": STATES,
        "nuclides": NUCLIDES,
        "history_variables": ["moderator-density"],
        "points": [
            {
                "efpd": e,
                "burnup": b,
                "state": s,
                "k": k,
                "burnups": [[b]],
                "densities": [[n]],
                "histories": [[[h]]],
            }
            for e, b, s, k, n, h in points
        ],
    }


def scaled(x, w):
    return {
        k: [w * a for a in v]
        if not isinstance(v[0], list)
        else [[w * a for a in row] for row in v]
        for k, v in x.items()
    }


def node_k(t, b, m, n, w):
    """テーブル `t` の燃焼度の点 `b`、減速材の点 `m` の定数 + 線形の補正（数密度 `n`）+ w U + (P / P_j) L の 2 群の k。P は熱群の νΣf の
    線形の補正（ノードの数密度の差）、P_j は同じ係数の履歴の分岐の数密度の差のもの。履歴の分岐の数密度が公称と同じ点は L を足さない"""
    p = t["points"].index([b, 0, m, 0])
    densities = t["micro"]["densities"][b]
    coefficients = r.coefficients(t["micro"]["states"][p])
    deltas = [x - y for x, y in zip(n, densities)]
    f = r.corrected(t["states"][p], coefficients, deltas)
    rod = t["histories"][0]["rods"][0]
    u = t["axes"][0]["points"][b]
    if w != 0.0:
        k = rod["burnup"].index(u)
        f = r.added(f, scaled(rod["untracked"][k], w))
        branch = [x - y for x, y in zip(rod["densities"][k], densities)]
        if any(d != 0.0 for d in branch):
            own = sum(c["nu_fission"][-1] * d for c, d in zip(coefficients, deltas))
            of_branch = sum(
                c["nu_fission"][-1] * d for c, d in zip(coefficients, branch)
            )
            f = r.added(f, scaled(rod["linearization"][k], own / of_branch))
    return r.k_inf(f)


# 第 1 項のテーブルで条件 j の燃焼をした core の 20 MWd/kg の数密度（第 1 項の N̄ᵢ^j = [5.2e-4, 2.6e-5] から +0.5% と −2%）
SIGMA20 = [5.2e-4 * 1.005, 2.6e-5 * 0.98]


def sigma_core():
    """第 1 項のテーブルで、条件 j で一生を燃焼させた core（R の点 0 と 20 の k は 1.2 と 0.95）"""
    n0 = [7.0e-4, 0.0]
    return core(
        [
            (0.0, 0.0, "branch", 1.2, n0, 0.661),
            (5.0, 10.0, "branch", 1.0, [6.0e-4, 1.0e-5], 0.661),
            (10.0, 20.0, "branch", 0.95, SIGMA20, 0.661),
        ]
    )


SUMMARY = {
    "rows": [{"burnup": 0.0, "applied_pcm": 0.0}, {"burnup": 20.0, "applied_pcm": -6.0}]
}
N20 = [5.25e-4, 2.5e-5]


class Compare(unittest.TestCase):
    def cores(self):
        t = table()
        n0 = [7.0e-4, 0.0]
        nominal = core(
            [
                (0.0, 0.0, "nominal", node_k(t, 0, 0, n0, 0.0), n0, 0.743),
                (5.0, 10.0, "nominal", 1.0, [6.0e-4, 1.0e-5], 0.743),
                (10.0, 20.0, "nominal", 0.9, N20, 0.743),
            ]
        )
        branch = core(
            [
                (0.0, 0.0, "branch", node_k(t, 0, 1, n0, 1.0), n0, 0.661),
                (5.0, 10.0, "branch", 1.0, [6.0e-4, 1.0e-5], 0.661),
                (10.0, 20.0, "branch", node_k(t, 2, 1, N20, 1.0), N20, 0.661),
            ]
        )
        # 10 MWd/kg まで公称、その後を条件 j: 20 MWd/kg の H = (0.743 × 10 + 0.661 × 10) / 20
        h = (0.743 * 10.0 + 0.661 * 10.0) / 20.0
        w = (h - 0.743) / (0.661 - 0.743)
        switch = core(
            [
                (0.0, 0.0, "nominal", node_k(t, 0, 0, n0, 0.0), n0, 0.743),
                (5.0, 10.0, "branch", 1.0, [6.0e-4, 1.0e-5], 0.743),
                (10.0, 20.0, "branch", node_k(t, 2, 1, N20, w) * (1.0 + 1e-7), N20, h),
            ]
        )
        return t, nominal, branch, switch, w

    def test_the_branch_core_is_compared_with_the_first_term_and_the_nominal_with_the_table(
        self,
    ):
        t, nominal, branch, switch, _ = self.cores()
        s = c.compare(
            t, branch_table(), nominal, branch, switch, sigma_core(), SUMMARY, True
        )
        self.assertEqual(s["format"], "xlc-compare-core-history")
        self.assertTrue(s["identical"])
        rows = s["branch_rows"]
        # R の燃焼度の点（0、20）だけ
        self.assertEqual([x["burnup"] for x in rows], [0.0, 20.0])
        last = rows[1]
        k_table = branch_table()["states"][1]["k_inf"]
        self.assertAlmostEqual(
            last["difference_pcm"],
            r.rho(branch["points"][2]["k"]) - r.rho(k_table),
            places=9,
        )
        # 公称の燃焼の core と、公称のテーブルの公称の状態の k∞
        self.assertAlmostEqual(
            last["nominal_difference_pcm"],
            r.rho(0.9) - r.rho(t["states"][3]["k_inf"]),
            places=9,
        )
        # core の k は、格子点の定数 + 線形の補正 + R の 2 群の k
        for x in rows:
            self.assertAlmostEqual(x["reconstruct_pcm"], 0.0, places=9)

    def test_the_branch_difference_is_split_by_the_tracked_densities(self):
        # core の数密度と第 1 項のテーブルの数密度の違いの効き（核種ごと）と、第 1 項のテーブルの微視的断面積の公称のテーブルの
        # 条件 j の値との相対差
        t, nominal, branch, switch, _ = self.cores()
        nj = branch_table()["micro"]["densities"][1]
        s = c.compare(
            t, branch_table(), nominal, branch, switch, sigma_core(), SUMMARY, True
        )
        last = s["branch_rows"][1]
        p = t["points"].index([2, 0, 1, 0])
        micro = t["micro"]
        coeff = r.coefficients(micro["states"][p])
        base = r.added(
            r.corrected(
                t["states"][p],
                coeff,
                [a - b for a, b in zip(nj, micro["densities"][2])],
            ),
            r.residual(t["histories"][0]["rods"][0], 1),
        )
        whole = r.corrected(base, coeff, [a - b for a, b in zip(N20, nj)])
        self.assertAlmostEqual(
            last["composition_pcm"],
            r.rho(r.k_inf(whole)) - r.rho(r.k_inf(base)),
            places=9,
        )
        names = [x["nuclide"] for x in last["nuclides"]]
        self.assertEqual(sorted(names), sorted(NUCLIDES))
        u235 = next(x for x in last["nuclides"] if x["nuclide"] == "U-235")
        self.assertAlmostEqual(u235["density_percent"], (N20[0] / nj[0] - 1.0) * 100.0)
        only = r.corrected(base, coeff, [N20[0] - nj[0], 0.0])
        self.assertAlmostEqual(
            u235["pcm"], r.rho(r.k_inf(only)) - r.rho(r.k_inf(base)), places=9
        )
        # 微視的断面積の相対差（第 1 項 / 公称の条件 j − 1、群ごと。捕獲は MT 102 の経路の和）
        first = branch_table()["micro"]["states"][1][0]
        own = micro["states"][p][0]
        self.assertAlmostEqual(
            u235["fission_percent"][1],
            (first["fission"][1] / own["fission"][1] - 1.0) * 100.0,
        )
        capture = lambda x: sum(
            ch[0] for m in x["reactions"] if m["mt"] == 102 for ch in m["channels"]
        )
        self.assertAlmostEqual(
            u235["capture_percent"][0], (capture(first) / capture(own) - 1.0) * 100.0
        )
        # 燃焼度 0 の数密度は同じ
        self.assertAlmostEqual(s["branch_rows"][0]["composition_pcm"], 0.0, places=12)

    def test_the_branch_cores_are_compared_with_the_first_term_densities(self):
        # 履歴の補正の core と、第 1 項のテーブルで燃焼させた core の、k の第 1 項の k∞ からの差と、数密度の第 1 項の N̄ᵢ^j からの
        # 相対差の絶対値の最大（N̄ᵢ^j が 0 の核種は除く）
        t, nominal, branch, switch, _ = self.cores()
        s = c.compare(
            t, branch_table(), nominal, branch, switch, sigma_core(), SUMMARY, True
        )
        first, last = s["branch_rows"]
        k_table = branch_table()["states"][1]["k_inf"]
        self.assertAlmostEqual(
            last["sigma_difference_pcm"], r.rho(0.95) - r.rho(k_table), places=9
        )
        self.assertEqual(last["sigma_k"], 0.95)
        self.assertAlmostEqual(last["sigma_density_percent"], -2.0, places=9)
        self.assertEqual(last["sigma_density_nuclide"], "Pu-240")
        # 履歴の補正の core は N20 = [5.25e-4, 2.5e-5]（U-235 +0.96%、Pu-240 −3.8%）
        nj = branch_table()["micro"]["densities"][1]
        self.assertAlmostEqual(
            last["density_percent"], (N20[1] / nj[1] - 1.0) * 100.0, places=9
        )
        self.assertEqual(last["density_nuclide"], "Pu-240")
        # 燃焼度 0 は Pu-240 が 0 なので U-235 だけ
        self.assertEqual(first["density_nuclide"], "U-235")
        self.assertAlmostEqual(first["density_percent"], 0.0, places=12)
        text = c.report(s)
        self.assertIn("第 1 項のテーブルで燃焼させた core", text)

    def test_the_switched_core_reports_the_weight_and_the_residual_effect(self):
        t, nominal, branch, switch, w = self.cores()
        s = c.compare(
            t, branch_table(), nominal, branch, switch, sigma_core(), SUMMARY, False
        )
        self.assertFalse(s["identical"])
        self.assertEqual(s["applied_max_pcm"], 6.0)
        rows = s["switch_rows"]
        self.assertEqual([x["burnup"] for x in rows], [0.0, 10.0, 20.0])
        self.assertEqual([x["state"] for x in rows], ["nominal", "branch", "branch"])
        self.assertEqual(rows[1]["weight"], 0.0)
        self.assertAlmostEqual(rows[2]["weight"], w, places=15)
        self.assertAlmostEqual(rows[2]["bound_pcm"], w * 6.0, places=12)
        # R の点の効き ρ(f + w R) − ρ(f) と、その点の誤差の上限
        self.assertAlmostEqual(
            rows[2]["residual_pcm"],
            r.rho(node_k(t, 2, 1, N20, w)) - r.rho(node_k(t, 2, 1, N20, 0.0)),
            places=9,
        )
        self.assertAlmostEqual(rows[2]["applied_bound_pcm"], w * 6.0, places=12)
        # core の k は 1e-7 だけ大きい
        self.assertAlmostEqual(
            rows[2]["reconstruct_pcm"],
            r.rho(node_k(t, 2, 1, N20, w)) - r.rho(switch["points"][2]["k"]),
            places=9,
        )
        self.assertLess(rows[2]["reconstruct_pcm"], -1e-3)
        # R の点でない燃焼点は効きを組み立てない
        self.assertIsNone(rows[1]["residual_pcm"])
        self.assertEqual(rows[0]["residual_pcm"], 0.0)
        self.assertIn("w", c.report(s))

    def test_cores_that_do_not_match_the_table_are_errors(self):
        t, nominal, branch, switch, _ = self.cores()
        other = dict(branch, nuclides=["U-235", "Pu-239"])
        with self.assertRaises(SystemExit):
            c.compare(
                t, branch_table(), nominal, other, switch, sigma_core(), SUMMARY, True
            )
        # R の点で、状態の減速材の温度がテーブルの格子点でない
        off = dict(
            branch, states=dict(STATES, branch=dict(STATES["branch"], moderator=590.0))
        )
        with self.assertRaises(SystemExit):
            c.compare(
                t, branch_table(), nominal, off, switch, sigma_core(), SUMMARY, True
            )
        # 公称の燃焼と条件 j の燃焼の R の点が食い違う
        shifted = dict(nominal, points=nominal["points"][:2])
        with self.assertRaises(SystemExit):
            c.compare(
                t, branch_table(), shifted, branch, switch, sigma_core(), SUMMARY, True
            )
        # 第 1 項のテーブルが R の点を解いていない
        short = branch_table()
        short["axes"][0]["points"] = [0.0, 15.0]
        with self.assertRaises(SystemExit):
            c.compare(t, short, nominal, branch, switch, sigma_core(), SUMMARY, True)

    def test_the_residual_summary_must_cover_the_residual_points_of_the_table(self):
        t, nominal, branch, switch, _ = self.cores()
        missing = {"rows": [SUMMARY["rows"][0]]}
        with self.assertRaises(SystemExit):
            c.compare(
                t, branch_table(), nominal, branch, switch, sigma_core(), missing, True
            )
        # 燃焼度は相対 1e-9 で照らす
        near = {
            "rows": [
                dict(x, burnup=x["burnup"] * (1.0 + 1e-12)) for x in SUMMARY["rows"]
            ]
        }
        s = c.compare(
            t, branch_table(), nominal, branch, switch, sigma_core(), near, True
        )
        self.assertAlmostEqual(
            s["switch_rows"][2]["applied_bound_pcm"],
            s["switch_rows"][2]["weight"] * 6.0,
        )
        # 要約は同じ R を持つテーブルから作ったもの
        c.check_summary({"inputs": {c.FILES["residual_table"]: "abc"}}, "abc")
        with self.assertRaises(SystemExit):
            c.check_summary({"inputs": {c.FILES["residual_table"]: "abc"}}, "abd")

    def test_reports_of_runs_without_a_weight_or_late_points_are_written(self):
        t, nominal, branch, _, _ = self.cores()
        s = c.compare(
            t, branch_table(), nominal, branch, nominal, sigma_core(), SUMMARY, True
        )
        s["branch_rows"] = [x for x in s["branch_rows"] if x["burnup"] < 15.0]
        text = c.report(s)
        self.assertIn("w", text)


if __name__ == "__main__":
    unittest.main()
