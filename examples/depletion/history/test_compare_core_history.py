"""compare_core_history.py（履歴の補正の core の確認: 公称の出力の一致、条件 j の k、途中で条件を変えた燃焼の w と R の効き）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_core_history.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_history as c
import compare_residual as r
from test_compare_residual import NUCLIDES, axes, few, residual, sigma

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
                "rods": [{"rod": "out", "burnup": [0.0, 20.0], "residuals": rs}],
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
            "states": [[sigma(9, False), sigma(9, False)]] * 2,
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


def node_k(t, b, m, n, w):
    """テーブル `t` の燃焼度の点 `b`、減速材の点 `m` の定数 + 線形の補正（数密度 `n`）+ w R の 2 群の k"""
    p = t["points"].index([b, 0, m, 0])
    deltas = [x - y for x, y in zip(n, t["micro"]["densities"][b])]
    f = r.corrected(t["states"][p], r.coefficients(t["micro"]["states"][p]), deltas)
    rod = t["histories"][0]["rods"][0]
    u = t["axes"][0]["points"][b]
    if w != 0.0:
        scaled = {
            k: [w * x for x in v]
            if not isinstance(v[0], list)
            else [[w * x for x in row] for row in v]
            for k, v in rod["residuals"][rod["burnup"].index(u)].items()
        }
        f = r.added(f, scaled)
    return r.k_inf(f)


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
        s = c.compare(t, branch_table(), nominal, branch, switch, SUMMARY, True)
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

    def test_the_switched_core_reports_the_weight_and_the_residual_effect(self):
        t, nominal, branch, switch, w = self.cores()
        s = c.compare(t, branch_table(), nominal, branch, switch, SUMMARY, False)
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
            c.compare(t, branch_table(), nominal, other, switch, SUMMARY, True)
        # R の点で、状態の減速材の温度がテーブルの格子点でない
        off = dict(
            branch, states=dict(STATES, branch=dict(STATES["branch"], moderator=590.0))
        )
        with self.assertRaises(SystemExit):
            c.compare(t, branch_table(), nominal, off, switch, SUMMARY, True)
        # 第 1 項のテーブルが R の点を解いていない
        short = branch_table()
        short["axes"][0]["points"] = [0.0, 15.0]
        with self.assertRaises(SystemExit):
            c.compare(t, short, nominal, branch, switch, SUMMARY, True)


if __name__ == "__main__":
    unittest.main()
