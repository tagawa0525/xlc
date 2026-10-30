"""compare_residual.py（履歴の残り R の大きさ、R の定義の検査、状態への依存、第 1 項の再現）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_residual.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_residual as c

NUCLIDES = ["U-235", "Pu-240"]


def few(s: float) -> dict:
    """2 群の少数群定数。`s` で値を変える"""
    x = lambda v: v * (1.0 + s / 7.0)
    state = {
        "total": [x(0.55), x(1.3)],
        "absorption": [x(0.01), x(0.1)],
        "nu_fission": [x(0.007), x(0.14)],
        "fission": [x(0.0028), x(0.057)],
        "chi": [1.0, 0.0],
        "diffusion": [x(1.4), x(0.4)],
        "scatter": [[x(0.52), x(0.017)], [x(0.0015), x(1.2)]],
        "adf": [[1.0, 1.0]] * 6,
    }
    state["k_inf"] = c.k_inf(state)
    return state


def sigma(s: float, derivative: bool) -> dict:
    x = lambda v: v * (1.0 + s / 5.0)
    out = {
        "total": [x(11.0), x(130.0)],
        "nu_fission": [x(1.1), x(170.0)],
        "fission": [x(0.4), x(70.0)],
        "scatter": [[x(9.0), x(0.1)], [x(0.01), x(8.0)]],
    }
    if derivative:
        out["derivative"] = {
            "total": [x(9.0), x(110.0)],
            "nu_fission": [x(1.0), x(-150.0)],
            "fission": [x(0.35), x(-60.0)],
            "scatter": [[x(8.5), x(0.12)], [x(0.02), x(7.0)]],
        }
    return out


def axes(burnups, temperatures, densities, nominal):
    return [
        {"name": "burnup", "points": burnups},
        {"name": "boron", "points": [1300.0], "b10": 0.199, "nominal": 1300.0},
        {
            "name": "moderator",
            "points": temperatures,
            "density": densities,
            "nominal": nominal,
        },
        {"name": "fuel", "points": [565.0], "nominal": 565.0},
    ]


def residual(branch: dict, linear: dict) -> dict:
    keys = ("total", "absorption", "nu_fission", "fission", "chi", "diffusion")
    out = {k: [a - b for a, b in zip(branch[k], linear[k])] for k in keys}
    for k in ("scatter", "adf"):
        out[k] = [[a - b for a, b in zip(p, q)] for p, q in zip(branch[k], linear[k])]
    return out


class Compare(unittest.TestCase):
    def tables(self):
        """公称（燃焼度 0、10、20。減速材の公称 565 K と条件 600 K。条件は 0 と 20 で解く）、条件の分岐、逆向き、履歴"""
        points = [[0, 0, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [2, 0, 0, 0], [2, 0, 1, 0]]
        states = [few(p) for p in range(5)]
        micro = {
            "nuclides": NUCLIDES,
            "densities": [[7.0e-4, 0.0], [6.0e-4, 1.0e-5], [5.0e-4, 2.0e-5]],
            "states": [[sigma(p, False), sigma(2 * p, True)] for p in range(5)],
        }
        branch_states = [few(1.5), few(4.5)]
        branch_densities = [[7.0e-4, 0.0], [5.2e-4, 2.6e-5]]
        reverse_states = [few(0.2), few(3.3)]
        # R は分岐の定数 − 公称の条件の状態の線形の補正の後（lattice histories と同じ式）
        rs = []
        for k, (p, u) in enumerate(((1, 0), (4, 2))):
            deltas = [a - b for a, b in zip(branch_densities[k], micro["densities"][u])]
            linear = c.corrected(states[p], c.coefficients(micro["states"][p]), deltas)
            rs.append(residual(branch_states[k], linear))
        nominal = {
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

        def single(states, temperature, density):
            return {
                "format": "xlc-table",
                "axes": axes([0.0, 20.0], [temperature], [density], temperature),
                "points": [[0, 0, 0, 0], [1, 0, 0, 0]],
                "states": states,
                "micro": {
                    "nuclides": NUCLIDES,
                    "densities": branch_densities,
                    "states": [[sigma(9, False), sigma(9, False)]] * 2,
                },
            }

        branch = single(branch_states, 600.0, 0.661)
        reverse = single(reverse_states, 565.0, 0.743)
        history = {
            "format": "xlc-depletion-history",
            "points": [
                {"burnup": 0.0, "k_inf": branch_states[0]["k_inf"]},
                {"burnup": 10.0, "k_inf": 1.0},
                {
                    "burnup": 20.0,
                    "k_inf": 1.0 / (1.0 / branch_states[1]["k_inf"] + 2e-6),
                },
            ],
        }
        return nominal, branch, reverse, history

    def test_the_residual_reproduces_the_branch_and_its_effect_is_measured(self):
        nominal, branch, reverse, history = self.tables()
        s = c.compare(nominal, branch, reverse, history)
        self.assertEqual(s["format"], "xlc-compare-history-residual")
        self.assertEqual([r["burnup"] for r in s["rows"]], [0.0, 20.0])
        for k, r in enumerate(s["rows"]):
            # 公称の条件の状態 + 線形の補正 + R は分岐の定数（2 群の k の差は丸めだけ）
            self.assertAlmostEqual(r["reproduce_constants_pcm"], 0.0, places=8)
            f0 = nominal["states"][(1, 4)[k]]
            fb = branch["states"][k]
            self.assertAlmostEqual(
                r["effect_pcm"], c.rho(c.k_inf(fb)) - c.rho(c.k_inf(f0)), places=9
            )
            # 効きは線形の補正と R に分かれる
            self.assertAlmostEqual(
                r["effect_pcm"], r["linear_pcm"] + r["residual_pcm"], places=8
            )
            self.assertAlmostEqual(
                r["dependence_pcm"],
                r["residual_pcm"] - r["reverse_residual_pcm"],
                places=9,
            )
        # 第 1 項の分岐は履歴の k∞ を再現する（20 で履歴の 1/k が 2e-6 大きい = 分岐が +0.2 pcm）
        self.assertAlmostEqual(s["rows"][0]["first_term_pcm"], 0.0, places=9)
        self.assertAlmostEqual(s["rows"][1]["first_term_pcm"], 0.2, places=6)

    def test_the_reverse_residual_uses_the_nominal_state_and_the_same_formula(self):
        nominal, branch, reverse, history = self.tables()
        s = c.compare(nominal, branch, reverse, history)
        r = s["rows"][1]
        f0 = nominal["states"][3]
        deltas = [
            a - b
            for a, b in zip(
                branch["micro"]["densities"][1], nominal["micro"]["densities"][2]
            )
        ]
        linear = c.corrected(f0, c.coefficients(nominal["micro"]["states"][3]), deltas)
        expected = c.rho(c.k_inf(reverse["states"][1])) - c.rho(c.k_inf(linear))
        self.assertAlmostEqual(r["reverse_residual_pcm"], expected, places=9)
        # 条件の R の公称の状態での効きと、それを公称の状態で使った誤差（= 公称の状態での R の効き − 逆向きの R の効き）
        self.assertAlmostEqual(
            r["residual_nominal_pcm"],
            c.rho(
                c.k_inf(
                    c.added(linear, nominal["histories"][0]["rods"][0]["residuals"][1])
                )
            )
            - c.rho(c.k_inf(linear)),
            places=9,
        )
        self.assertAlmostEqual(
            r["applied_pcm"],
            r["residual_nominal_pcm"] - r["reverse_residual_pcm"],
            places=8,
        )
        with_r = c.added(linear, nominal["histories"][0]["rods"][0]["residuals"][1])
        self.assertAlmostEqual(
            r["applied_pcm"],
            c.rho(c.k_inf(with_r)) - c.rho(c.k_inf(reverse["states"][1])),
            places=9,
        )

    def test_the_whole_difference_at_the_nominal_state_is_measured_against_the_reverse(
        self,
    ):
        # 定数の差の全体 f(分岐の組成、条件 j) − f(公称の組成、条件 j) を、公称の状態の定数に足した誤差（ADR-0215 の案 S1 の R'）
        nominal, branch, reverse, history = self.tables()
        s = c.compare(nominal, branch, reverse, history)
        r = s["rows"][1]
        whole = residual(branch["states"][1], nominal["states"][4])
        expected = c.rho(c.k_inf(c.added(nominal["states"][3], whole))) - c.rho(
            c.k_inf(reverse["states"][1])
        )
        self.assertAlmostEqual(r["applied_whole_pcm"], expected, places=9)
        self.assertNotAlmostEqual(r["applied_whole_pcm"], r["applied_pcm"], places=3)

    def test_tables_that_do_not_match_are_errors(self):
        nominal, branch, reverse, history = self.tables()
        nominal["histories"] = []
        with self.assertRaises(SystemExit):
            c.compare(nominal, branch, reverse, history)
        nominal, branch, reverse, history = self.tables()
        branch["axes"][0]["points"] = [0.0, 10.0]
        with self.assertRaises(SystemExit):
            c.compare(nominal, branch, reverse, history)
        # 逆向きは同じ履歴の組成
        nominal, branch, reverse, history = self.tables()
        reverse["micro"]["densities"] = [[7.0e-4, 0.0], [5.3e-4, 2.6e-5]]
        with self.assertRaises(SystemExit):
            c.compare(nominal, branch, reverse, history)
        # 逆向きは公称の条件
        nominal, branch, reverse, history = self.tables()
        reverse["axes"] = axes([0.0, 20.0], [600.0], [0.661], 600.0)
        with self.assertRaises(SystemExit):
            c.compare(nominal, branch, reverse, history)

    def test_the_report_gives_the_ranges(self):
        nominal, branch, reverse, history = self.tables()
        text = c.report(c.compare(nominal, branch, reverse, history))
        self.assertIn("R の効き", text)
        self.assertIn("状態への依存", text)
        self.assertIn("第 1 項", text)


if __name__ == "__main__":
    unittest.main()
