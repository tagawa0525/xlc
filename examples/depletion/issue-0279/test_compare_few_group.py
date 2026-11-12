"""compare_few_group.py の比較（ADR-0247 の段 4。#0279）。"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_few_group as cfg

NAMES = [
    "total-1", "total-2", "nu-fission-1", "nu-fission-2", "fission-1", "fission-2",
    "scatter-1-1", "scatter-1-2", "scatter-2-1", "scatter-2-2",
]  # fmt: skip

# 上方散乱のない 2 群の定数（χ = (1, 0)）
VALUES = {
    "total-1": 0.5,
    "total-2": 1.2,
    "nu-fission-1": 0.006,
    "nu-fission-2": 0.12,
    "fission-1": 0.0025,
    "fission-2": 0.05,
    "scatter-1-1": 0.45,
    "scatter-1-2": 0.02,
    "scatter-2-1": 0.0,
    "scatter-2-2": 1.1,
}


def analytic_k() -> float:
    """上方散乱がなければ k∞ = (νΣf1 + νΣf2 Σ12 / Σa2) / Σr1（Σr1 = Σt1 − Σ11、Σa2 = Σt2 − Σ22）"""
    v = VALUES
    removal = v["total-1"] - v["scatter-1-1"]
    absorption = v["total-2"] - v["scatter-2-2"]
    return (
        v["nu-fission-1"] + v["nu-fission-2"] * v["scatter-1-2"] / absorption
    ) / removal


def steps(central: tuple[float, float]) -> dict:
    coarse, fine = central
    return {
        "central": [[coarse] * 10, [fine] * 10],
        "richardson": [fine + (fine - coarse) / 3.0] * 10,
    }


def report(
    adjoint: float, direct: tuple[float, float], shielded: dict | None = None
) -> dict:
    k = analytic_k()
    nuclide = {"name": "U-238", "relative": [adjoint] * 10, "direct": steps(direct)}
    if shielded is not None:
        nuclide["shielded"] = shielded
    return {
        "k": k,
        "epsilon": 0.1,
        "flux_tolerance": 1e-11,
        "nuclides": [
            {"name": "U-238", "density": 7e-3, "direct": {"richardson": -0.18}},
        ],
        "few_group": {
            "responses": NAMES,
            "values": [VALUES[n] for n in NAMES],
            "chi": [1.0, 0.0],
            "nuclides": [nuclide],
        },
    }


class TwoGroup(unittest.TestCase):
    def test_k_infinity_without_upscatter_is_analytic(self):
        self.assertAlmostEqual(
            cfg.two_group_k(VALUES, [1.0, 0.0]), analytic_k(), places=14
        )

    def test_gradient_of_the_thermal_production(self):
        # ∂k∞/∂νΣf2 = Σ12 / (Σa2 Σr1)
        g = cfg.k_gradient(VALUES, [1.0, 0.0])
        expected = VALUES["scatter-1-2"] / ((1.2 - 1.1) * (0.5 - 0.45))
        self.assertAlmostEqual(g["nu-fission-2"] / expected, 1.0, places=8)


class Explicit(unittest.TestCase):
    def test_agreement_within_the_step_error_passes(self):
        result = cfg.compare(report(-0.059266, (-0.05925, -0.059262)))
        self.assertTrue(all(r["ok"] for r in result["explicit"]))

    def test_difference_beyond_every_limit_fails(self):
        result = cfg.compare(report(-0.0600, (-0.059266, -0.059266)))
        self.assertFalse(any(r["ok"] for r in result["explicit"]))


class Implicit(unittest.TestCase):
    def shielded(self, first: float, direct: float) -> dict:
        return {
            "first_order": steps((first, first)),
            "direct": steps((direct, direct)),
            "k_first_order": {"central": [-0.17, -0.17], "richardson": -0.17},
            "k_direct": {"central": [-0.17, -0.17], "richardson": -0.17},
        }

    def test_implicit_part_is_the_reshielded_minus_the_fixed_direct_perturbation(self):
        result = cfg.compare(report(-0.05, (-0.05, -0.05), self.shielded(-0.03, -0.03)))
        (n,) = result["implicit"]
        self.assertTrue(all(r["ok"] for r in n["rows"]))
        for r in n["rows"]:
            self.assertAlmostEqual(r["implicit"], 0.02)
        # k の直接の摂動の差 (−0.17 + 0.18) を 1 割の変化の pcm に
        self.assertAlmostEqual(n["k_pcm"], 1e5 * 0.1 * 0.01 / analytic_k() ** 2)

    def test_first_order_away_from_the_reshielded_perturbation_fails(self):
        result = cfg.compare(report(-0.05, (-0.05, -0.05), self.shielded(-0.02, -0.03)))
        (n,) = result["implicit"]
        self.assertFalse(any(r["ok"] for r in n["rows"]))

    def test_pair_estimate_weights_the_implicit_part_by_the_density_change(self):
        implicit = [{"name": "U-238", "k_pcm": 900.0}, {"name": "Pu-240", "k_pcm": 8.0}]
        terms = cfg.pair_estimate(implicit, {"U-238": -0.005, "Pu-240": 0.03})
        self.assertEqual([t[0] for t in terms], ["U-238", "Pu-240"])
        self.assertAlmostEqual(terms[0][1], 900.0 / 0.1 * -0.005)


class Inputs(unittest.TestCase):
    def test_other_than_two_groups_is_refused(self):
        r = report(-0.05, (-0.05, -0.05))
        r["few_group"]["responses"] = r["few_group"]["responses"][:9]
        with self.assertRaises(SystemExit):
            cfg.compare(r)

    def test_a_report_without_direct_perturbations_is_refused(self):
        r = report(-0.05, (-0.05, -0.05))
        del r["epsilon"]
        with self.assertRaises(SystemExit):
            cfg.compare(r)


class Relative(unittest.TestCase):
    def test_exactly_zero_derivatives_agree(self):
        self.assertEqual(cfg.relative(0.0, 0.0), 0.0)
        self.assertEqual(cfg.relative(1e-9, 0.0), float("inf"))
        self.assertAlmostEqual(cfg.relative(1.1, 1.0), 0.1)

    def test_step_error_of_the_first_order_differences_widens_the_bound(self):
        # I1 の刻みの誤差（ε と ε/2 の差）が I3 の側より大きければ、I1 の側の見積もりで許す
        first = steps((-0.0300, -0.0306))
        direct = steps((-0.0306, -0.0306))
        shielded = {
            "first_order": first,
            "direct": direct,
            "k_first_order": {"central": [-0.17, -0.17], "richardson": -0.17},
            "k_direct": {"central": [-0.17, -0.17], "richardson": -0.17},
        }
        result = cfg.compare(report(-0.05, (-0.05, -0.05), shielded))
        (n,) = result["implicit"]
        self.assertTrue(all(r["ok"] for r in n["rows"]))


class ImplicitK(unittest.TestCase):
    def test_first_order_k_away_from_the_reshielded_k_fails(self):
        shielded = {
            "first_order": steps((-0.03, -0.03)),
            "direct": steps((-0.03, -0.03)),
            "k_first_order": {"central": [-0.16, -0.16], "richardson": -0.16},
            "k_direct": {"central": [-0.17, -0.17], "richardson": -0.17},
        }
        r = report(-0.05, (-0.05, -0.05), shielded)
        r["k_tolerance"] = 1e-12
        (n,) = cfg.compare(r)["implicit"]
        self.assertFalse(n["k_ok"])
        shielded["k_first_order"] = shielded["k_direct"]
        (n,) = cfg.compare(r)["implicit"]
        self.assertTrue(n["k_ok"])


if __name__ == "__main__":
    unittest.main()
