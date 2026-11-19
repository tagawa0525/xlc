"""untransported_absorption.py（輸送しない核種の吸収の割合）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_untransported_absorption.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import untransported_absorption as ua


def reaction(mt, *channels):
    return {"mt": mt, "channels": list(channels)}


def axis(name, points, nominal=None):
    a = {"name": name, "points": points}
    if nominal is not None:
        a["nominal"] = nominal
    return a


class Absorption(unittest.TestCase):
    def test_neutron_loss_counts_capture_and_fission_and_subtracts_the_extra_neutrons(
        self,
    ):
        nuclide = {
            "fission": 2.0,
            "reactions": [
                reaction(102, 3.0, 1.0),
                reaction(103, 0.5),
                reaction(107, 0.25),
                reaction(16, 0.1),
                reaction(17, 0.01),
                reaction(37, 0.001),
            ],
        }
        self.assertAlmostEqual(
            ua.absorption(nuclide), 2.0 + 4.0 + 0.5 + 0.25 - 0.1 - 0.02 - 0.003
        )

    def test_reactions_without_a_neutron_balance_are_errors(self):
        with self.assertRaises(SystemExit):
            ua.absorption({"fission": 0.0, "reactions": [reaction(28, 1.0)]})


class NominalStates(unittest.TestCase):
    def table(self, points):
        return {
            "format": "xlc-table",
            "cooling": {},
            "axes": [
                axis("burnup", [0.0, 10.0]),
                axis("boron", [0.0, 650.0], 650.0),
                axis("fuel", [900.0]),
            ],
            "points": points,
        }

    def test_each_burnup_takes_the_state_at_the_nominal_points_of_the_other_axes(self):
        points = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
        burnups, states = ua.nominal_states(self.table(points))
        self.assertEqual(burnups, [0.0, 10.0])
        self.assertEqual(states, [1, 3])

    def test_a_burnup_without_the_nominal_state_is_an_error(self):
        with self.assertRaises(SystemExit):
            ua.nominal_states(self.table([[0, 1, 0], [1, 0, 0]]))

    def test_an_axis_of_several_points_needs_a_nominal_point(self):
        table = self.table([[0, 0, 0], [1, 0, 0]])
        del table["axes"][1]["nominal"]
        with self.assertRaises(SystemExit):
            ua.nominal_states(table)


class TwoGroupFlux(unittest.TestCase):
    def test_flux_solves_the_removal_with_the_fission_spectrum_as_the_source(self):
        state = {
            "total": [0.5, 1.0],
            "scatter": [[0.4, 0.05], [0.0, 0.9]],
            "chi": [1.0, 0.0],
        }
        phi = ua.two_group_flux(state)
        # 高速群: (0.5 − 0.4) φ₁ = 1、熱群: (1.0 − 0.9) φ₂ = 0.05 φ₁
        self.assertAlmostEqual(phi[0], 10.0)
        self.assertAlmostEqual(phi[1], 5.0)


def history(points, nuclides):
    """燃焼の領域 1 つ（体積 2）の履歴。`points` は (時間, 燃焼度, {核種: 数密度})"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": nuclides,
        "transportable": [n for n in nuclides if n != "X"],
        "regions": [{"volume": 2.0}],
        "points": [
            {"time": t, "burnup": b, "densities": [[d.get(n, 0.0) for n in nuclides]]}
            for t, b, d in points
        ],
    }


class CapturedRatios(unittest.TestCase):
    def test_a_stable_nuclide_made_at_a_constant_rate_saturates_at_its_own_removal(
        self,
    ):
        # U-235 は σφ = 1e-2 /s で減る。X は崩壊だけの履歴で一定の生成 P = 1 /s（体積 2 で 2）。σₐ,X φ = 0.1 /s
        phi = 1e-2 / 50.0
        times = [0.0, 10.0, 20.0, 40.0]
        points = [
            (t, t, {"U-235": 1.0 * math.exp(-1e-2 * t), "X": 1.0 * t}) for t in times
        ]
        h = history(points, ["U-235", "X"])
        burnups = [0.0, 40.0]
        sigmas = {"U-235": [50.0, 50.0], "X": [0.1 / phi, 0.1 / phi]}
        ratios = ua.captured_ratios(h, [0.0, 10.0, 20.0, 40.0], sigmas, ["X"], {})
        for t in times[1:]:
            expected = (1.0 - math.exp(-0.1 * t)) / 0.1
            self.assertAlmostEqual(ratios["X"][t], expected / t, places=10)
        self.assertEqual(sorted(ratios["X"]), [10.0, 20.0, 40.0])
        del burnups

    def test_a_decaying_nuclide_keeps_its_decay_and_adds_the_capture(self):
        # λ = 0.05 /s、σₐ φ = 0.05 /s。崩壊だけの履歴は P / λ (1 − e^{−λt})
        decay, capture = 0.05, 0.05
        half_life = math.log(2.0) / decay
        phi = 1e-2 / 50.0
        times = [0.0, 5.0, 30.0]
        points = [
            (
                t,
                t,
                {"U-235": math.exp(-1e-2 * t), "X": (1 - math.exp(-decay * t)) / decay},
            )
            for t in times
        ]
        h = history(points, ["U-235", "X"])
        sigmas = {"U-235": [50.0], "X": [capture / phi]}
        ratios = ua.captured_ratios(h, [5.0, 30.0], sigmas, ["X"], {"X": half_life})
        removal = decay + capture
        for t in times[1:]:
            expected = (1 - math.exp(-removal * t)) / removal
            decayed = (1 - math.exp(-decay * t)) / decay
            self.assertAlmostEqual(ratios["X"][t], expected / decayed, places=10)


class Ranking(unittest.TestCase):
    def point(self, burnup, k_inf, shares):
        return {
            "burnup": burnup,
            "k_inf": k_inf,
            "untransported_all": {
                n: {"fraction": f, "captured_ratio": r} for n, (f, r) in shares.items()
            },
        }

    def test_nuclides_rank_by_their_largest_reactivity_over_burnup_with_the_rest_summed(
        self,
    ):
        points = [
            self.point(
                10.0, 1.0, {"A": (1e-4, 1.0), "B": (5e-4, 0.5), "C": (1e-5, 1.0)}
            ),
            self.point(
                60.0, 0.5, {"A": (4e-4, 0.5), "B": (1e-4, 1.0), "C": (2e-5, 1.0)}
            ),
        ]
        result = ua.ranking(points)
        rows = result["ranking"]
        # A: max(10, 40) = 40 pcm、B: max(25, 20) = 25、C: max(1, 4) = 4
        self.assertEqual([r["nuclide"] for r in rows], ["A", "B", "C"])
        self.assertAlmostEqual(rows[0]["max_pcm"], 40.0)
        self.assertEqual(rows[0]["burnup"], 60.0)
        self.assertAlmostEqual(rows[1]["max_pcm"], 25.0)
        self.assertEqual(rows[1]["burnup"], 10.0)
        self.assertAlmostEqual(rows[0]["rest_pcm"], 29.0)
        self.assertAlmostEqual(rows[1]["rest_pcm"], 4.0)
        self.assertAlmostEqual(rows[2]["rest_pcm"], 0.0)
        # 残りの和が 50 pcm 以下になるのは 1 つ目から、5 pcm 以下は 2 つ目から
        self.assertEqual(result["counts"], {"50.0": 1, "5.0": 2})

    def test_nuclides_xs_does_not_transport_are_marked(self):
        points = [self.point(10.0, 1.0, {"Pm-147": (1e-4, 1.0), "Xe-131": (1e-5, 1.0)})]
        rows = ua.ranking(points)["ranking"]
        self.assertEqual([r["xs_excluded"] for r in rows], [True, False])


class Prune(unittest.TestCase):
    def test_nuclides_below_the_floor_leave_the_json_with_their_count_and_sum(self):
        points = [
            {"untransported_all": {"A": {"fraction": 1.0}, "B": {"fraction": 2.0}}}
        ]
        ranked = {
            "ranking": [
                {"nuclide": "A", "max_pcm": 3.0, "rest_pcm": 0.004},
                {"nuclide": "B", "max_pcm": 0.004, "rest_pcm": 0.0},
            ],
            "counts": {},
        }
        pruned = ua.prune(points, ranked, 0.01)
        self.assertEqual(list(pruned["points"][0]["untransported_all"]), ["A"])
        self.assertEqual([r["nuclide"] for r in pruned["ranking"]], ["A"])
        self.assertEqual(pruned["pruned"], {"count": 1, "sum_pcm": 0.004})


class Category(unittest.TestCase):
    def test_tracked_nuclides_come_first_then_the_transported_ones(self):
        transportable = {"U-235", "Rh-103"}
        self.assertEqual(ua.category("U-235", transportable), "tracked")
        self.assertEqual(ua.category("Rh-103", transportable), "transported")
        self.assertEqual(ua.category("Xe-131", transportable), "untransported")


if __name__ == "__main__":
    unittest.main()
