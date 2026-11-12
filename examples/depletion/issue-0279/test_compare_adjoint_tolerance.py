"""compare_adjoint_tolerance.py の比較（ADR-0247 の段 4。#0279）。"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_adjoint_tolerance as cat
import compare_few_group as cfg
from test_compare_few_group import NAMES, VALUES, analytic_k


def report(tolerance: float, derivatives: dict[str, list[float]]) -> dict:
    return {
        "k": analytic_k(),
        "few_group": {
            "responses": NAMES,
            "values": [VALUES[n] for n in NAMES],
            "chi": [1.0, 0.0],
            "tolerance": tolerance,
            "sweeps": [10] * len(NAMES),
            "seconds": {"generalized": 1.0},
            "nuclides": [
                {"name": name, "relative": relative}
                for name, relative in derivatives.items()
            ],
        },
    }


ZERO = [0.0] * len(NAMES)


def shifted(i: int, by: float) -> list[float]:
    relative = list(ZERO)
    relative[i] = by
    return relative


class Differences(unittest.TestCase):
    def test_largest_difference_is_taken_relative_to_the_constant(self):
        coarse = report(1e-7, {"U-238": shifted(0, 1e-8), "Pu-239": shifted(7, 1e-8)})
        fine = report(1e-10, {"U-238": ZERO, "Pu-239": ZERO})
        largest = cat.differences(coarse, fine)["largest"]
        # 同じ差なら、total-1 の 0.5 より値の小さい scatter-1-2 の 0.02 で大きい
        self.assertEqual(largest["name"], "Pu-239")
        self.assertEqual(largest["response"], "scatter-1-2")
        self.assertAlmostEqual(largest["relative"], 1e-8 / 0.02)

    def test_a_derivative_of_a_zero_constant_is_refused(self):
        # VALUES の scatter-2-1 は 0 で、定数に対する差が決まらない
        coarse = report(1e-7, {"U-238": shifted(NAMES.index("scatter-2-1"), 1e-8)})
        fine = report(1e-10, {"U-238": ZERO})
        with self.assertRaises(SystemExit):
            cat.differences(coarse, fine)

    def test_reactivity_of_a_nuclide_sums_the_gradient_times_the_difference(self):
        i = NAMES.index("nu-fission-2")
        j = NAMES.index("total-2")
        coarse_u = list(ZERO)
        coarse_u[i], coarse_u[j] = 1e-6, 1e-6
        coarse = report(1e-7, {"U-238": coarse_u})
        fine = report(1e-10, {"U-238": ZERO})
        (row,) = cat.differences(coarse, fine)["nuclides"]
        g = cfg.k_gradient(VALUES, [1.0, 0.0])
        scale = 1e5 * cfg.FRACTION / analytic_k() ** 2
        self.assertAlmostEqual(
            row["pcm"] / (scale * (g["nu-fission-2"] + g["total-2"]) * 1e-6), 1.0
        )
        self.assertAlmostEqual(
            row["absolute"]
            / (scale * (abs(g["nu-fission-2"]) + abs(g["total-2"])) * 1e-6),
            1.0,
        )

    def test_results_of_different_states_are_refused(self):
        coarse = report(1e-7, {"U-238": ZERO})
        fine = report(1e-10, {"U-238": ZERO})
        fine["k"] *= 1.0 + 1e-9
        with self.assertRaises(SystemExit):
            cat.differences(coarse, fine)

    def test_results_of_different_nuclides_are_refused(self):
        coarse = report(1e-7, {"U-238": ZERO})
        fine = report(1e-10, {"Pu-239": ZERO})
        with self.assertRaises(SystemExit):
            cat.differences(coarse, fine)


if __name__ == "__main__":
    unittest.main()
