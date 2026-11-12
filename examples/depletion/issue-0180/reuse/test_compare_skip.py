"""compare_skip.py の 1 次の見込みと 2 群の定数（#0180）。

実行: cd examples/depletion/issue-0180/reuse && python3 -m unittest test_compare_skip
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_skip import bound, estimate, homogenized, shielding_burnup, two_group


class Estimate(unittest.TestCase):
    def test_the_estimate_sums_the_effects_times_the_relative_shift_of_the_shielding(
        self,
    ):
        # 自己遮蔽の組成が状態より 2 % 多い核種は効きの 0.02 倍、少ない核種は負の向き
        total, terms = estimate(
            {"a": 100.0, "b": -50.0},
            {"a": 1.02, "b": 0.9},
            {"a": 1.0, "b": 1.0},
        )
        self.assertAlmostEqual(total, 2.0 + 5.0, places=12)
        self.assertEqual([n for n, _ in terms], ["b", "a"])

    def test_nuclides_absent_from_both_compositions_add_nothing(self):
        total, terms = estimate({"a": 100.0}, {"a": 0.0}, {"a": 0.0})
        self.assertEqual((total, terms), (0.0, []))

    def test_a_nuclide_only_in_the_shielding_is_an_error(self):
        with self.assertRaises(SystemExit):
            estimate({"a": 100.0}, {"a": 1.0}, {"a": 0.0})


class Bound(unittest.TestCase):
    def test_the_bound_sums_the_terms_without_their_cancellation(self):
        self.assertEqual(bound([("a", -3.0), ("b", 2.0)]), 5.0)
        self.assertEqual(bound([]), 0.0)


class ShieldingBurnup(unittest.TestCase):
    def test_the_points_come_from_the_history_of_the_input_that_made_the_table(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "skip.toml"
            path.write_text(
                "[history]\nburnup = [0, 1, 2]\nshielding-burnup = [0, 2]\n"
            )
            self.assertEqual(shielding_burnup(path), [0, 2])
            path.write_text("[history]\nburnup = [0, 1, 2]\n")
            with self.assertRaises(SystemExit):
                shielding_burnup(path)


class Homogenized(unittest.TestCase):
    def test_densities_are_volume_weighted_and_negative_ones_are_zero(self):
        history = {
            "nuclides": ["x"],
            "regions": [{"volume": 1.0}, {"volume": 3.0}],
            "points": [{"burnup": 2.0, "densities": [[4.0], [-1.0]]}],
        }
        self.assertEqual(homogenized(history, 2.0), {"x": 1.0})


class TwoGroup(unittest.TestCase):
    def test_constants_take_the_change_times_the_scale(self):
        state = {
            "total": [1.0, 2.0],
            "nu_fission": [0.1, 0.2],
            "scatter": [[0.5, 0.25], [0.0, 1.5]],
        }
        change = {
            "total": [10.0, 0.0],
            "nu_fission": [0.0, 1.0],
            "scatter": [[0.0, 2.0], [0.0, 0.0]],
            "fission": [9.0, 9.0],
        }
        values = two_group(state, change, 0.1)
        self.assertAlmostEqual(values["total-1"], 2.0)
        self.assertAlmostEqual(values["nu-fission-2"], 0.3)
        self.assertAlmostEqual(values["scatter-1-2"], 0.45)
        self.assertEqual(values["scatter-2-1"], 0.0)
        self.assertEqual(two_group(state)["scatter-1-1"], 0.5)


if __name__ == "__main__":
    unittest.main()
