"""inputs.py（#0293 のピンセルの入力）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_inputs.py'
"""

import sys
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inputs


class Doubled(unittest.TestCase):
    def test_midpoints_are_inserted_before_every_point_including_the_first(self):
        self.assertEqual(inputs.doubled([0.1, 0.5, 1]), [0.05, 0.1, 0.3, 0.5, 0.75, 1])

    def test_the_production_points_double_and_keep_every_original_point(self):
        points = inputs.doubled(inputs.BURNUP_40)
        self.assertEqual(len(points), 80)
        self.assertTrue(set(inputs.BURNUP_40) <= set(points))
        self.assertEqual(points, sorted(points))


class Inputs(unittest.TestCase):
    def test_the_two_depletions_differ_only_in_the_burnup_points(self):
        for fuel in inputs.FUELS:
            coarse = tomllib.loads(inputs.depletion_input(fuel, 40))
            fine = tomllib.loads(inputs.depletion_input(fuel, 80))
            self.assertEqual(len(fine["depletion"].pop("burnup")), 80)
            self.assertEqual(len(coarse["depletion"].pop("burnup")), 40)
            fine["depletion"].pop("output")
            coarse["depletion"].pop("output")
            self.assertEqual(coarse, fine)

    def test_the_two_libraries_differ_only_in_the_h1_temperatures(self):
        a = tomllib.loads(inputs.library_input("a"))
        b = tomllib.loads(inputs.library_input("b"))
        self.assertEqual(a["nuclide"].pop("H-1")["temperatures"], [550, 573.6, 600])
        self.assertEqual(b["nuclide"].pop("H-1")["temperatures"], [550, 600])
        a.pop("output")
        b.pop("output")
        self.assertEqual(a, b)
        self.assertIn(inputs.SAB_TEMPERATURE, a["nuclide"]["O-16"]["temperatures"])

    def test_the_sab_branches_put_only_the_water_at_the_evaluated_temperature(self):
        for variant in ("a", "b"):
            text = inputs.branch_input(
                "u",
                "",
                f"out/depletion-{variant}.msgpack",
                inputs.SAB_TEMPERATURE,
                "",
                [],
            )
            branches = tomllib.loads(text)["branches"]
            self.assertEqual(branches["moderator"]["temperature"], [573.6])
            self.assertEqual(branches["moderator"]["materials"], ["water"])


class Links(unittest.TestCase):
    def test_library_b_takes_the_rebuilt_light_nuclides_of_library_a(self):
        self.assertIsNone(inputs.link_target("b", "H-1"))
        self.assertEqual(
            inputs.link_target("b", "O-16"),
            inputs.HERE / "sab" / "nuclides-a" / "O-16.msgpack",
        )
        self.assertIsNone(inputs.link_target("a", "O-16"))
        self.assertEqual(
            inputs.link_target("a", "U-235"),
            inputs.HERE.parent / "xs-nuclides" / "U-235.msgpack",
        )


if __name__ == "__main__":
    unittest.main()
