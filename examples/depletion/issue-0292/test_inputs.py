"""inputs.py（#0292 の燃料ごとの入力）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_inputs.py'
"""

import sys
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import inputs


def history(densities):
    """燃焼の領域 1 つ、燃焼点 0 と 2 の履歴。`densities` は燃焼点ごとの核種の数密度"""
    return {
        "nuclides": ["O-16", "U-235", "Xe-131", "Cs-133"],
        "transportable": ["O-16", "U-235", "Cs-133"],
        "points": [
            {"burnup": b, "densities": [d]} for b, d in zip([0.0, 2.0], densities)
        ],
    }


class DifferingNuclides(unittest.TestCase):
    def test_transported_nuclides_that_differ_at_a_branch_point_are_replaced(self):
        nominal = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        varied = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.5, 3.5, 4.5]])
        # Xe-131 は違うが輸送しない。O-16 は同じ
        self.assertEqual(
            inputs.differing_nuclides(nominal, varied, [2.0]), ["U-235", "Cs-133"]
        )

    def test_differences_only_outside_the_branch_points_do_not_count(self):
        nominal = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        varied = history([[1.5, 2.5, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        self.assertEqual(inputs.differing_nuclides(nominal, varied, [2.0]), [])

    def test_branch_points_missing_from_the_histories_are_errors(self):
        h = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        with self.assertRaises(SystemExit):
            inputs.differing_nuclides(h, h, [6.0])


class Inputs(unittest.TestCase):
    def test_branch_inputs_hold_one_state_at_the_given_density(self):
        text = inputs.branch_input(
            "gd", "x", "nominal.msgpack", [0, 2], ["U-235"], inputs.NOMINAL_DENSITY
        )
        parsed = tomllib.loads(text)
        self.assertEqual(parsed["branches"]["moderator"]["density"], [0.743])
        self.assertEqual(parsed["material"]["water"]["water"]["density"], 0.743)
        self.assertEqual(parsed["history"]["micro"], ["U-235"])
        self.assertEqual(
            parsed["branches"]["fuel"]["materials"], ["fuel", *inputs.GD_MATERIALS]
        )
        self.assertEqual(parsed["library"], "../../../p2/depletion.msgpack")

    def test_depletion_inputs_burn_every_fuel_material(self):
        parsed = tomllib.loads(
            inputs.depletion_input("gd", inputs.VARIED_DENSITY, "density")
        )
        self.assertEqual(
            parsed["depletion"]["materials"], ["fuel", *inputs.GD_MATERIALS]
        )
        self.assertEqual(parsed["material"]["water"]["nuclides"]["H-1"], 4.41459e-02)
        self.assertEqual(parsed["library"], "../../p2/depletion.msgpack")


if __name__ == "__main__":
    unittest.main()
