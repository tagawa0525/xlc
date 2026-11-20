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

    def test_histories_with_other_burnup_points_are_errors(self):
        nominal = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        varied = history([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]])
        varied["points"][1]["burnup"] = 3.0
        with self.assertRaises(SystemExit):
            inputs.differing_nuclides(nominal, varied, [2.0])

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


# 原子量（g/mol）。重量割合と密度を確かめるだけなので 4 桁
MASS = {"U": 238.03, "Gd": 157.25, "O": 15.999}


def element(name):
    return name.split("-")[0]


def oxide_masses(nuclides):
    """UO2 と Gd2O3 の質量（数密度 × 原子量の和。酸素は化学量論で分ける）"""
    u = sum(v for n, v in nuclides.items() if element(n) == "U")
    gd = sum(v for n, v in nuclides.items() if element(n) == "Gd")
    return (
        u * (MASS["U"] + 2 * MASS["O"]),
        gd * (MASS["Gd"] + 1.5 * MASS["O"]),
    )


class CheckCell(unittest.TestCase):
    """確認のセル（ADR-0270 の確定に足りない測定の (2)）: Gd のピンの Gd を 2 倍にした gd"""

    def materials(self, fuel):
        return tomllib.loads(
            inputs.depletion_input(fuel, inputs.NOMINAL_DENSITY, "nominal")
        )["material"]

    def test_the_gadolinia_doubles_to_ten_weight_percent_at_the_same_pellet_mass(self):
        base, check = (
            self.materials("gd")["gd0"]["nuclides"],
            self.materials("gd10")["gd0"]["nuclides"],
        )
        for n, v in base.items():
            if element(n) == "Gd":
                self.assertAlmostEqual(check[n] / v, 2.0, places=4)
        uo2, gd2o3 = oxide_masses(base)
        self.assertAlmostEqual(gd2o3 / (uo2 + gd2o3), 0.05, places=3)
        uo2, gd2o3 = oxide_masses(check)
        self.assertAlmostEqual(gd2o3 / (uo2 + gd2o3), 0.10, places=3)
        self.assertAlmostEqual(
            sum(oxide_masses(check)) / sum(oxide_masses(base)), 1.0, places=4
        )
        # 酸素は UO2 と Gd2O3 の化学量論
        u = sum(v for n, v in check.items() if element(n) == "U")
        gd = sum(v for n, v in check.items() if element(n) == "Gd")
        self.assertAlmostEqual(check["O-16"] / (2 * u + 1.5 * gd), 1.0, places=4)

    def test_the_check_cell_is_the_gd_cell_but_the_gd_pin(self):
        base, check = self.materials("gd"), self.materials("gd10")
        self.assertEqual(base.keys(), check.keys())
        for m in base:
            if m in inputs.GD_MATERIALS:
                self.assertNotEqual(base[m], check[m])
            else:
                self.assertEqual(base[m], check[m])
        a = tomllib.loads(inputs.depletion_input("gd", inputs.NOMINAL_DENSITY, "n"))
        b = tomllib.loads(inputs.depletion_input("gd10", inputs.NOMINAL_DENSITY, "n"))
        del a["material"], b["material"]
        self.assertEqual(a, b)

    def test_the_check_cell_is_not_a_ranked_fuel(self):
        self.assertNotIn("gd10", inputs.FUELS)
        self.assertIn("gd10", inputs.CHECK_FUELS)


if __name__ == "__main__":
    unittest.main()
