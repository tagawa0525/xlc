"""p9_timing_inputs.py（Problem 9 の形の全炉心の燃焼の時間を測るケースを書く）。

実行: python3 -m unittest discover -s examples/depletion/p9-timing -p 'test_p9_timing_inputs.py'
"""

import sys
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p9_timing_inputs as inputs


class Inputs(unittest.TestCase):
    def test_the_committed_files_are_the_written_ones(self):
        for name, text in inputs.files().items():
            self.assertEqual((inputs.HERE / name).read_text(), text, name)

    def test_lattice_layers_take_the_2A_table_and_the_others_stay_p5(self):
        p5 = tomllib.loads(
            (inputs.HERE / inputs.P5 / "core-hzp-assemblies.toml").read_text()
        )
        ours = tomllib.loads(inputs.assemblies())
        self.assertEqual(ours["assembly"], p5["assembly"])
        self.assertEqual(ours["design"].keys(), p5["design"].keys())
        for name, design in p5["design"].items():
            self.assertEqual(ours["design"][name]["bounds"], design["bounds"])
            for theirs, layer in zip(
                design["layers"], ours["design"][name]["layers"], strict=True
            ):
                if theirs.startswith("lattice/"):
                    self.assertEqual(layer, inputs.TABLE_2A)
                else:
                    self.assertEqual(layer, str(inputs.P5 / theirs))

    def test_the_cycle_keeps_the_p5_plant_and_loading(self):
        p5 = tomllib.loads(
            (inputs.HERE / inputs.P5 / "core-hzp-cycle.toml").read_text()
        )
        ours = tomllib.loads(inputs.cycle())
        self.assertEqual(ours["loading"], p5["loading"])
        self.assertEqual(ours["plant"], str(inputs.P5 / p5["plant"]))

    def test_the_first_step_gives_the_boron_and_the_second_searches_within_both_tables(
        self,
    ):
        case = tomllib.loads(inputs.case(1.05))
        first, second = (case["state"][name] for name, *_ in inputs.STEPS)
        self.assertEqual(first["boron"], inputs.GIVEN_BORON)
        self.assertEqual(second["search"]["boron"], [1100.0, 1300.0])
        self.assertEqual(second["search"]["k"], 1.05)
        self.assertEqual(
            [s["state"] for s in case["depletion"]["steps"]], ["case-2", "case-3"]
        )
        self.assertEqual(case["depletion"]["micro"]["xenon"], "equilibrium")

    def test_without_a_target_both_steps_give_the_boron(self):
        case = tomllib.loads(inputs.case(None))
        for name, *_ in inputs.STEPS:
            self.assertEqual(case["state"][name]["boron"], inputs.GIVEN_BORON)
            self.assertNotIn("search", case["state"][name])

    def test_states_take_the_bank_D_of_table_P9_2_and_the_moderator_of_the_table(self):
        case = tomllib.loads(inputs.case(1.05))
        self.assertEqual(case["state"]["case-2"]["banks"]["D"], 192)
        self.assertEqual(case["state"]["case-3"]["banks"]["D"], 219)
        for state in case["state"].values():
            self.assertEqual(state["moderator"], 565)
            self.assertEqual(
                {b: v for b, v in state["banks"].items() if b != "D"},
                dict.fromkeys(["A", "B", "C", "SA", "SB", "SC", "SD"], 230),
            )


if __name__ == "__main__":
    unittest.main()
