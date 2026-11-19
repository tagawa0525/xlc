"""background_inputs.py（同じ組成で輸送しない核種の扱いだけを変えたブランチ計算の入力）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_background_inputs.py'
"""

import sys
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import background_inputs


def history(transportable, dilute, densities):
    return {
        "format": "xlc-depletion-history",
        "nuclides": ["O-16", "U-238", "Xe-136", "Eu-155"],
        "regions": [{"volume": 0.5}],
        "transportable": transportable,
        "dilute_absorbers": dilute,
        "points": [{"burnup": 60.0, "k_inf": 0.8, "densities": [densities]}],
    }


class WithLists(unittest.TestCase):
    def test_the_composition_is_kept_and_the_lists_are_replaced(self):
        composition = history(
            ["O-16", "U-238", "Xe-136"], ["Eu-155"], [1.0, 2.0, 3.0, 4.0]
        )
        lists = history(["O-16", "U-238"], ["Xe-136", "Eu-155"], [9.0, 9.0, 9.0, 9.0])
        mixed = background_inputs.with_lists(composition, lists)
        self.assertEqual(mixed["transportable"], ["O-16", "U-238"])
        self.assertEqual(mixed["dilute_absorbers"], ["Xe-136", "Eu-155"])
        self.assertEqual(mixed["points"][0]["densities"], [[1.0, 2.0, 3.0, 4.0]])
        self.assertEqual(composition["transportable"], ["O-16", "U-238", "Xe-136"])

    def test_histories_of_other_nuclides_are_an_error(self):
        composition = history(["O-16"], [], [1.0, 2.0, 3.0, 4.0])
        lists = history(["O-16"], [], [1.0, 2.0, 3.0, 4.0])
        lists["nuclides"] = ["O-16", "U-238", "Xe-136", "Sm-149"]
        with self.assertRaisesRegex(SystemExit, "nuclides"):
            background_inputs.with_lists(composition, lists)


class CaseInput(unittest.TestCase):
    def test_the_input_solves_the_burnup_with_the_library_groups_and_every_transported_nuclide(
        self,
    ):
        text = background_inputs.case_input(
            "u",
            "form",
            "issue-0292/dilute/library-13/depletion.msgpack",
            [1e7, 4.0, 1e-5],
            ["O-16", "U-238"],
        )
        data = tomllib.loads(text)
        self.assertEqual(
            data["library"], "../../../issue-0292/dilute/library-13/depletion.msgpack"
        )
        self.assertEqual(
            data["few-group"],
            {"output": "out/form.msgpack", "bounds": [1e7, 4.0, 1e-5]},
        )
        self.assertEqual(data["history"]["file"], "out/form-history.msgpack")
        self.assertEqual(data["history"]["burnup"], [background_inputs.BURNUP])
        self.assertEqual(data["history"]["micro"], ["O-16", "U-238"])
        self.assertEqual(data["branches"]["moderator"]["density"], [0.743])
        self.assertTrue(text.startswith("# #0292"))
        self.assertIn("background_inputs.py", text.splitlines()[0])


if __name__ == "__main__":
    unittest.main()
