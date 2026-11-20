"""read_survey.py の要約（版ごとの `xs check` の報告を、燃焼チェーンの核種とその外に分けて数える）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_read_survey.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import read_survey


def nuclide(name, z, a, state=0):
    return {"name": name, "z": z, "a": a, "state": state}


CHAIN = {
    "format": "xlc-chain",
    "nuclides": [
        nuclide("U-235", 92, 235),
        nuclide("Am-242m", 95, 242, 1),
        nuclide("Xe-135", 54, 135),
    ],
}


def evaluation(file, error=None):
    e = {"endf": f"../../references/endf/ENDF-B-VII.0_neutrons/{file}"}
    if error is not None:
        e["error"] = error
    return e


REPORT = {
    "evaluation": [
        evaluation("n-092_U_235.endf"),
        evaluation("n-095_Am_242m1.endf", "line 12: LF=11 is not supported"),
        evaluation("n-023_V_000.endf", "line 155: NIS=2 isotopes are not supported"),
        evaluation("n-030_Zn_000.endf"),
    ]
}


class NuclideKeyTest(unittest.TestCase):
    def test_reads_the_element_mass_and_state_from_the_file_name(self):
        self.assertEqual(read_survey.nuclide_key("a/n-095_Am_242m1.endf"), (95, 242, 1))
        self.assertEqual(read_survey.nuclide_key("n-092_U_235.endf"), (92, 235, 0))

    def test_a_natural_element_has_mass_zero(self):
        self.assertEqual(read_survey.nuclide_key("n-023_V_000.endf"), (23, 0, 0))

    def test_rejects_a_name_of_another_form(self):
        with self.assertRaises(SystemExit):
            read_survey.nuclide_key("U235.endf")


class ReasonTest(unittest.TestCase):
    def test_drops_only_the_leading_line_number(self):
        self.assertEqual(
            read_survey.reason("line 12: LF=11 is not supported"),
            "LF=11 is not supported",
        )
        self.assertEqual(read_survey.reason("600 K: line 3: x"), "600 K: line 3: x")


class SummarizeTest(unittest.TestCase):
    def setUp(self):
        self.summary = read_survey.summarize([("VII.0", REPORT, CHAIN)])
        self.version = self.summary["versions"]["VII.0"]

    def test_splits_the_files_by_the_chain(self):
        self.assertEqual(
            self.version["in_chain"], {"evaluations": 2, "readable": 1, "stopped": 1}
        )
        self.assertEqual(
            self.version["outside_chain"],
            {"evaluations": 2, "readable": 1, "stopped": 1},
        )

    def test_counts_the_chain_nuclides_without_a_file(self):
        self.assertEqual(self.version["chain_nuclides"], 3)
        self.assertEqual(self.version["chain_without_file"], 1)

    def test_lists_each_stopped_file_with_its_chain_name(self):
        self.assertEqual(
            self.version["stopped"],
            [
                {
                    "file": "n-095_Am_242m1.endf",
                    "nuclide": "Am-242m",
                    "reason": "LF=11 is not supported",
                },
                {
                    "file": "n-023_V_000.endf",
                    "nuclide": None,
                    "reason": "NIS=2 isotopes are not supported",
                },
            ],
        )

    def test_counts_the_files_by_reason_and_version(self):
        self.assertEqual(
            self.summary["reasons"],
            {
                "LF=11 is not supported": {"VII.0": ["Am-242m"]},
                "NIS=2 isotopes are not supported": {"VII.0": ["n-023_V_000.endf"]},
            },
        )

    def test_gives_the_watched_nuclides_their_result_per_version(self):
        watched = {w["nuclide"]: w["versions"] for w in self.summary["watched"]}
        self.assertEqual(watched["Am-242m"], {"VII.0": "LF=11 is not supported"})
        # ファイルが無ければ None
        self.assertEqual(watched["U-237"], {"VII.0": None})

    def test_rejects_a_report_that_lists_a_file_twice(self):
        report = {"evaluation": REPORT["evaluation"] + [evaluation("n-092_U_235.endf")]}
        with self.assertRaises(SystemExit):
            read_survey.summarize([("VII.0", report, CHAIN)])


if __name__ == "__main__":
    unittest.main()
