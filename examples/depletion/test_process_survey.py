"""process_survey.py（版ごとの中性子の評価のすべてを `xs nuclide` に通した結果の入力、報告、要約）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_process_survey.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_survey


def nuclide(name, z, a, state=0):
    return {"name": name, "z": z, "a": a, "state": state}


CHAIN = {
    "format": "xlc-chain",
    "nuclides": [
        nuclide("H-1", 1, 1),
        nuclide("C-12", 6, 12),
        nuclide("C-13", 6, 13),
        nuclide("U-235", 92, 235),
        nuclide("Am-242m", 95, 242, 1),
        nuclide("Rb-87", 37, 87),
    ],
}

DIRECTORY = "../../references/endf/ENDF-B-VII.0_neutrons"


def evaluation(file, error=None):
    e = {"endf": f"{DIRECTORY}/{file}"}
    if error is not None:
        e["error"] = error
    return e


SURVEY = {
    "evaluation": [
        evaluation("n-000_n_001.endf"),
        evaluation("n-001_H_001.endf"),
        evaluation("n-006_C_000.endf"),
        evaluation("n-037_Rb_087.endf", "line 1193: last NBT must equal NP"),
        evaluation("n-092_U_235.endf"),
        evaluation("n-095_Am_242m1.endf"),
        evaluation("n-099_Es_254.endf"),
    ]
}

H1 = {
    "endf": "n-001_H_001.endf",
    "temperatures": [550, 600],
    "thermal": "tsl-HinH2O.endf",
    "sab-upper": "evaluated",
}
BASE = {"H-1": H1, "U-235": {"endf": "n-092_U_235.endf", "temperatures": [565, 900]}}


class InputEntriesTest(unittest.TestCase):
    def setUp(self):
        self.entries, self.left_out = process_survey.input_entries(
            SURVEY, CHAIN, ["C-12", "C-13"], BASE
        )

    def test_names_the_chain_evaluations_by_the_chain_at_one_temperature(self):
        self.assertEqual(
            self.entries["U-235"],
            {"endf": "n-092_U_235.endf", "temperatures": [process_survey.TEMPERATURE]},
        )
        self.assertEqual(self.entries["Am-242m"]["endf"], "n-095_Am_242m1.endf")

    def test_keeps_the_evaluations_the_survey_stops(self):
        self.assertIn("Rb-87", self.entries)

    def test_takes_an_s_alpha_beta_nuclide_from_the_base_input_as_it_is(self):
        self.assertEqual(self.entries["H-1"], H1)

    def test_gives_a_natural_element_its_first_listed_isotope(self):
        self.assertEqual(
            self.entries["C-12"],
            {"endf": "n-006_C_000.endf", "temperatures": [process_survey.TEMPERATURE]},
        )
        self.assertNotIn("C-13", self.entries)

    def test_names_an_evaluation_outside_the_chain_by_its_file(self):
        self.assertEqual(self.entries["Es-254"]["endf"], "n-099_Es_254.endf")

    def test_leaves_out_the_free_neutron_with_the_reason(self):
        self.assertEqual([f for f, _ in self.left_out], ["n-000_n_001.endf"])
        self.assertEqual(len(self.entries), 6)

    def test_keeps_the_survey_order(self):
        self.assertEqual(
            list(self.entries), ["H-1", "C-12", "Rb-87", "U-235", "Am-242m", "Es-254"]
        )

    def test_a_natural_element_without_a_listed_isotope_is_an_error(self):
        with self.assertRaises(SystemExit):
            process_survey.input_entries(SURVEY, CHAIN, [], BASE)


class NameFromFileTest(unittest.TestCase):
    def test_writes_the_first_metastable_state_as_m(self):
        self.assertEqual(
            process_survey.name_from_file("n-095_Am_242m1.endf"), "Am-242m"
        )
        self.assertEqual(
            process_survey.name_from_file("n-052_Te_127m2.endf"), "Te-127m2"
        )
        self.assertEqual(process_survey.name_from_file("n-092_U_235.endf"), "U-235")


class ParseTimeTest(unittest.TestCase):
    def test_reads_a_finished_run(self):
        self.assertEqual(
            process_survey.parse_time("1.15 s 1.10 user 13092 kB exit 0\n"),
            {"seconds": 1.15, "user_seconds": 1.10, "max_rss_kb": 13092, "exit": 0},
        )

    def test_reads_a_failed_run(self):
        t = process_survey.parse_time(
            "Command exited with non-zero status 1\n0.48 s 0.47 user 14952 kB exit 1\n"
        )
        self.assertEqual(t["exit"], 1)
        self.assertNotIn("signal", t)

    def test_reads_a_run_killed_by_a_signal(self):
        t = process_survey.parse_time(
            "Command terminated by signal 6\n33.73 s 33.48 user 52208 kB exit 0\n"
        )
        self.assertEqual(t["signal"], 6)

    def test_reads_the_line_of_the_combined_log_with_its_name(self):
        self.assertEqual(
            process_survey.parse_time("O-16 1.40 s 1.37 user 123540 kB", named=True),
            ("O-16", {"seconds": 1.40, "user_seconds": 1.37, "max_rss_kb": 123540}),
        )

    def test_another_form_is_an_error(self):
        with self.assertRaises(SystemExit):
            process_survey.parse_time("garbage\n")


WARN_A = (
    "2026-11-20T08:28:26.695235Z  WARN nuclide{name=U-234}: xlc_xs::pointwise: 22 reactions start below "
    "their kinematic threshold; the largest shift is 1.0809e0 eV in MT17"
)
WARN_B = (
    "2026-11-20T08:28:26.695235Z  WARN nuclide{name=U-234}: xlc_xs::endf::mf2: line 572: 27 of 27 SLBW "
    "resonances with L=0 have J [0.5]"
)
CLOSE = "2026-11-20T11:39:28.979951Z  INFO nuclide{name=U-234}: xs: close time.busy=3.30ms time.idle=3.55µs"


class WarningKindTest(unittest.TestCase):
    def test_keeps_the_target_and_the_words_without_the_numbers(self):
        self.assertEqual(
            process_survey.warning_kind(WARN_A),
            "xlc_xs::pointwise: N reactions start below their kinematic threshold; the largest shift is N eV in MTN",
        )
        self.assertEqual(
            process_survey.warning_kind(WARN_B),
            "xlc_xs::endf::mf2: N of N SLBW resonances with L=N have J [N]",
        )

    def test_another_line_is_not_a_warning(self):
        self.assertIsNone(process_survey.warning_kind(CLOSE))


class EvaluationTest(unittest.TestCase):
    def test_a_finished_run_counts_its_warnings_by_kind(self):
        e = process_survey.evaluation(
            "U-234",
            "n-092_U_234.endf",
            f"{WARN_A}\n{WARN_A}\n{WARN_B}\n{CLOSE}" + "\n",
            "2.00 s 3.00 user 1000 kB exit 0\n",
        )
        self.assertEqual(e["exit"], 0)
        self.assertNotIn("error", e)
        self.assertEqual(
            e["warnings"],
            {
                process_survey.warning_kind(WARN_A): 2,
                process_survey.warning_kind(WARN_B): 1,
            },
        )
        self.assertEqual((e["seconds"], e["max_rss_kb"]), (2.0, 1000))

    def test_a_failed_run_keeps_the_error_without_the_name_and_the_line(self):
        e = process_survey.evaluation(
            "Rb-87",
            "n-037_Rb_087.endf",
            CLOSE + "\nerror: Rb-87: line 1193: last NBT must equal NP\n",
            "Command exited with non-zero status 1\n0.00 s 0.00 user 6884 kB exit 1\n",
        )
        self.assertEqual(e["exit"], 1)
        self.assertEqual(e["error"], "last NBT must equal NP")

    def test_a_run_past_the_time_limit_is_an_error(self):
        e = process_survey.evaluation(
            "U-238",
            "n-092_U_238.endf",
            CLOSE + "\n",
            "Command exited with non-zero status 124\n2.00 s 1.99 user 13372 kB exit 124\n",
        )
        self.assertEqual(e["error"], process_survey.TIMED_OUT)

    def test_a_run_killed_by_a_signal_keeps_the_last_line_of_its_log(self):
        e = process_survey.evaluation(
            "U-238",
            "n-092_U_238.endf",
            CLOSE
            + "\nmemory allocation of 16777216 bytes failed\nnote: run with `RUST_BACKTRACE=1`\n",
            "Command terminated by signal 6\n33.73 s 33.48 user 52208 kB exit 0\n",
        )
        self.assertEqual(
            e["error"], "signal 6 (memory allocation of 16777216 bytes failed)"
        )

    def test_a_failure_without_an_error_line_is_an_error(self):
        with self.assertRaises(SystemExit):
            process_survey.evaluation(
                "U-238",
                "n-092_U_238.endf",
                CLOSE + "\n",
                "0.1 s 0.1 user 1 kB exit 1\n",
            )


CLOSE_O16 = "2026-11-20T08:28:28.092625Z  INFO nuclide{name=O-16}: xs: close time.busy=1.34s time.idle=4.07µs"
# Slurm 999 の xs.log の形: 核種ごとの標準エラーと、GNU time の `<核種> %e s %U user %M kB` の行が混ざる
COMBINED = (
    f"{WARN_A}\n{CLOSE_O16}\nO-16 1.40 s 1.37 user 123540 kB\n"
    f"{CLOSE}\nU-234 2.00 s 3.00 user 1000 kB\n"
)


class CombinedLogTest(unittest.TestCase):
    def test_splits_the_finished_runs_by_nuclide(self):
        runs = process_survey.combined_runs(COMBINED)
        self.assertEqual(set(runs), {"O-16", "U-234"})
        e = runs["U-234"]
        self.assertEqual(e["exit"], 0)
        self.assertEqual(e["warnings"], {process_survey.warning_kind(WARN_A): 1})
        self.assertEqual(runs["O-16"]["max_rss_kb"], 123540)

    def test_an_error_line_is_an_error(self):
        with self.assertRaises(SystemExit):
            process_survey.combined_runs(COMBINED + "error: U-238: x\n")

    def test_a_run_without_its_close_line_is_an_error(self):
        with self.assertRaises(SystemExit):
            process_survey.combined_runs(COMBINED + "U-238 2.00 s 3.00 user 1000 kB\n")


def run(name, file, error=None, seconds=1.0, rss=1000, warnings=None):
    e = {
        "nuclide": name,
        "endf": file,
        "exit": 0 if error is None else 1,
        "seconds": seconds,
        "user_seconds": seconds,
        "max_rss_kb": rss,
        "warnings": warnings or {},
    }
    if error is not None:
        e["error"] = error
    return e


REPORT = {
    "format": "xlc-process-survey-report",
    "evaluations": [
        run("H-1", "n-001_H_001.endf", seconds=200.0, warnings={"w: a": 1}),
        run("C-12", "n-006_C_000.endf", rss=50000),
        run("Rb-87", "n-037_Rb_087.endf", "last NBT must equal NP"),
        run(
            "U-235",
            "n-092_U_235.endf",
            "negative cross section at N eV",
            warnings={"w: a": 2, "w: b": 1},
        ),
        run("Am-242m", "n-095_Am_242m1.endf"),
        run("Es-254", "n-099_Es_254.endf"),
    ],
}

CAUSES = [
    (r"last NBT", "TAB1 の補間の区切り", "endf::tab1", "#0316"),
]


class SummarizeTest(unittest.TestCase):
    def setUp(self):
        self.summary = process_survey.summarize(
            [("VII.0", REPORT, SURVEY, CHAIN)], CAUSES
        )
        self.version = self.summary["versions"]["VII.0"]

    def test_splits_the_evaluations_by_the_chain(self):
        self.assertEqual(
            self.version["in_chain"],
            {"evaluations": 4, "passed": 2, "stopped": 2, "stopped_reading": 1},
        )
        self.assertEqual(
            self.version["outside_chain"],
            {
                "evaluations": 3,
                "processed": 2,
                "passed": 2,
                "stopped": 0,
                "left_out": 1,
            },
        )

    def test_lists_each_stopped_evaluation_with_its_stage_and_cause(self):
        self.assertEqual(
            self.version["stopped"],
            [
                {
                    "nuclide": "Rb-87",
                    "file": "n-037_Rb_087.endf",
                    "chain": True,
                    "stage": "reading",
                    "reason": "last NBT must equal NP",
                    "cause": "TAB1 の補間の区切り",
                    "where": "endf::tab1",
                    "issue": "#0316",
                },
                {
                    "nuclide": "U-235",
                    "file": "n-092_U_235.endf",
                    "chain": True,
                    "stage": "processing",
                    "reason": "negative cross section at N eV",
                    "cause": None,
                    "where": None,
                    "issue": None,
                },
            ],
        )

    def test_counts_the_evaluations_with_each_warning_kind(self):
        self.assertEqual(
            self.summary["warnings"], {"w: a": {"VII.0": 2}, "w: b": {"VII.0": 1}}
        )

    def test_lists_the_slowest_and_the_largest_runs(self):
        self.assertEqual(self.version["slowest"][0]["nuclide"], "H-1")
        self.assertEqual(self.version["largest"][0]["nuclide"], "C-12")

    def test_a_survey_stop_the_run_passes_is_an_error(self):
        report = {
            **REPORT,
            "evaluations": [
                run("Rb-87", "n-037_Rb_087.endf") if e["nuclide"] == "Rb-87" else e
                for e in REPORT["evaluations"]
            ],
        }
        with self.assertRaises(SystemExit):
            process_survey.summarize([("VII.0", report, SURVEY, CHAIN)], CAUSES)

    def test_an_evaluation_missing_from_the_report_is_an_error(self):
        report = {**REPORT, "evaluations": REPORT["evaluations"][1:]}
        with self.assertRaises(SystemExit):
            process_survey.summarize([("VII.0", report, SURVEY, CHAIN)], CAUSES)


class AcrossVersionsTest(unittest.TestCase):
    def test_lists_a_nuclide_whose_time_or_memory_differs_tenfold(self):
        slow = {
            **REPORT,
            "evaluations": [
                {**e, "seconds": 20.0} if e["nuclide"] == "Am-242m" else e
                for e in REPORT["evaluations"]
            ],
        }
        summary = process_survey.summarize(
            [("VII.0", REPORT, SURVEY, CHAIN), ("VII.1", slow, SURVEY, CHAIN)], CAUSES
        )
        self.assertEqual(
            summary["across"],
            [
                {
                    "nuclide": "Am-242m",
                    "seconds": {"VII.0": 1.0, "VII.1": 20.0},
                    "max_rss_kb": {"VII.0": 1000, "VII.1": 1000},
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
