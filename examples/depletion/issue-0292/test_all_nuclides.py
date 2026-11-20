"""all_nuclides.py（#0292、データのある全核種を輸送する試行のライブラリの入力）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_all_nuclides.py'
"""

import sys
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import all_nuclides

GLOBAL = """output = "depletion.msgpack"
endf-directory = "../e"
# 書いておく理由
sigma0 = [1e10, 1]
"""
BASE = (
    "# 注釈\n"
    + GLOBAL
    + """
[nuclide.U-235]
endf = "n-092_U_235.endf"
temperatures = [560, 565, 1350]

[nuclide.H-1]
endf = "n-001_H_001.endf"
temperatures = [550, 573.6]
thermal = "tsl-HinH2O.endf"
free-gas-upper = 9118

# 放射化の核種
[nuclide]
Xe-131 = { endf = "n-054_Xe_131.endf", temperatures = [560, 565], activation = true }
Th-232 = { endf = "n-090_Th_232.endf", temperatures = [560, 565], activation = true }
"""
)
CHAIN = {
    "format": "xlc-chain",
    "nuclides": [
        {"name": n, "z": z, "a": a, "state": s}
        for n, z, a, s in [
            ("H-1", 1, 1, 0),
            ("He-4", 2, 4, 0),
            ("Xe-131", 54, 131, 0),
            ("Th-232", 90, 232, 0),
            ("U-235", 92, 235, 0),
            ("Am-242m", 95, 242, 1),
        ]
    ],
}


def survey(errors=None):
    errors = errors or {}
    files = [
        "n-000_n_001.endf",
        "n-001_H_001.endf",
        "n-002_He_004.endf",
        "n-054_Xe_131.endf",
        "n-090_Th_232.endf",
        "n-092_U_235.endf",
        "n-095_Am_242m1.endf",
    ]
    out = []
    for f in files:
        e = {"endf": f"../e/{f}"}
        if f in errors:
            e["error"] = errors[f]
        out.append(e)
    return {"evaluation": out}


EXCLUDED = {"Th-232": ("#0312", "核分裂の中性子の分布が MF6 MT18")}


def entries(errors=None, excluded=EXCLUDED):
    stopped = {"n-090_Th_232.endf": "line 9: fission needs MF1 MT452 and MF5 MT18"}
    return all_nuclides.library_entries(
        tomllib.loads(BASE)["nuclide"],
        survey(stopped if errors is None else errors),
        CHAIN,
        excluded,
    )


class Library(unittest.TestCase):
    def test_every_chain_evaluation_is_transported_at_the_cell_temperature(self):
        e = entries()
        self.assertEqual(
            list(e), ["U-235", "H-1", "Xe-131", "Th-232", "He-4", "Am-242m"]
        )
        for name in ["U-235", "Xe-131", "He-4", "Am-242m"]:
            self.assertEqual(e[name]["temperatures"], [565])
            self.assertNotIn("activation", e[name])
        self.assertEqual(e["Am-242m"]["endf"], "n-095_Am_242m1.endf")

    def test_the_thermal_scatterer_keeps_its_temperatures_and_keys(self):
        h = entries()["H-1"]
        self.assertEqual(h["temperatures"], [550, 573.6])
        self.assertEqual(h["thermal"], "tsl-HinH2O.endf")
        self.assertEqual(h["free-gas-upper"], 9118)

    def test_an_excluded_activation_nuclide_stays_activation(self):
        th = entries()["Th-232"]
        self.assertTrue(th["activation"])
        self.assertEqual(th["temperatures"], [565])

    def test_with_no_exclusions_an_activation_nuclide_of_the_base_is_transported(self):
        th = entries(errors={}, excluded={})["Th-232"]
        self.assertNotIn("activation", th)
        self.assertEqual(th["temperatures"], [565])

    def test_an_excluded_nuclide_without_an_activation_entry_is_left_out(self):
        excluded = {**EXCLUDED, "He-4": ("#0311", "理由")}
        e = entries(
            errors={
                "n-090_Th_232.endf": "line 9: fission needs MF1 MT452 and MF5 MT18",
                "n-002_He_004.endf": "line 3: neutron LAW=6 is not supported",
            },
            excluded=excluded,
        )
        self.assertNotIn("He-4", e)

    def test_the_exclusions_must_be_the_stopped_evaluations(self):
        with self.assertRaises(SystemExit):
            entries(errors={})
        with self.assertRaises(SystemExit):
            entries(
                errors={
                    "n-090_Th_232.endf": "line 9: x",
                    "n-092_U_235.endf": "line 1: y",
                }
            )

    def test_a_base_nuclide_without_a_chain_evaluation_stops(self):
        base = tomllib.loads(BASE)["nuclide"]
        base["C-0"] = {"endf": "n-006_C_000.endf", "temperatures": [565]}
        with self.assertRaises(SystemExit):
            all_nuclides.library_entries(
                base,
                survey({"n-090_Th_232.endf": "line 9: x"}),
                CHAIN,
                EXCLUDED,
            )

    def test_the_control_transports_only_what_the_base_transports(self):
        control = all_nuclides.control_entries(
            tomllib.loads(BASE)["nuclide"], entries()
        )
        self.assertNotIn("activation", control["U-235"])
        self.assertNotIn("activation", control["H-1"])
        for name in ["Xe-131", "He-4", "Am-242m", "Th-232"]:
            self.assertTrue(control[name]["activation"])
            self.assertEqual(control[name]["temperatures"], [565])


SCALE = {"unresolved-negative-background": "scale"}


class Declarations(unittest.TestCase):
    def test_declared_items_reach_the_transported_entries(self):
        e = all_nuclides.declare(entries(), {"Am-242m": SCALE, "U-235": SCALE})
        self.assertEqual(e["Am-242m"]["unresolved-negative-background"], "scale")
        self.assertEqual(e["U-235"]["unresolved-negative-background"], "scale")
        self.assertEqual(list(e), list(entries()))
        self.assertNotIn("unresolved-negative-background", e["He-4"])

    def test_activation_entries_take_the_reading_declarations_but_not_the_shielding(
        self,
    ):
        # 放射化の核種は遮蔽しない（xs は宣言を書くと止まる）ので、非分離共鳴の読み替えだけを足さない
        corrections = {"corrections": [{"line": 1, "mf": 5, "mt": 17}]}
        e = all_nuclides.declare(entries(), {"Th-232": {**SCALE, **corrections}})
        self.assertEqual(e["Th-232"]["corrections"], corrections["corrections"])
        self.assertNotIn("unresolved-negative-background", e["Th-232"])

    def test_a_declaration_outside_the_entries_is_an_error(self):
        with self.assertRaises(SystemExit):
            all_nuclides.declare(entries(), {"Es-253": SCALE})

    def test_a_declaration_conflicting_with_the_base_entry_is_an_error(self):
        e = entries()
        e["U-235"] = {**e["U-235"], "unresolved-negative-background": "other"}
        with self.assertRaises(SystemExit):
            all_nuclides.declare(e, {"U-235": SCALE})
        # 同じ値なら試行の入力の宣言をそのまま使う
        e["U-235"]["unresolved-negative-background"] = "scale"
        self.assertEqual(all_nuclides.declare(e, {"U-235": SCALE})["U-235"], e["U-235"])

    def test_the_control_leaves_out_the_shielding_of_nuclides_it_makes_activation(self):
        e = all_nuclides.declare(entries(), {"Am-242m": SCALE, "U-235": SCALE})
        control = all_nuclides.control_entries(tomllib.loads(BASE)["nuclide"], e)
        self.assertTrue(control["Am-242m"]["activation"])
        self.assertNotIn("unresolved-negative-background", control["Am-242m"])
        self.assertEqual(control["U-235"]["unresolved-negative-background"], "scale")

    def test_lists_of_tables_render_inline(self):
        e = {"Rb-87": {"endf": "rb.endf", "corrections": [{"line": 1194, "mf": 5}]}}
        text = all_nuclides.render(BASE, ["見出し"], e)
        self.assertEqual(tomllib.loads(text)["nuclide"], e)


class Render(unittest.TestCase):
    def test_the_written_input_reads_back_with_the_global_keys(self):
        e = entries()
        text = all_nuclides.render(BASE, ["見出し"], e)
        data = tomllib.loads(text)
        self.assertEqual(data["output"], "depletion.msgpack")
        self.assertEqual(data["sigma0"], [1e10, 1])
        self.assertEqual(data["nuclide"], e)
        self.assertTrue(text.startswith("# 見出し\n"))
        # 共通の項目の間の注釈は残し、ファイルの先頭の注釈は見出しに替える
        self.assertIn("# 書いておく理由\nsigma0", text)
        self.assertNotIn("# 注釈", text)

    def test_the_lists_of_nuclides_split_by_the_libraries(self):
        p2 = tomllib.loads(BASE)["nuclide"]
        trial = tomllib.loads(BASE)["nuclide"]
        del trial["Xe-131"]["activation"]
        lists = all_nuclides.nuclide_lists(entries(), trial, p2)
        self.assertEqual(lists["added-to-p2"], ["Xe-131", "He-4", "Am-242m"])
        self.assertEqual(lists["added-to-trial"], ["He-4", "Am-242m"])
        # p2 のライブラリに放射化の核種としても無い核種（89 核種の履歴の吸収の割合が無い）
        self.assertEqual(lists["unranked"], ["He-4", "Am-242m"])


class Depletion(unittest.TestCase):
    def test_the_depletion_input_keeps_the_trial_input_but_its_head(self):
        trial = '# 見出し\nlibrary = "../library/depletion.msgpack"\npitch = 1.26\n'
        text = all_nuclides.depletion_input(trial)
        self.assertNotEqual(text.splitlines()[0], "# 見出し")
        self.assertEqual(text.splitlines()[1:], trial.splitlines()[1:])

    def test_a_trial_input_of_another_library_stops(self):
        with self.assertRaises(SystemExit):
            all_nuclides.depletion_input(
                '# 見出し\nlibrary = "../../p2/depletion.msgpack"\n'
            )


class Notes(unittest.TestCase):
    def test_the_notes_say_which_exclusions_stay_activation(self):
        base = tomllib.loads(BASE)["nuclide"]
        excluded = {**EXCLUDED, "Ca-40": ("#0264", "荷電粒子の対")}
        notes = all_nuclides.exclusion_notes(base, excluded)
        self.assertIn("放射化の核種として残す", notes[0])
        self.assertIn("#0312", notes[0])
        self.assertNotIn("放射化", notes[1])

    def test_with_no_exclusions_the_note_says_none(self):
        base = tomllib.loads(BASE)["nuclide"]
        self.assertEqual(all_nuclides.exclusion_notes(base, {}), ["なし"])


class Jobs(unittest.TestCase):
    def test_the_library_job_processes_every_nuclide_and_the_control_its_activation_ones(
        self,
    ):
        e = entries()
        control = all_nuclides.control_entries(tomllib.loads(BASE)["nuclide"], e)
        jobs = all_nuclides.job_lists(e, control)
        self.assertEqual(jobs["library"], list(e))
        self.assertEqual(jobs["control"], ["Xe-131", "Th-232", "He-4", "Am-242m"])


class Duplicates(unittest.TestCase):
    def test_two_evaluations_of_one_chain_nuclide_stop(self):
        s = survey({"n-090_Th_232.endf": "line 9: x"})
        s["evaluation"].append({"endf": "../e/n-092_U_235.endf"})
        with self.assertRaises(SystemExit):
            all_nuclides.chain_evaluations(s, CHAIN)


if __name__ == "__main__":
    unittest.main()
