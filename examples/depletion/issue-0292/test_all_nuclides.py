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
        trial_version=True,
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

    def test_a_base_nuclide_without_a_chain_evaluation_stops_in_its_own_version(self):
        # 試行の入力と同じ版（VIII.0）では、試行の核種はどれもチェーンの評価のはず
        base = tomllib.loads(BASE)["nuclide"]
        base["C-0"] = {"endf": "n-006_C_000.endf", "temperatures": [565]}
        with self.assertRaises(SystemExit):
            all_nuclides.library_entries(
                base,
                survey({"n-090_Th_232.endf": "line 9: x"}),
                CHAIN,
                EXCLUDED,
                trial_version=True,
            )

    def test_a_base_nuclide_of_another_file_stops_in_its_own_version(self):
        base = tomllib.loads(BASE)["nuclide"]
        base["U-235"] = {**base["U-235"], "endf": "n-092_U_235b.endf"}
        with self.assertRaises(SystemExit):
            all_nuclides.library_entries(
                base,
                survey({"n-090_Th_232.endf": "line 9: x"}),
                CHAIN,
                EXCLUDED,
                trial_version=True,
            )

    def test_a_base_nuclide_without_a_chain_evaluation_is_left_out_of_another_version(
        self,
    ):
        # 試行の入力（VIII.0）の核種の評価が無い版がある。核種の組は版のチェーンと評価で決める
        base = tomllib.loads(BASE)["nuclide"]
        base["C-0"] = {"endf": "n-006_C_000.endf", "temperatures": [565]}
        e = all_nuclides.library_entries(
            base,
            survey({"n-090_Th_232.endf": "line 9: x"}),
            CHAIN,
            EXCLUDED,
            trial_version=False,
        )
        self.assertNotIn("C-0", e)
        self.assertEqual(all_nuclides.missing_from_version(base, e), ["C-0"])

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


def viii0():
    return all_nuclides.VERSIONS["endfb80"]


class Depletion(unittest.TestCase):
    def test_the_depletion_input_keeps_the_trial_input_but_its_head(self):
        trial = '# 見出し\nlibrary = "../library/depletion.msgpack"\npitch = 1.26\n'
        text = all_nuclides.depletion_input(trial, viii0())
        self.assertNotEqual(text.splitlines()[0], "# 見出し")
        self.assertEqual(text.splitlines()[1:], trial.splitlines()[1:])

    def test_the_check_cell_input_is_written_like_the_trial_inputs(self):
        # 確認のセル（gd10）は試行の入力を持たないので inputs.py の燃焼の入力から書く。同じ道で gd を書くと、試行の
        # gd の入力から書いたものと見出しのほかは同じ
        trial = all_nuclides.depletion_input(
            (all_nuclides.HERE / "trial" / "gd" / "nominal.toml").read_text(), viii0()
        )
        self.assertEqual(
            all_nuclides.check_depletion_input("gd").splitlines()[1:],
            trial.splitlines()[1:],
        )
        text = all_nuclides.check_depletion_input("gd10")
        self.assertIn("gd10", text.splitlines()[0])
        self.assertEqual(
            text.splitlines()[1], 'library = "../library/depletion.msgpack"'
        )

    def test_a_trial_input_of_another_library_stops(self):
        with self.assertRaises(SystemExit):
            all_nuclides.depletion_input(
                '# 見出し\nlibrary = "../../p2/depletion.msgpack"\n', viii0()
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


# 天然の C の評価だけを持ち、Am-242m と Th-232 の評価が無い版（VII.0 の形）
NATURAL_CHAIN = {
    "format": "xlc-chain",
    "nuclides": CHAIN["nuclides"]
    + [
        {"name": n, "z": 6, "a": a, "state": 0} for n, a in [("C-12", 12), ("C-13", 13)]
    ],
}


def other_version(errors=None):
    files = [
        "n-000_n_001.endf",
        "n-001_H_001.endf",
        "n-002_He_004.endf",
        "n-006_C_000.endf",
        "n-054_Xe_131.endf",
        "n-092_U_235.endf",
    ]
    errors = errors or {}
    return {
        "evaluation": [
            {"endf": f"../e7/{f}", **({"error": errors[f]} if f in errors else {})}
            for f in files
        ]
    }


def version_entries(**kwargs):
    return all_nuclides.library_entries(
        tomllib.loads(BASE)["nuclide"],
        kwargs.pop("survey", other_version()),
        NATURAL_CHAIN,
        {},
        natural=["C-12", "C-13"],
        trial_version=False,
        **kwargs,
    )


class Versions(unittest.TestCase):
    def test_a_natural_evaluation_goes_to_each_listed_isotope(self):
        e = all_nuclides.chain_evaluations(
            other_version(), NATURAL_CHAIN, ["C-12", "C-13"]
        )
        self.assertEqual(e["C-12"], ("n-006_C_000.endf", None))
        self.assertEqual(e["C-13"], ("n-006_C_000.endf", None))
        self.assertNotIn("n-001", e)

    def test_a_natural_evaluation_without_a_listed_isotope_stops(self):
        with self.assertRaises(SystemExit):
            all_nuclides.chain_evaluations(other_version(), NATURAL_CHAIN, [])

    def test_the_version_of_the_trial_input_must_be_given(self):
        # 試行の入力と同じ版かで、試行の核種の食い違いを止めるか入れないかが変わるので、既定を置かない
        with self.assertRaises(TypeError):
            all_nuclides.library_entries(
                tomllib.loads(BASE)["nuclide"], other_version(), NATURAL_CHAIN, {}
            )

    def test_a_listed_isotope_outside_the_chain_stops_as_such(self):
        with self.assertRaisesRegex(SystemExit, "not in the chain"):
            all_nuclides.chain_evaluations(
                other_version(), NATURAL_CHAIN, ["C-12", "C-13", "C-14"]
            )

    def test_a_listed_isotope_without_a_natural_evaluation_stops(self):
        # ADR-0262: 天然の元素の評価の無い元素の核種を挙げたらエラー（黙って核種を落とさない）
        s = other_version()
        s["evaluation"] = [e for e in s["evaluation"] if "_C_000" not in e["endf"]]
        with self.assertRaises(SystemExit):
            all_nuclides.chain_evaluations(s, NATURAL_CHAIN, ["C-12", "C-13"])

    def test_the_version_decides_the_nuclides_and_their_files(self):
        e = version_entries()
        # 試行の入力の核種のうち版に評価のある核種を試行の順に、続いて版のほかの評価を走査の順に
        self.assertEqual(list(e), ["U-235", "H-1", "Xe-131", "He-4", "C-12", "C-13"])
        self.assertEqual(e["C-13"], {"endf": "n-006_C_000.endf", "temperatures": [565]})
        self.assertNotIn("activation", e["Xe-131"])

    def test_the_base_nuclides_without_an_evaluation_of_the_version_are_listed(self):
        base = tomllib.loads(BASE)["nuclide"]
        self.assertEqual(
            all_nuclides.missing_from_version(base, version_entries()), ["Th-232"]
        )

    def test_the_thermal_scatterer_takes_the_temperatures_of_the_version(self):
        h = version_entries(thermal_temperatures=[550, 600])["H-1"]
        self.assertEqual(h["temperatures"], [550, 600])
        self.assertEqual(h["thermal"], "tsl-HinH2O.endf")
        self.assertEqual(h["free-gas-upper"], 9118)

    def test_a_reading_error_that_a_correction_declares_is_not_an_exclusion(self):
        s = other_version({"n-054_Xe_131.endf": "line 1193: last NBT must equal NP"})
        with self.assertRaises(SystemExit):
            version_entries(survey=s)
        e = version_entries(survey=s, resolved={"Xe-131"})
        self.assertIn("Xe-131", e)

    def test_the_depletion_input_of_a_version_reads_its_chain(self):
        trial = (all_nuclides.HERE / "trial" / "u" / "nominal.toml").read_text()
        text = all_nuclides.depletion_input(trial, all_nuclides.VERSIONS["endfb70"])
        data = tomllib.loads(text)
        self.assertEqual(data["depletion"]["chain"], "../../../endfb70/chain.msgpack")
        self.assertEqual(data["library"], "../library/depletion.msgpack")
        self.assertIn("ENDF/B-VII.0", text.splitlines()[0])
        ours = all_nuclides.depletion_input(trial, all_nuclides.VERSIONS["endfb80"])
        self.assertEqual(
            tomllib.loads(ours)["depletion"]["chain"], "../../../chain.msgpack"
        )
        del data["depletion"]["chain"]
        theirs = tomllib.loads(ours)
        del theirs["depletion"]["chain"]
        self.assertEqual(data, theirs)

    def test_a_cell_nuclide_outside_the_library_stops(self):
        text = "[material.fuel]\nnuclides = { U-235 = 1e-3, Am-242m = 1e-9 }\n"
        all_nuclides.check_cell_nuclides(text, {"U-235": {}, "Am-242m": {}})
        with self.assertRaises(SystemExit):
            all_nuclides.check_cell_nuclides(text, {"U-235": {}})


class VersionFiles(unittest.TestCase):
    """版ごとの設定が、その版の走査、チェーン、宣言、熱散乱のファイルを指すこと"""

    def test_each_version_reads_its_own_files(self):
        names = {
            "endfb70": "VII.0",
            "endfb71": "VII.1",
            "endfb80": "VIII.0",
            "endfb81": "VIII.1",
        }
        self.assertEqual(sorted(all_nuclides.VERSIONS), sorted(names))
        for key, version in all_nuclides.VERSIONS.items():
            with self.subTest(key):
                self.assertEqual(version["label"], f"ENDF/B-{names[key]}")
                neutrons = f"ENDF-B-{names[key]}_neutrons"
                chain = tomllib.loads(version["chain_input"].read_text())
                self.assertEqual(Path(chain["neutron-directory"]).name, neutrons)
                files = {
                    Path(e["endf"]).parent.name
                    for e in tomllib.loads(version["survey"].read_text())["evaluation"]
                }
                self.assertEqual(files, {neutrons})
                self.assertEqual(version["declarations"].parent.name, names[key])
                self.assertTrue(version["chain"].name == "chain.msgpack")
                self.assertTrue(
                    version["chain"].parent == version["chain_input"].parent
                )

    def test_only_viii0_is_the_version_of_the_trial_input(self):
        for key, version in all_nuclides.VERSIONS.items():
            self.assertEqual(version["key"], key)
            self.assertEqual(version["trial_version"], key == "endfb80")

    def test_the_head_needs_one_thermal_scatterer(self):
        for entries in (
            {},
            {
                "H-1": {"thermal": "t", "temperatures": [550]},
                "H-2": {"thermal": "t", "temperatures": [550]},
            },
        ):
            with self.assertRaises(SystemExit):
                all_nuclides.library_head(
                    all_nuclides.VERSIONS["endfb81"], ["なし"], [], entries
                )

    def test_the_head_counts_the_thermal_points_of_the_entries(self):
        entries = {"H-1": {"thermal": "t", "temperatures": [550, 573.6, 600]}}
        head = all_nuclides.library_head(
            all_nuclides.VERSIONS["endfb81"], ["なし"], [], entries
        )
        self.assertIn("3 点", head[2])

    def test_the_vii_versions_drop_the_thermal_temperatures_they_lack(self):
        # tsl-HinH2O.endf に 573.6 K と 623.6 K が無い（docs/environment.md）
        for key in ("endfb70", "endfb71"):
            self.assertEqual(
                all_nuclides.VERSIONS[key]["thermal_temperatures"], [550, 600]
            )
        for key in ("endfb80", "endfb81"):
            self.assertIsNone(all_nuclides.VERSIONS[key]["thermal_temperatures"])

    def test_the_viii0_library_input_is_written_as_before(self):
        # VIII.0 は 70048a610 の試行（Slurm 1009 のライブラリ）の入力のまま
        version = all_nuclides.VERSIONS["endfb80"]
        self.assertEqual(version["directory"], "all")
        text, _ = all_nuclides.library_text(version)
        self.assertEqual(
            text, (all_nuclides.HERE / "all" / "library" / "xs.toml").read_text()
        )


if __name__ == "__main__":
    unittest.main()
