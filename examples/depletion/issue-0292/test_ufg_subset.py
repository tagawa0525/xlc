"""ufg_subset.py（#0292、ADR-0268 の UFG で解く核種の組）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_ufg_subset.py'
"""

import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ufg_subset


class Listed(unittest.TestCase):
    def test_the_list_drops_the_unkept_added_nuclides_and_adds_the_structure(self):
        names = ufg_subset.listed(
            ["O-16", "U-235", "Xe-131", "Cs-134"], ["Xe-131", "Cs-134"], ["Xe-131"]
        )
        self.assertEqual(names[:3], ["O-16", "U-235", "Xe-131"])
        self.assertNotIn("Cs-134", names)
        for n in ["H-1", "B-10", "Zr-90", "Sn-124", "Fe-56"]:
            self.assertIn(n, names)
        self.assertEqual(names.count("O-16"), 1)

    def test_all_writes_no_list(self):
        self.assertIsNone(ufg_subset.listed(["U-235"], [], None))

    def test_kept_nuclides_must_be_added_and_added_ones_transported(self):
        with self.assertRaises(SystemExit):
            ufg_subset.listed(["U-235", "Xe-131"], ["Xe-131"], ["U-235"])
        with self.assertRaises(SystemExit):
            ufg_subset.listed(["U-235"], ["Xe-131"], [])


class Logs(unittest.TestCase):
    def test_shielding_spans_are_read_in_seconds(self):
        log = (
            "2026 INFO shielding:mix: xlc_lattice::shielding: close time.busy=8.88ms time.idle=2.29µs\n"
            "2026 INFO shielding: xlc_lattice::shielding: close time.busy=22.2s time.idle=2.26µs\n"
            "2026 INFO shielding: xlc_lattice::shielding: close time.busy=950ms time.idle=2.26µs\n"
        )
        seconds = ufg_subset.shielding_seconds(log)
        self.assertEqual(len(seconds), 2)
        for a, b in zip(seconds, [22.2, 0.95]):
            self.assertAlmostEqual(a, b, places=12)


def one_state_table(burnup, k):
    """燃焼点ごとに 1 状態のブランチ計算のテーブル"""
    return {
        "axes": [
            {"name": "burnup", "points": burnup},
            {"name": "boron", "points": [1300.0]},
        ],
        "points": [[i, 0] for i in range(len(burnup))],
        "states": [{"k_inf": x} for x in k],
    }


class Table(unittest.TestCase):
    def table(self, burnup, k):
        return one_state_table(burnup, k)

    def test_k_is_read_by_burnup(self):
        self.assertEqual(
            ufg_subset.k_by_burnup(self.table([0.0, 10.0], [1.2, 1.1])),
            {0.0: 1.2, 10.0: 1.1},
        )

    def test_tables_with_branches_are_refused(self):
        t = self.table([0.0], [1.2])
        t["axes"][1]["points"] = [600.0, 1300.0]
        with self.assertRaises(SystemExit):
            ufg_subset.k_by_burnup(t)


def write_case(fuel_dir, case, k, timed, listed=None):
    """組の入力（`listed` は `ufg-nuclides`）、テーブル（out/）と、自己遮蔽の区間の時間のログ（out/ と、排他のジョブの out-timed/）"""
    line = "" if listed is None else f"ufg-nuclides = {listed}\n".replace("'", '"')
    (fuel_dir / f"{case}.toml").write_text(f'library = "x"\n{line}')
    table = one_state_table([0.0, 10.0], k)
    (fuel_dir / "out" / f"{case}.msgpack").write_bytes(msgpack.packb(table))
    span = "INFO shielding: xlc_lattice::shielding: close time.busy={} time.idle=1µs\n"
    (fuel_dir / "out" / f"{case}.log").write_text(
        "".join(span.format("99s") for _ in k)
    )
    (fuel_dir / "out-timed" / f"{case}.log").write_text(
        "".join(span.format(t) for t in timed)
    )


class Compare(unittest.TestCase):
    """k は out/ のテーブルから、自己遮蔽の時間は排他のジョブのログ（out-timed/）から読む"""

    def test_seconds_come_from_the_timed_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            fuel_dir = root / "ufg-subset" / "u"
            (fuel_dir / "out").mkdir(parents=True)
            (fuel_dir / "out-timed").mkdir()
            (root / "trial").mkdir()
            (root / "trial" / "added-65.txt").write_text("Xe-131\nPm-147\n")
            for case in ufg_subset.SETS:
                write_case(fuel_dir, case, [1.2, 1.1], ["4s", "8s"])
            write_case(fuel_dir, "s4", [1.2, 1.1], ["2s", "6s"], ["U-235", "Xe-131"])
            d = ufg_subset.compare("u", root, "trial")
            s4 = next(c for c in d["cases"] if c["case"] == "s4")
            self.assertEqual(s4["seconds"], [2.0, 6.0])
            self.assertEqual(s4["ratio"], [0.5, 0.75])
            # UFG で解く足した核種の数は、解いた入力の `ufg-nuclides` から数える
            self.assertEqual(s4["ufg_added"], 1)
            self.assertIsNone(d["cases"][0]["ufg_added"])
            paths = [i["path"] for i in d["inputs"]]
            self.assertIn("issue-0292/ufg-subset/u/out/s4.msgpack", paths)
            self.assertIn("issue-0292/ufg-subset/u/out-timed/s4.log", paths)
            self.assertIn("issue-0292/ufg-subset/u/s4.toml", paths)


class Report(unittest.TestCase):
    def summary(self):
        def case(name, pcm, seconds):
            return {
                "case": name,
                "pcm": pcm,
                "seconds": seconds,
                "ratio": [s / a for s, a in zip(seconds, [4.0, 8.0])],
            }

        return {
            "fuel": "u",
            "burnup": [0.0, 10.0],
            "cases": [case("all", [0.0, 0.0], [4.0, 8.0]),
                      case("s4", [0.04, -12.345], [2.0, 6.0])],
        }  # fmt: skip

    def test_the_k_table_has_a_row_per_fuel_and_subset_but_all(self):
        text = ufg_subset.k_table([self.summary()])
        self.assertIn("| 0 ", text.splitlines()[0] + " ")
        rows = [line for line in text.splitlines()[2:]]
        self.assertEqual(len(rows), 1)
        self.assertIn("U s4", rows[0])
        self.assertIn("| 0.0 ", rows[0])
        self.assertIn("−12.3", rows[0])

    def test_fuels_are_named_as_in_the_record(self):
        summary = self.summary()
        summary["fuel"] = "gd"
        self.assertIn("Gd s4", ufg_subset.k_table([summary]))

    def test_an_untimed_set_has_a_dash_in_the_time_table(self):
        summary = self.summary()
        summary["cases"][1]["seconds"] = None
        summary["cases"][1]["ratio"] = None
        rows = ufg_subset.time_table([summary]).splitlines()[2:]
        self.assertIn("U s4", rows[1])
        self.assertIn("—", rows[1])

    def test_the_time_table_gives_seconds_and_the_ratio_to_all(self):
        text = ufg_subset.time_table([self.summary()])
        rows = text.splitlines()[2:]
        self.assertEqual(len(rows), 2)
        self.assertIn("U all", rows[0])
        self.assertIn("2.00（0.50）", rows[1])
        self.assertIn("6.00（0.75）", rows[1])


class ThresholdSets(unittest.TestCase):
    def dilution(self, values):
        return {"nuclides": [{"nuclide": n, "max_pcm": v} for n, v in values]}

    def test_a_set_takes_the_nuclides_above_its_threshold_in_any_fuel(self):
        sets = ufg_subset.threshold_sets(
            [
                self.dilution([("A", 5.0), ("B", 0.5), ("C", 0.05)]),
                self.dilution([("B", -2.0), ("C", 0.2), ("D", 0.01)]),
            ],
            (1.0, 0.1),
        )
        self.assertEqual(list(sets), ["t1", "t0.1"])
        self.assertEqual(sets["t1"], ["A", "B"])
        # 燃料を通じた最大の絶対値の大きい順
        self.assertEqual(sets["t0.1"], ["A", "B", "C"])

    def test_the_sets_of_another_ranking_take_its_prefix(self):
        sets = ufg_subset.threshold_sets([self.dilution([("A", 5.0)])], (1.0,), "a")
        self.assertEqual(sets, {"a1": ["A"]})

    def test_the_listed_rankings_read_only_their_files(self):
        # 全核種の組成の順位（a）は、t の組の all を解いた後でなければ出せない。t の組の入力はその前に書ける
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "all").mkdir()
            for fuel in ufg_subset.FUELS:
                (root / "all" / f"dilution-{fuel}.json").write_text(
                    json.dumps(self.dilution([("A", 5.0)]))
                )
            sets = ufg_subset.sets("all", root, ("t",))
            self.assertEqual(sets["t1"], ["A"])
            self.assertNotIn("a1", sets)

    def test_the_trial_series_keeps_its_fixed_sets(self):
        self.assertIs(ufg_subset.sets("trial", Path("unused")), ufg_subset.SETS)


# 全核種の系列の比べ方の組（名前 → UFG で解く足した核種）
ALL_CASES = {"all": None, "none": [], "t1": ["A"]}


class AllSeries(unittest.TestCase):
    def test_the_history_is_found_from_the_input_directory(self):
        for series, path in [
            ("trial", "r/trial/u/nominal.msgpack"),
            ("all", "r/all/u/nominal.msgpack"),
        ]:
            found = ufg_subset.history_file(Path("r"), ufg_subset.SERIES[series], "u")
            self.assertEqual(found.resolve(), Path(path).resolve())

    def test_the_inputs_read_the_all_nuclides_library_and_history(self):
        text = ufg_subset.case_input(
            "u", "t1", ["U-235", "H-1"], 1, ufg_subset.SERIES["all"]
        )
        self.assertIn('library = "../../all/library/depletion.msgpack"', text)
        self.assertIn('file = "../../all/u/nominal.msgpack"', text)
        self.assertIn('ufg-nuclides = ["U-235", "H-1"]', text)
        self.assertIn("全核種", text.splitlines()[0])

    def test_the_control_reads_the_control_library_and_the_trial_history(self):
        text = ufg_subset.control_input("mox")
        self.assertIn('library = "../../all/control/depletion.msgpack"', text)
        self.assertIn('file = "../../trial/mox/nominal.msgpack"', text)
        self.assertNotIn("ufg-nuclides", text)

    def write_series(self, root, transported):
        fuel_dir = root / "all-ufg" / "u"
        (fuel_dir / "out").mkdir(parents=True)
        (fuel_dir / "out-timed").mkdir()
        (root / "ufg-subset" / "u" / "out").mkdir(parents=True, exist_ok=True)
        (root / "all" / "u").mkdir(parents=True, exist_ok=True)
        (root / "all" / "library").mkdir()
        for case, kept in ALL_CASES.items():
            listed = None if kept is None else ["U-235", *kept]
            write_case(fuel_dir, case, [1.2, 1.1], ["4s", "8s"], listed)
        (root / "all" / "added-to-p2.txt").parent.mkdir(parents=True, exist_ok=True)
        (fuel_dir / "control.toml").write_text('library = "c"\n')
        (root / "ufg-subset" / "u" / "all.toml").parent.mkdir(
            parents=True, exist_ok=True
        )
        (root / "ufg-subset" / "u" / "all.toml").write_text('library = "t"\n')
        (fuel_dir / "out" / "control.msgpack").write_bytes(
            msgpack.packb(one_state_table([0.0, 10.0], [1.2, 1.0]))
        )
        (root / "ufg-subset" / "u" / "out" / "all.msgpack").write_bytes(
            msgpack.packb(one_state_table([0.0, 10.0], [1.2, 1.1]))
        )
        (root / "all" / "library" / "xs.toml").write_text(
            "[nuclide]\n"
            'U-235 = { endf = "u", temperatures = [565] }\n'
            'A = { endf = "a", temperatures = [565] }\n'
            'H-1 = { endf = "h", temperatures = [550, 600] }\n'
            'Th-232 = { endf = "t", temperatures = [565], activation = true }\n'
        )
        # 燃焼の入力: 燃焼させる物質 fuel の温度（被覆管の温度は燃焼の領域でないので数えない）
        (root / "all" / "u" / "nominal.toml").write_text(
            "[material.fuel]\ntemperature = 565\n"
            "[material.zircaloy]\ntemperature = 600\n"
            '[depletion]\nmaterials = ["fuel"]\n'
        )
        (root / "all" / "u" / "nominal.msgpack").write_bytes(
            msgpack.packb({"transportable": transported})
        )
        (root / "all" / "added-to-p2.txt").write_text("A\n")

    def test_the_comparison_gives_the_library_change_from_the_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["A", "U-235"])
            d = ufg_subset.compare("u", root, "all", ALL_CASES)
            self.assertEqual([c["case"] for c in d["cases"]], list(ALL_CASES))
            self.assertEqual([c["ufg_added"] for c in d["cases"]], [None, 0, 1])
            paths = [i["path"] for i in d["inputs"]]
            self.assertIn("issue-0292/all-ufg/u/control.toml", paths)
            self.assertIn("issue-0292/ufg-subset/u/all.toml", paths)
            change = d["library_change_pcm"]
            self.assertAlmostEqual(change[0], 0.0)
            self.assertAlmostEqual(change[1], 1e5 * (1 / 1.1 - 1 / 1.0))

    def test_an_untimed_set_reads_no_timed_log(self):
        # 排他のジョブでまだ解いていない組（区分 decision の k だけ）は、時間を書かず、ログを入力に挙げない
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["A", "U-235"])
            (root / "all-ufg" / "u" / "out-timed" / "t1.log").unlink()
            d = ufg_subset.compare("u", root, "all", ALL_CASES, untimed={"t1"})
            t1 = next(c for c in d["cases"] if c["case"] == "t1")
            self.assertIsNone(t1["seconds"])
            self.assertIsNone(t1["ratio"])
            paths = [i["path"] for i in d["inputs"]]
            self.assertNotIn("issue-0292/all-ufg/u/out-timed/t1.log", paths)
            self.assertIn("issue-0292/all-ufg/u/out/t1.msgpack", paths)

    def test_an_untimed_name_must_be_a_set(self):
        # 書き違えた名前を黙って無視すると、残っている古いログの時間が測った時間として入る
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["A", "U-235"])
            with self.assertRaises(SystemExit):
                ufg_subset.compare("u", root, "all", ALL_CASES, untimed={"t01"})

    def test_the_reference_set_must_be_timed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["A", "U-235"])
            with self.assertRaises(SystemExit):
                ufg_subset.compare("u", root, "all", ALL_CASES, untimed={"all"})

    def test_a_history_transporting_other_nuclides_than_the_library_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["U-235"])
            with self.assertRaises(SystemExit):
                ufg_subset.compare("u", root, "all", ALL_CASES)

    def test_the_transported_nuclides_are_read_from_the_library_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, [])
            self.assertEqual(
                ufg_subset.library_transported(root / "all" / "library" / "xs.toml"),
                ["U-235", "A", "H-1"],
            )

    def test_the_history_leaves_out_nuclides_without_the_burnt_temperature(self):
        # H-1 は燃焼させる物質の温度 565 K の点を持たない（S(α, β) の内挿の 4 点）。燃焼は H-1 を崩壊だけで追い、履歴の
        # `transportable` に入れない
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_series(root, ["A", "U-235"])
            self.assertEqual(
                ufg_subset.burnt_transported(root, ufg_subset.SERIES["all"], "u"),
                ["U-235", "A"],
            )


class VersionSeries(unittest.TestCase):
    """版ごとの全核種の系列（all_nuclides.py の `VERSIONS`）: 基準の組 `all` と、その版の順位で選んだ ADR-0270 の組を、版のライブラリ、
    履歴、チェーンで解く"""

    def test_each_other_version_has_a_series_of_its_places(self):
        for key in ("endfb70", "endfb71", "endfb81"):
            with self.subTest(key):
                config = ufg_subset.SERIES[f"all-{key}"]
                self.assertEqual(config["directory"], f"all-ufg-{key}")
                self.assertEqual(config["library_input"], f"all-{key}/library/xs.toml")
                self.assertEqual(config["added"], f"all-{key}/added-to-p2.txt")

    def write_rankings(self, root, key, values):
        """版の置き場の全核種の組成の自己遮蔽の目安（燃料ごとに同じ値）"""
        place = root / f"all-{key}"
        place.mkdir(parents=True, exist_ok=True)
        for fuel in ufg_subset.FUELS:
            (place / f"dilution-all-{fuel}.json").write_text(
                json.dumps(
                    {"nuclides": [{"nuclide": n, "max_pcm": v} for n, v in values]}
                )
            )

    def test_a_version_series_solves_all_and_the_set_chosen_by_its_own_ranking(self):
        # ADR-0270 の閾値の組を、VIII.0 の順位でなく、その版の全核種の組成の目安で選ぶ
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_rankings(root, "endfb70", [("A", 0.5), ("B", 0.05)])
            # VIII.0 の順位（all/）は読まない
            (root / "all").mkdir()
            self.assertEqual(
                ufg_subset.sets("all-endfb70", root), {"all": None, "a0.1": ["A"]}
            )

    def test_a_version_series_names_its_own_rankings_as_sources(self):
        root = Path("r")
        self.assertEqual(
            ufg_subset.set_sources("all-endfb81", root),
            [root / "all-endfb81" / f"dilution-all-{f}.json" for f in ufg_subset.FUELS],
        )

    def test_a_version_series_is_compared_without_timed_logs(self):
        # 版の系列は区分 decision のジョブで解き、時間を測らない（out-timed/ が無い）
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_rankings(root, "endfb70", [("A", 0.5)])
            place = root / "all-endfb70"
            (place / "library").mkdir()
            (place / "u").mkdir()
            (place / "added-to-p2.txt").write_text("A\n")
            (place / "library" / "xs.toml").write_text(
                "[nuclide]\n"
                'U-235 = { endf = "u", temperatures = [565] }\n'
                'A = { endf = "a", temperatures = [565] }\n'
            )
            (place / "u" / "nominal.toml").write_text(
                '[material.fuel]\ntemperature = 565\n[depletion]\nmaterials = ["fuel"]\n'
            )
            (place / "u" / "nominal.msgpack").write_bytes(
                msgpack.packb({"transportable": ["U-235", "A"]})
            )
            fuel_dir = root / "all-ufg-endfb70" / "u"
            (fuel_dir / "out").mkdir(parents=True)
            (fuel_dir / "out-timed").mkdir()
            write_case(fuel_dir, "all", [1.2, 1.1], ["4s", "8s"])
            write_case(fuel_dir, "a0.1", [1.2, 1.0], ["2s", "4s"], ["U-235", "A"])
            shutil.rmtree(fuel_dir / "out-timed")
            d = ufg_subset.compare("u", root, "all-endfb70")
            self.assertEqual([c["case"] for c in d["cases"]], ["all", "a0.1"])
            for c in d["cases"]:
                self.assertIsNone(c["seconds"])
                self.assertIsNone(c["ratio"])
            self.assertAlmostEqual(d["cases"][1]["pcm"][1], 1e5 * (1 / 1.1 - 1 / 1.0))
            self.assertEqual(d["cases"][1]["ufg_added"], 1)
            self.assertNotIn("library_change_pcm", d)
            paths = [i["path"] for i in d["inputs"]]
            self.assertIn("issue-0292/all-endfb70/dilution-all-u.json", paths)
            self.assertIn("issue-0292/all-ufg-endfb70/u/out/a0.1.msgpack", paths)

    def test_a_version_input_reads_the_library_history_and_chain_of_its_version(self):
        text = ufg_subset.case_input(
            "gd", "all", None, 0, ufg_subset.SERIES["all-endfb70"]
        )
        self.assertIn('library = "../../all-endfb70/library/depletion.msgpack"', text)
        self.assertIn('file = "../../all-endfb70/gd/nominal.msgpack"', text)
        self.assertIn('chain = "../../../endfb70/chain.msgpack"', text)
        self.assertNotIn('chain = "../../../chain.msgpack"', text)
        self.assertNotIn("ufg-nuclides", text)
        self.assertIn("ENDF/B-VII.0", text.splitlines()[0])

    def test_the_viii0_series_keeps_the_production_chain(self):
        text = ufg_subset.case_input("u", "all", None, 0, ufg_subset.SERIES["all"])
        self.assertIn('chain = "../../../chain.msgpack"', text)


class CheckFuel(unittest.TestCase):
    """確認のセル（gd10）: 3 つの燃料の順位で選んだ組を、選んだ組成と別の履歴で解く。順位に入れず、時間を測らない"""

    def write_check(self, root):
        fuel_dir = root / "all-ufg" / "gd10"
        (fuel_dir / "out").mkdir(parents=True)
        (root / "all" / "gd10").mkdir(parents=True)
        (root / "all" / "library").mkdir()
        (root / "all" / "added-to-p2.txt").write_text("A\n")
        # 組は 3 つの燃料の順位から（a0.1 = [A]）
        for name in ufg_subset.RANKINGS.values():
            for fuel in ufg_subset.FUELS:
                (root / "all" / name.format(fuel=fuel)).write_text(
                    json.dumps({"nuclides": [{"nuclide": "A", "max_pcm": 5.0}]})
                )
        (root / "all" / "library" / "xs.toml").write_text(
            "[nuclide]\n"
            'U-235 = { endf = "u", temperatures = [565] }\n'
            'A = { endf = "a", temperatures = [565] }\n'
        )
        (root / "all" / "gd10" / "nominal.toml").write_text(
            '[material.fuel]\ntemperature = 565\n[depletion]\nmaterials = ["fuel"]\n'
        )
        (root / "all" / "gd10" / "nominal.msgpack").write_bytes(
            msgpack.packb({"transportable": ["U-235", "A"]})
        )
        for case, k, listed in [
            ("all", [1.2, 1.1], None),
            ("a0.1", [1.2, 1.0], ["U-235", "A"]),
        ]:
            line = (
                "" if listed is None else f"ufg-nuclides = {listed}\n".replace("'", '"')
            )
            (fuel_dir / f"{case}.toml").write_text(f'library = "x"\n{line}')
            (fuel_dir / "out" / f"{case}.msgpack").write_bytes(
                msgpack.packb(one_state_table([0.0, 10.0], k))
            )
        return fuel_dir

    def test_the_check_fuel_solves_all_and_the_chosen_set(self):
        self.assertEqual(ufg_subset.CHECK_SETS, ("all", "a0.1"))
        self.assertNotIn("gd10", ufg_subset.FUELS)

    def test_the_check_fuel_reads_no_timed_log_and_no_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            self.write_check(root)
            d = ufg_subset.compare("gd10", root, "all")
            self.assertEqual([c["case"] for c in d["cases"]], ["all", "a0.1"])
            for c in d["cases"]:
                self.assertIsNone(c["seconds"])
                self.assertIsNone(c["ratio"])
            self.assertAlmostEqual(d["cases"][1]["pcm"][1], 1e5 * (1 / 1.1 - 1 / 1.0))
            self.assertEqual(d["cases"][1]["ufg_added"], 1)
            self.assertNotIn("library_change_pcm", d)
            paths = [i["path"] for i in d["inputs"]]
            self.assertFalse(any("out-timed" in p for p in paths))
            self.assertIn("issue-0292/all/dilution-all-gd.json", paths)

    def test_the_check_fuel_inputs_read_its_history(self):
        text = ufg_subset.case_input(
            "gd10", "a0.1", ["U-235", "A"], 1, ufg_subset.SERIES["all"]
        )
        self.assertIn('file = "../../all/gd10/nominal.msgpack"', text)
        self.assertIn('rows = ["pin pin pin", "pin gd  pin", "pin pin pin"]', text)

    def test_the_check_fuel_is_named_in_the_report(self):
        summary = {
            "fuel": "gd10",
            "burnup": [0.0],
            "cases": [{"case": "a0.1", "pcm": [0.5], "seconds": None, "ratio": None}],
        }
        self.assertIn("Gd 10 wt% a0.1", ufg_subset.k_table([summary]))


class Sources(unittest.TestCase):
    def test_the_dilution_sets_name_their_json_sources(self):
        root = Path("r")
        self.assertEqual(ufg_subset.set_sources("trial", root), [])
        self.assertEqual(
            ufg_subset.set_sources("all", root),
            [root / "all" / f"dilution-{f}.json" for f in ufg_subset.FUELS]
            + [root / "all" / f"dilution-all-{f}.json" for f in ufg_subset.FUELS],
        )
        self.assertEqual(
            ufg_subset.set_sources("all", root, ("a",)),
            [root / "all" / f"dilution-all-{f}.json" for f in ufg_subset.FUELS],
        )
        with self.assertRaises(SystemExit):
            ufg_subset.set_sources("all", root, ("x",))


class LibraryTable(unittest.TestCase):
    def test_the_report_has_the_library_change_per_fuel(self):
        summaries = [
            {"fuel": "u", "burnup": [0.0, 10.0], "library_change_pcm": [0.0, -0.17]},
            {"fuel": "gd", "burnup": [0.0, 10.0], "library_change_pcm": [0.01, 0.3]},
        ]
        rows = ufg_subset.library_table(summaries).splitlines()[2:]
        self.assertEqual(len(rows), 2)
        self.assertIn("U", rows[0])
        self.assertIn("−0.2", rows[0])
        self.assertIn("Gd", rows[1])


if __name__ == "__main__":
    unittest.main()
