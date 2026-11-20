"""ufg_subset.py（#0292、ADR-0268 の UFG で解く核種の組）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_ufg_subset.py'
"""

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


class Sources(unittest.TestCase):
    def test_the_dilution_sets_name_their_json_sources(self):
        root = Path("r")
        self.assertEqual(ufg_subset.set_sources("trial", root), [])
        self.assertEqual(
            ufg_subset.set_sources("all", root),
            [root / "all" / f"dilution-{f}.json" for f in ufg_subset.FUELS],
        )


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
