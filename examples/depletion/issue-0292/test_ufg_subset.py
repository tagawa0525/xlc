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


class Compare(unittest.TestCase):
    """k は out/ のテーブルから、自己遮蔽の時間は排他のジョブのログ（out-timed/）から読む"""

    def write_case(self, fuel_dir, case, k, timed):
        table = one_state_table([0.0, 10.0], k)
        (fuel_dir / "out" / f"{case}.msgpack").write_bytes(msgpack.packb(table))
        span = (
            "INFO shielding: xlc_lattice::shielding: close time.busy={} time.idle=1µs\n"
        )
        (fuel_dir / "out" / f"{case}.log").write_text(
            "".join(span.format("99s") for _ in k)
        )
        (fuel_dir / "out-timed" / f"{case}.log").write_text(
            "".join(span.format(t) for t in timed)
        )

    def test_seconds_come_from_the_timed_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "issue-0292"
            fuel_dir = root / "ufg-subset" / "u"
            (fuel_dir / "out").mkdir(parents=True)
            (fuel_dir / "out-timed").mkdir()
            for case in ufg_subset.SETS:
                self.write_case(fuel_dir, case, [1.2, 1.1], ["4s", "8s"])
            self.write_case(fuel_dir, "s4", [1.2, 1.1], ["2s", "6s"])
            d = ufg_subset.compare("u", root)
            s4 = next(c for c in d["cases"] if c["case"] == "s4")
            self.assertEqual(s4["seconds"], [2.0, 6.0])
            self.assertEqual(s4["ratio"], [0.5, 0.75])
            paths = [i["path"] for i in d["inputs"]]
            self.assertIn("issue-0292/ufg-subset/u/out/s4.msgpack", paths)
            self.assertIn("issue-0292/ufg-subset/u/out-timed/s4.log", paths)


if __name__ == "__main__":
    unittest.main()
