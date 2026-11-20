"""ufg_subset.py（#0292、ADR-0268 の UFG で解く核種の組）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_ufg_subset.py'
"""

import sys
import unittest
from pathlib import Path

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


class Table(unittest.TestCase):
    def table(self, burnup, k):
        return {
            "axes": [
                {"name": "burnup", "points": burnup},
                {"name": "boron", "points": [1300.0]},
            ],
            "points": [[i, 0] for i in range(len(burnup))],
            "states": [{"k_inf": x} for x in k],
        }

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


if __name__ == "__main__":
    unittest.main()
