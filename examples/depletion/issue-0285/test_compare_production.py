"""compare_production.py のログの読み方、状態の種類の分け方、テーブルの差（#0285、ADR-0250）。

実行: cd examples/depletion/issue-0285 && python3 -m unittest test_compare_production
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_production import breakdown, differences, kind, log

LOG = (
    "t  INFO shielding:dancoff: xlc_lattice::shielding: close time.busy=2.20s time.idle=1µs\n"
    "t  INFO shielding: xlc_lattice::shielding: close time.busy=30.0s time.idle=1µs\n"
    "t  INFO ratio-shielding: xlc_lattice::shielding: close time.busy=20.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc::symmetric: symmetry tracks=1\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=6.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3000ms time.idle=1µs\n"
    't  INFO lattice: derivative state=3 nuclide="Pu-240" iterations=14 seconds=4.0\n'
    "t  INFO lattice: lattice state=3 iterations=25 tracks=1 seconds=70.0 total_seconds=75\n"
    "t  INFO ratio-shielding: xlc_lattice::shielding: close time.busy=15.0s time.idle=1µs\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3.0s time.idle=1µs\n"
    "t  INFO lattice: lattice state=10 iterations=13 tracks=1 seconds=20.0 total_seconds=95\n"
    "t  INFO moc{regions=1 tracks=1}: xlc_lattice::moc: close time.busy=3.0s time.idle=1µs\n"
    "t  INFO lattice: lattice state=17 iterations=13 tracks=1 seconds=8.0 total_seconds=103\n"
    "t  INFO lattice: branches shieldings=1 ratios=2 derivative_shieldings=0 total_seconds=104.5\n"
)


class Log(unittest.TestCase):
    def test_spans_are_summed_per_state_and_kinds_follow_the_shielding(self):
        parsed = log(LOG)
        first, ratio, remix = parsed["states"]
        # 自己遮蔽の中の区間（dancoff）は数えず、MOC の区間は ms も秒にする
        self.assertEqual(first["shielding"], [30.0])
        self.assertEqual(first["moc"], [6.0, 3.0])
        self.assertEqual(first["derivative"], [4.0])
        self.assertEqual(
            [kind(s) for s in parsed["states"]], ["ufg+ratio", "ratio", "remix"]
        )
        self.assertEqual((ratio["state"], remix["seconds"]), (10, 8.0))

    def test_the_breakdown_counts_and_sums_by_kind(self):
        b = breakdown(log(LOG))
        self.assertEqual(b["total_seconds"], 104.5)
        self.assertEqual(b["state_seconds"], 98.0)
        self.assertEqual((b["shieldings"], b["ratios"]), (1, 2))
        self.assertEqual(b["kinds"]["ratio"]["ratio"], [1, 15.0])
        self.assertEqual(b["kinds"]["remix"]["moc"], [1, 3.0])

    def test_an_unfinished_run_stops(self):
        with self.assertRaises(SystemExit):
            log(LOG.rsplit("t  INFO lattice: branches", 1)[0])


def state(k: float, absorption: list) -> dict:
    return {
        "k_inf": k,
        "absorption": absorption,
        "nu_fission": [0.01, 0.1],
        "kappa_fission": [1.0, 2.0],
        "diffusion": [1.4, 0.4],
        "scatter": [[0.5, 0.02], [0.0, 1.0]],
        "adf": [1.0, 1.0],
    }


class Differences(unittest.TestCase):
    POINTS = ((0, 2, 2, 3), (0, 2, 6, 3))

    def test_the_worst_state_is_named_by_its_point(self):
        base = (list(self.POINTS), [state(1.1, [0.01, 0.1]), state(1.0, [0.01, 0.1])])
        variant = (
            list(self.POINTS),
            [state(1.1, [0.01, 0.1]), state(1.0002, [0.01, 0.101])],
        )
        d = differences(variant, base)
        self.assertAlmostEqual(d["worst_k"]["pcm"], 20.0)
        self.assertEqual(d["worst_k"]["point"], self.POINTS[1])
        self.assertAlmostEqual(d["items"]["absorption"]["relative"], 0.001 / 0.101)
        self.assertEqual(d["items"]["nu_fission"]["relative"], 0.0)

    def test_tables_of_other_states_stop(self):
        base = (list(self.POINTS), [state(1.1, [0.01, 0.1])] * 2)
        other = ([self.POINTS[0], (0, 2, 5, 3)], base[1])
        with self.assertRaises(SystemExit):
            differences(other, base)


if __name__ == "__main__":
    unittest.main()
