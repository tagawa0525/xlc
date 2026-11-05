"""p9_timing.py（core の燃焼のログから時間の内訳を取る）。

実行: python3 -m unittest discover -s examples/depletion/p9-timing -p 'test_p9_timing.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p9_timing

# 燃焼点 0（探索なし）、ステップの中点の探索（固有値計算 2 回と平衡 Xe）と、その間のミクロ燃焼
LOG = """\
2026-11-05T15:22:59.236621Z  INFO xlc_core::depletion: micro depletion tracked=31 auxiliary=143 nodes=9843
2026-11-05T15:22:59.386358Z  INFO constants: xlc_core::case: close time.busy=150ms time.idle=3.45µs
2026-11-05T15:24:57.000000Z  INFO nem{iterations=4165}: xlc_core::nem: close time.busy=118s time.idle=2.40µs
2026-11-05T15:24:57.100000Z  INFO xlc_core::depletion: depletion point efpd=0.0 days=0.0 state="case-2" k=1.17 boron=1200.0 solves=1
2026-11-05T15:24:58.000000Z  INFO depletion step{efpd=32.0}:micro depletion: xlc_core::depletion: close time.busy=900ms time.idle=1.2µs
2026-11-05T15:24:58.200000Z  INFO depletion step{efpd=32.0}:search:constants: xlc_core::case: close time.busy=200ms time.idle=3µs
2026-11-05T15:25:58.200000Z  INFO depletion step{efpd=32.0}:search:nem{iterations=2000}: xlc_core::nem: close time.busy=60.0s time.idle=3µs
2026-11-05T15:25:58.400000Z  INFO depletion step{efpd=32.0}:search:constants: xlc_core::case: close time.busy=200ms time.idle=3µs
2026-11-05T15:26:28.400000Z  INFO depletion step{efpd=32.0}:search:nem{iterations=1000}: xlc_core::nem: close time.busy=30.0s time.idle=3µs
2026-11-05T15:26:28.500000Z  INFO depletion step{efpd=32.0}:search: xlc_core::depletion: close time.busy=90.5s time.idle=3µs
2026-11-05T15:26:28.600000Z  INFO depletion step{efpd=32.0}:xenon balance: xlc_core::depletion: close time.busy=500ms time.idle=3µs
2026-11-05T15:26:28.700000Z  INFO depletion step{efpd=32.0}: xlc_core::depletion: equilibrium xenon iterations=1 residual=0.0
2026-11-05T15:26:28.800000Z  INFO depletion step{efpd=32.0}: xlc_core::depletion: depletion midpoint k=1.0 boron=1150.0 solves=2
2026-11-05T15:26:29.000000Z  INFO depletion step{efpd=32.0}: xlc_core::depletion: close time.busy=91.9s time.idle=3µs
2026-11-05T15:26:29.100000Z  INFO core: core seconds=210.0
"""


class Summary(unittest.TestCase):
    def setUp(self):
        self.summary = p9_timing.summarize(LOG)

    def test_the_exclusive_parts_and_the_rest_add_up_to_the_total(self):
        s = self.summary["seconds"]
        self.assertAlmostEqual(s["total"], 210.0)
        self.assertAlmostEqual(s["eigenvalue"], 208.0)
        self.assertAlmostEqual(s["interpolation"], 0.55)
        self.assertAlmostEqual(s["micro"], 0.9)
        self.assertAlmostEqual(s["xenon"], 0.5)
        self.assertAlmostEqual(
            s["eigenvalue"] + s["interpolation"] + s["micro"] + s["xenon"] + s["other"],
            s["total"],
        )

    def test_the_search_time_holds_its_solves_and_interpolations(self):
        self.assertAlmostEqual(self.summary["seconds"]["search"], 90.5)

    def test_eigenvalue_solves_give_their_count_iterations_and_times(self):
        nem = self.summary["nem"]
        self.assertEqual(nem["count"], 3)
        self.assertEqual(nem["iterations"], [1000, 4165])
        self.assertEqual(nem["seconds"], [30.0, 118.0])

    def test_points_and_midpoints_give_their_solves(self):
        self.assertEqual(
            self.summary["points"], [{"efpd": 0.0, "boron": 1200.0, "solves": 1}]
        )
        self.assertEqual(self.summary["midpoints"], [{"boron": 1150.0, "solves": 2}])
        self.assertEqual(self.summary["xenon_iterations"], [1])

    def test_durations_read_every_unit(self):
        self.assertEqual(p9_timing.seconds("3.45µs"), 3.45e-6)
        self.assertEqual(p9_timing.seconds("150ms"), 0.15)
        self.assertEqual(p9_timing.seconds("12ns"), 1.2e-8)
        self.assertEqual(p9_timing.seconds("1.5s"), 1.5)

    def test_a_log_without_the_total_is_an_error(self):
        with self.assertRaises(ValueError):
            p9_timing.summarize(LOG.replace("core seconds=210.0", "core"))


if __name__ == "__main__":
    unittest.main()
