"""p9_timing.py（core の燃焼のログから時間の内訳を取る）。

実行: python3 -m unittest discover -s examples/depletion/p9-timing -p 'test_p9_timing.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p9_timing

# 燃焼点 0（探索なし）、ステップの中点（平衡 Xe）、ステップの終わりの点（平衡 Xe の後に別の状態で探索）、崩壊だけのステップの後の点
# （探索）。固有値計算の秒は span のフィールド（3 桁に丸めた time.busy でなく）
LOG = """\
2026-11-05T15:22:59.236621Z  INFO xlc_core::depletion: micro depletion tracked=31 auxiliary=143 nodes=9843
2026-11-05T15:22:59.386358Z  INFO constants: xlc_core::case: close time.busy=150ms time.idle=3.45µs
2026-11-05T15:24:57.000000Z  INFO nem{iterations=4165 seconds=118.25}: xlc_core::nem: close time.busy=118s time.idle=2.40µs
2026-11-05T15:24:57.100000Z  INFO xlc_core::depletion: depletion point efpd=0.0 days=0.0 state="case-2" k=1.17 boron=1200.0 solves=1
2026-11-05T15:24:58.000000Z  INFO depletion step{efpd=9.0}:micro depletion: xlc_core::depletion: close time.busy=900ms time.idle=1.2µs
2026-11-05T15:24:58.200000Z  INFO depletion step{efpd=9.0}:constants: xlc_core::case: close time.busy=200ms time.idle=3µs
2026-11-05T15:25:58.200000Z  INFO depletion step{efpd=9.0}:nem{iterations=2000 seconds=60.5}: xlc_core::nem: close time.busy=60.5s time.idle=3µs
2026-11-05T15:25:58.300000Z  INFO depletion step{efpd=9.0}:xenon balance: xlc_core::depletion: close time.busy=500ms time.idle=3µs
2026-11-05T15:25:58.350000Z  INFO depletion step{efpd=9.0}:xenon anderson: xlc_core::depletion: close time.busy=50ms time.idle=3µs
2026-11-05T15:25:58.400000Z  INFO depletion step{efpd=9.0}: xlc_core::depletion: equilibrium xenon iterations=2 residual=1e-7
2026-11-05T15:25:58.500000Z  INFO depletion step{efpd=9.0}: xlc_core::depletion: depletion midpoint k=1.13 boron=1200.0 solves=2
2026-11-05T15:25:58.600000Z  INFO depletion step{efpd=9.0}: xlc_core::depletion: close time.busy=61.9s time.idle=3µs
2026-11-05T15:26:28.600000Z  INFO nem{iterations=1000 seconds=30.0}: xlc_core::nem: close time.busy=30.0s time.idle=3µs
2026-11-05T15:26:28.700000Z  INFO xlc_core::depletion: equilibrium xenon iterations=1 residual=1e-7
2026-11-05T15:26:28.800000Z  INFO search:constants: xlc_core::case: close time.busy=150ms time.idle=3µs
2026-11-05T15:27:28.800000Z  INFO search:nem{iterations=1500 seconds=40.0}: xlc_core::nem: close time.busy=40.0s time.idle=3µs
2026-11-05T15:27:28.900000Z  INFO search: xlc_core::depletion: close time.busy=40.2s time.idle=3µs
2026-11-05T15:27:29.000000Z  INFO xlc_core::depletion: depletion point efpd=9.0 days=13.7 state="case-3" k=1.136 boron=1210.5 solves=2
2026-11-05T15:27:29.100000Z  INFO decay step{efpd=9.0 days=0.01}:micro decay: xlc_core::depletion: close time.busy=300ms time.idle=3µs
2026-11-05T15:27:29.200000Z  INFO decay step{efpd=9.0 days=0.01}: xlc_core::depletion: close time.busy=310ms time.idle=3µs
2026-11-05T15:27:29.300000Z  INFO search:nem{iterations=900 seconds=20.0}: xlc_core::nem: close time.busy=20.0s time.idle=3µs
2026-11-05T15:27:29.400000Z  INFO search: xlc_core::depletion: close time.busy=20.1s time.idle=3µs
2026-11-05T15:27:29.500000Z  INFO xlc_core::depletion: depletion point efpd=9.0 days=13.71 state="case-3" k=1.136 boron=1205.0 solves=1
2026-11-05T15:27:29.600000Z  INFO core: core seconds=272.0
"""


class Summary(unittest.TestCase):
    def setUp(self):
        self.summary = p9_timing.summarize(LOG)

    def test_the_exclusive_parts_and_the_rest_add_up_to_the_total(self):
        s = self.summary["seconds"]
        self.assertAlmostEqual(s["total"], 272.0)
        self.assertAlmostEqual(s["eigenvalue"], 268.75)
        self.assertAlmostEqual(s["interpolation"], 0.5)
        self.assertAlmostEqual(s["micro"], 1.2)
        self.assertAlmostEqual(s["xenon"], 0.55)
        self.assertAlmostEqual(
            s["eigenvalue"] + s["interpolation"] + s["micro"] + s["xenon"] + s["other"],
            s["total"],
        )

    def test_the_search_time_holds_its_solves_and_interpolations(self):
        self.assertAlmostEqual(self.summary["seconds"]["search"], 60.3)

    def test_eigenvalue_solves_give_their_count_iterations_and_times(self):
        nem = self.summary["nem"]
        self.assertEqual(nem["count"], 5)
        self.assertEqual(nem["iterations"], [900, 4165])
        self.assertEqual(nem["seconds"], [20.0, 118.25])

    def test_points_and_midpoints_give_their_solves_in_order(self):
        self.assertEqual(
            self.summary["solves"],
            [
                {
                    "kind": "point",
                    "efpd": 0.0,
                    "boron": 1200.0,
                    "solves": 1,
                    "xenon": None,
                },
                {
                    "kind": "midpoint",
                    "efpd": 9.0,
                    "boron": 1200.0,
                    "solves": 2,
                    "xenon": 2,
                },
                # ステップの終わりの点: 前の状態の平衡 Xe の解き直しの後に探索した
                {
                    "kind": "point",
                    "efpd": 9.0,
                    "boron": 1210.5,
                    "solves": 2,
                    "xenon": 1,
                },
                {
                    "kind": "point",
                    "efpd": 9.0,
                    "boron": 1205.0,
                    "solves": 1,
                    "xenon": None,
                },
            ],
        )

    def test_midpoints_read_the_step_efpd_among_other_fields(self):
        line = LOG.splitlines()[10].replace("step{efpd=9.0}", "step{index=1 efpd=9.0}")
        log = LOG.replace(LOG.splitlines()[10], line)
        self.assertEqual(p9_timing.summarize(log)["solves"][1]["efpd"], 9.0)

    def test_the_summary_names_its_format(self):
        self.assertEqual(self.summary["format"], "xlc-core-timing")

    def test_durations_read_every_unit(self):
        self.assertEqual(p9_timing.seconds("3.45µs"), 3.45e-6)
        self.assertEqual(p9_timing.seconds("150ms"), 0.15)
        self.assertEqual(p9_timing.seconds("12ns"), 1.2e-8)
        self.assertEqual(p9_timing.seconds("1.5s"), 1.5)

    def test_a_log_without_the_total_is_an_error(self):
        with self.assertRaises(ValueError):
            p9_timing.summarize(LOG.replace("core seconds=272.0", "core"))


if __name__ == "__main__":
    unittest.main()
