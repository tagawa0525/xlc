"""compare_forgetting.py（追わない核種の分 U を、指数の忘却の履歴の変数の重みで R(1) から内挿したときの外れ）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_forgetting.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_forgetting as c
import compare_index as ci
import compare_split as s
from test_compare_split import tables

# switch が条件を変えた燃焼度（fixture の参照の点は 0 と 20 MWd/kg）
START = 10.0


class Weight(unittest.TestCase):
    def test_an_infinite_time_constant_is_the_burnup_average(self):
        # τ → ∞ の H は燃焼度で重みを付けた平均（ADR-0215 の案 V1）で、w = (B − B₀) / B
        self.assertAlmostEqual(c.weight(60.0, 20.0, None), 40.0 / 60.0, delta=1e-15)
        self.assertAlmostEqual(c.weight(60.0, 20.0, 1e9), 40.0 / 60.0, delta=1e-7)

    def test_the_weight_is_the_exponential_average_of_the_switched_part(self):
        # H_τ(B) = ∫ x e^{−(B − B')/τ} dB' / ∫ e^{−(B − B')/τ} dB'（0 から B）。B₀ より後が条件 j
        tau, b, b0 = 15.0, 60.0, 20.0
        expected = (1.0 - math.exp(-(b - b0) / tau)) / (1.0 - math.exp(-b / tau))
        self.assertAlmostEqual(c.weight(b, b0, tau), expected, delta=1e-15)
        # 最近の燃焼に大きな重み: 線形の平均より大きい
        self.assertGreater(c.weight(b, b0, tau), (b - b0) / b)
        # 負の τ は古い燃焼に大きな重み
        self.assertLess(c.weight(b, b0, -tau), (b - b0) / b)

    def test_histories_switched_at_burnup_zero_have_the_weight_one(self):
        for tau in (None, 5.0, -5.0):
            self.assertAlmostEqual(c.weight(30.0, 0.0, tau), 1.0, delta=1e-15)


class Compare(unittest.TestCase):
    def test_gaps_are_the_reference_less_the_part_interpolated_at_the_weight(self):
        t = tables()
        result = c.compare(t, START, [None, 30.0, -30.0])
        _, one = s.split(t["nominal"], t["one"], t["one_hybrid"], "R(1)")
        _, switch = s.split(t["nominal"], t["switch"], t["switch_hybrid"], "switch")
        x, o = switch[-1], one[-1]
        for row in result["rows"]:
            w = c.weight(20.0, START, row["tau"])
            self.assertEqual([p["burnup"] for p in row["points"]], [20.0])
            self.assertAlmostEqual(row["points"][0]["weight"], w, delta=1e-15)
            self.assertAlmostEqual(
                row["points"][0]["gap_pcm"],
                ci.gap(x["branch"], x["untracked"], o["untracked"], w),
                delta=1e-12,
            )
            self.assertEqual(
                row["largest_pcm"], max(abs(p["gap_pcm"]) for p in row["points"])
            )

    def test_the_best_time_constant_has_the_smallest_largest_gap(self):
        result = c.compare(tables(), START, [None, 30.0, -30.0, 5.0])
        best = min(result["rows"], key=lambda r: r["largest_pcm"])
        self.assertEqual(result["best"], best["tau"])

    def test_points_up_to_the_switch_are_left_out(self):
        t = tables()
        result = c.compare(t, 20.0, [None])
        self.assertEqual(result["rows"][0]["points"], [])


class Report(unittest.TestCase):
    def test_report_names_the_best_time_constant_and_its_gaps(self):
        text = c.report(c.compare(tables(), START, [None, 30.0, -30.0]))
        self.assertIn("τ", text)
        self.assertIn("最小", text)
        self.assertNotRegex(text, r"(?<![A-Za-z])-[0-9]")


if __name__ == "__main__":
    unittest.main()
