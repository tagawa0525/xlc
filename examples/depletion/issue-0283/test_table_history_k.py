"""table_history_k.py の比較（テーブルの状態の k∞ と、同じ燃焼度の格子計算の燃焼の履歴の k∞ の反応度の差）。

実行: python3 -m unittest discover -s examples/depletion/issue-0283 -p 'test_table_history_k.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import table_history_k


def history():
    """格子計算の履歴（燃焼度 0、0.1、20）"""
    return {
        "points": [
            {"burnup": b, "k_inf": k} for b, k in [(0.0, 1.2), (0.1, 1.15), (20.0, 1.0)]
        ]
    }


def table(burnups, ks, boron=(1300.0,)):
    """燃焼度の点 `burnups` × ホウ素の点 `boron` のテーブル。状態の並びは燃焼度が外"""
    points = [[i, j, 0, 0] for i in range(len(burnups)) for j in range(len(boron))]
    return {
        "split": 1,
        "axes": [
            {"name": "burnup", "points": list(burnups)},
            {"name": "boron", "points": list(boron)},
            {"name": "moderator", "points": [565.0]},
            {"name": "fuel", "points": [565.0]},
        ],
        "points": points,
        "states": [{"k_inf": k} for k in ks],
    }


class Rows(unittest.TestCase):
    def test_reactivity_is_the_history_inverse_minus_the_table_inverse(self):
        rows = table_history_k.rows(history(), table([0.0, 20.0], [1.1, 1.0 + 1e-5]))
        self.assertEqual([r[0] for r in rows], [0.0, 20.0])
        self.assertEqual([r[1] for r in rows], [1.1, 1.0 + 1e-5])
        self.assertEqual([r[2] for r in rows], [1.2, 1.0])
        self.assertAlmostEqual(rows[0][3], (1 / 1.2 - 1 / 1.1) * 1e5, places=9)
        # テーブルの k∞ が大きければ正
        self.assertGreater(rows[1][3], 0.0)
        self.assertAlmostEqual(rows[1][3], (1.0 - 1 / (1.0 + 1e-5)) * 1e5, places=9)

    def test_burnups_off_the_history_points_are_errors(self):
        with self.assertRaises(SystemExit):
            table_history_k.rows(history(), table([0.0, 0.5], [1.1, 1.0]))

    def test_tables_with_several_states_per_burnup_are_errors(self):
        # 燃焼度ごとにホウ素 2 状態あると、公称でない状態も公称の履歴と比べてしまう
        with self.assertRaises(SystemExit):
            table_history_k.rows(
                history(), table([0.0], [1.1, 1.0], boron=(600.0, 1300.0))
            )

    def test_split_tables_are_errors(self):
        split = table([0.0], [1.1])
        split["split"] = 2
        with self.assertRaises(SystemExit):
            table_history_k.rows(history(), split)


if __name__ == "__main__":
    unittest.main()
