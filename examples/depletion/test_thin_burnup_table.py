"""thin_burnup_table.py（燃焼度の軸を持つテーブルから燃焼度の点を選んだテーブル。ADR-0192）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_thin_burnup_table.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import thin_burnup_table


def table() -> dict:
    """燃焼度 3 点 × ホウ素 2 点、部分ノード 1 つ。状態の k_inf に (燃焼度, ホウ素) を符号化する"""
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "split": 1,
        "heavy_metal": 1258.0,
        "axes": [
            {"name": "burnup", "points": [0.0, 10.0, 20.0]},
            {"name": "boron", "points": [1000.0, 1300.0], "b10": 0.199},
        ],
        "points": [[b, p] for b in range(3) for p in range(2)],
        "states": [
            {"k_inf": b + p / 1e4} for b in [0.0, 10.0, 20.0] for p in [1000.0, 1300.0]
        ],
    }


class Thin(unittest.TestCase):
    def test_keeps_the_states_of_the_chosen_burnups_in_order(self):
        out = thin_burnup_table.thin(table(), [0.0, 20.0])
        self.assertEqual(out["axes"][0]["points"], [0.0, 20.0])
        self.assertEqual([s["k_inf"] for s in out["states"]], [0.1, 0.13, 20.1, 20.13])
        self.assertEqual(out["points"], [[0, 0], [0, 1], [1, 0], [1, 1]])
        # ほかの項目はそのまま
        self.assertEqual(out["heavy_metal"], 1258.0)
        self.assertEqual(out["axes"][1], table()["axes"][1])

    def test_split_tables_keep_every_sub_node_of_a_state(self):
        t = table()
        # 2 × 2 の部分ノード: 燃焼度 3 点 × ホウ素 2 点 × 4 = 24 状態
        t["split"] = 2
        t["states"] = [{"k_inf": float(i)} for i in range(24)]
        out = thin_burnup_table.thin(t, [10.0])
        self.assertEqual(
            [s["k_inf"] for s in out["states"]], [float(i) for i in range(8, 16)]
        )

    def test_burnups_not_in_the_table_or_out_of_order_are_errors(self):
        for burnup in ([0.0, 15.0], [20.0, 0.0], [10.0, 10.0], []):
            with self.assertRaises(SystemExit, msg=str(burnup)):
                thin_burnup_table.thin(table(), burnup)

    def test_state_counts_that_do_not_fill_the_points_are_errors(self):
        t = table()
        t["states"] = t["states"][:-1]
        with self.assertRaises(SystemExit):
            thin_burnup_table.thin(t, [0.0])

    def test_points_solved_at_some_burnups_keep_their_burnups(self):
        # 公称でない状態（ホウ素 1000 ppm）を燃焼度 0 と 20 だけで解いたテーブル（ADR-0200）。燃焼度 10 を除くと、その状態は
        # 0 と 20 で残る
        t = table()
        t["points"] = [[0, 0], [0, 1], [1, 1], [2, 0], [2, 1]]
        t["states"] = [{"k_inf": float(i)} for i in range(5)]
        out = thin_burnup_table.thin(t, [0.0, 20.0])
        self.assertEqual(out["points"], [[0, 0], [0, 1], [1, 0], [1, 1]])
        self.assertEqual([s["k_inf"] for s in out["states"]], [0.0, 1.0, 3.0, 4.0])

    def test_tables_without_a_leading_burnup_axis_are_errors(self):
        t = table()
        t["axes"].reverse()
        with self.assertRaises(SystemExit):
            thin_burnup_table.thin(t, [0.0])


if __name__ == "__main__":
    unittest.main()
