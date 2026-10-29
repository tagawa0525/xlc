"""core_history.py（core の燃焼の履歴を読む比較のツールの共通の部分）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_core_history.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core_history import single_state

STATE = {"boron": 1300.0, "b10": 0.199, "moderator": 565.0, "fuel": 565.0, "banks": {}}


class SingleState(unittest.TestCase):
    def test_one_state_histories_give_their_state(self):
        self.assertEqual(
            single_state({"states": {"nominal": dict(STATE)}, "points": []}), STATE
        )

    def test_histories_of_several_states_are_errors(self):
        # ステップごとに状態の変わる燃焼（ADR-0218）は、状態が 1 つの燃焼を比べるツールで比べない
        two = {"states": {"in": dict(STATE), "out": dict(STATE)}, "points": []}
        with self.assertRaises(SystemExit):
            single_state(two)
        with self.assertRaises(SystemExit):
            single_state({"states": {}, "points": []})

    def test_histories_off_the_rated_power_are_errors(self):
        # 出力の割合と崩壊だけのステップ（ADR-0218）のある燃焼は、定格のままの燃焼を比べるツールで比べない
        for power in [0.5, 0.0]:
            history = {
                "states": {"nominal": dict(STATE)},
                "points": [
                    {"efpd": 0.0},
                    {"efpd": 1.0, "power": 1.0},
                    {"efpd": 2.0, "power": power},
                ],
            }
            with self.assertRaises(SystemExit):
                single_state(history)
        history["points"][2]["power"] = 1.0
        self.assertEqual(single_state(history), STATE)


if __name__ == "__main__":
    unittest.main()
