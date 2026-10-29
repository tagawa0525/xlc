"""small_core_cases.py（小さな炉心の燃焼のケースを基準のサイクルから書く）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_small_core_cases.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import small_core_cases


class Cases(unittest.TestCase):
    def test_the_committed_cases_are_the_written_ones(self):
        for method, split in small_core_cases.CASES:
            path = small_core_cases.BRANCHES / f"core-small-{method}-{split}.toml"
            self.assertEqual(
                path.read_text(), small_core_cases.case(method, split), path.name
            )

    def test_splits_keep_the_base_points_and_the_decay_step(self):
        base = small_core_cases.steps(1)
        fine = small_core_cases.steps(4)
        # 出力のあるステップは 4 倍、崩壊だけのステップは分けない
        self.assertEqual(len(fine), 4 * (len(base) - 1) + 1)
        self.assertEqual(sum("decay-days" in line for line in fine), 1)
        # 基準のステップの終わりの EFPD とコメント、状態は最初の分けたステップだけ
        ends = [line for line in fine if "#" in line and "decay" not in line]
        self.assertEqual(
            [line.split("#")[0] for line in ends],
            [
                line.split("#")[0]
                .replace(', state = "aro"', "")
                .replace(', state = "rodded"', "")
                for line in base
                if "decay" not in line
            ],
        )
        self.assertIn('state = "aro"', fine[0])
        self.assertNotIn("state", fine[1])
