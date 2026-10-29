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
        for name, text in small_core_cases.restart_files().items():
            path = small_core_cases.BRANCHES / name
            self.assertEqual(path.read_text(), text, name)

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


def step_lines(text: str) -> list[str]:
    """ケースの `steps` の行"""
    body = text.split("steps = [\n")[1].split("\n]")[0]
    return body.splitlines()


def efpd(line: str) -> float | None:
    if "efpd = " not in line:
        return None
    return float(line.split("efpd = ")[1].split(",")[0].split(" ")[0])


class Restarts(unittest.TestCase):
    files = small_core_cases.restart_files()

    def test_the_halves_join_to_the_uncut_steps(self):
        # 切った前と後の半分のステップは、切らない燃焼のステップに戻る（後の半分の EFPD は切った点から数え直す。ADR-0220）
        full = step_lines(self.files["core-small-restart-full.toml"])
        self.assertIn("decay-days = 30.0", full[-1])
        for cut, n in small_core_cases.CUTS.items():
            first = step_lines(self.files[f"core-small-restart-{cut}1.toml"])
            second = step_lines(self.files[f"core-small-restart-{cut}2.toml"])
            self.assertEqual(len(first), n)
            self.assertEqual(len(first) + len(second), len(full))
            start = small_core_cases.cycle_steps(small_core_cases.RESTART_CYCLE)[n - 1][
                "efpd"
            ]
            for got, want in zip(second, full[n:]):
                if efpd(want) is None:
                    self.assertIsNone(efpd(got))
                else:
                    self.assertAlmostEqual(efpd(got) + start, efpd(want), delta=1e-9)
            # 後の半分の最初のステップは、切った点で効いている状態を書く
            self.assertIn("state = ", second[0])

    def test_cases_hold_only_the_states_their_steps_use(self):
        # どのステップも指さない状態は燃焼の入力のエラー（ADR-0218）
        a1 = self.files["core-small-restart-a1.toml"]
        b2 = self.files["core-small-restart-b2.toml"]
        self.assertNotIn("[state.rodded]", a1)
        self.assertNotIn("[state.rodded]", b2)
        self.assertIn('state = "aro"', step_lines(b2)[0])
        self.assertIn("[state.rodded]", self.files["core-small-restart-a2.toml"])

    def test_second_halves_restart_from_the_first_half_histories(self):
        for cut in small_core_cases.CUTS:
            cycle = self.files[f"core-small-restart-{cut}2-cycle.toml"]
            self.assertIn(
                f'restart = {{ history = "core-small-restart-{cut}1.msgpack" }}', cycle
            )
            self.assertIn(
                f'cycle = "core-small-restart-{cut}2-cycle.toml"',
                self.files[f"core-small-restart-{cut}2.toml"],
            )

    def test_the_reload_restarts_from_the_uncut_history_with_a_fresh_assembly(self):
        cycle = self.files["core-small-reload-cycle.toml"]
        self.assertIn(
            'restart = { history = "core-small-restart-full.msgpack", fresh = ["A05"] }',
            cycle,
        )
        self.assertIn('"A04  A02@90  ."', cycle)
        self.assertIn('"A05  A01     ."', cycle)
        self.assertNotIn("A03", cycle.split("loading")[1])
        # サイクルの EFPD は 0 から、燃焼度の増分 3 MWd/kg まで
        lines = step_lines(self.files["core-small-reload.toml"])
        self.assertAlmostEqual(efpd(lines[-1]), 3000 / 38.41, delta=1e-9)
