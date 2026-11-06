"""compare_core_threads.py の要約（固有値計算を解く順とスレッドの数だけを変えた core の燃焼の差、履歴のバイトの一致、時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_threads.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_threads
import msgpack
from test_compare_core_guess import core, log


def run(history: dict, iterations, threads: int, total: float, midpoint=900.0):
    """履歴のバイトと、スレッドのプールの行を先頭に置いたログ"""
    text = log(iterations, midpoints=((1.0, midpoint, 1),), total=total)
    pool = (
        f'2026-11-06T00:00:00Z  INFO core: thread pool threads={threads} from="input"\n'
    )
    return msgpack.packb(history), pool + text


class Compare(unittest.TestCase):
    def test_differences_are_the_red_black_minus_the_lexicographic_solves(self):
        lexicographic = run(core(), [60, 20, 15, 40], 1, 2.0)
        red_black = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4)
        runs = [
            run(red_black, [50, 20, 10, 30], 1, 1.8, midpoint=900.0 + 2e-4),
            run(red_black, [50, 20, 10, 30], 4, 0.9, midpoint=900.0 + 2e-4),
        ]
        c = compare_core_threads.compare(lexicographic, runs)
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(
            [(r["order"], r["threads"], r["same_bytes"]) for r in c["runs"]],
            [
                ("lexicographic", 1, None),
                ("red-black", 1, True),
                ("red-black", 4, True),
            ],
        )
        self.assertEqual(
            c["runs"][2] | {"order": None},
            {
                "order": None,
                "threads": 4,
                "same_bytes": True,
                "solves": 4,
                "iterations": 110,
                "eigenvalue_seconds": 1.0,
                "seconds": 0.9,
                "per_solve": {"iterations": [10, 50], "seconds": [0.25, 0.25]},
            },
        )

    def test_red_black_histories_of_other_bytes_are_marked(self):
        lexicographic = run(core(), [60, 20, 15, 40], 1, 2.0)
        runs = [
            run(core(), [50, 20, 10, 30], 1, 1.8),
            run(core(k=1.0 + 1e-15), [50, 20, 10, 30], 16, 0.5),
        ]
        c = compare_core_threads.compare(lexicographic, runs)
        self.assertEqual([r["same_bytes"] for r in c["runs"]], [None, True, False])

    def test_red_black_runs_need_a_thread_pool_line(self):
        history, text = run(core(), [60, 20, 15, 40], 1, 2.0)
        bare = text.split("\n", 1)[1]
        with self.assertRaises(SystemExit):
            compare_core_threads.compare(
                run(core(), [1, 1, 1, 1], 1, 1.0), [(history, bare)]
            )

    def test_a_base_without_a_thread_pool_line_is_the_bin_core_before(self):
        # ADR-0234 の前の bin/core はスレッドのプールの行を出さない（固有値計算は 1 スレッドの直方体の並びの順）
        history, text = run(core(), [60, 20, 15, 40], 1, 2.0)
        bare = text.split("\n", 1)[1]
        c = compare_core_threads.compare(
            (history, bare), [run(core(), [50, 20, 10, 30], 16, 0.5)]
        )
        self.assertEqual(
            [(r["order"], r["threads"]) for r in c["runs"]],
            [("before", None), ("red-black", 16)],
        )
        text = "\n".join(compare_core_threads.report({"cases": {"p9": c}}))
        self.assertIn("ADR-0234 の前の bin/core: 固有値計算 4 回", text)

    def test_report_gives_the_time_per_thread_count(self):
        c = compare_core_threads.compare(
            run(core(), [60, 20, 15, 40], 1, 2.0),
            [
                run(core(), [50, 20, 10, 30], 1, 1.8),
                run(core(), [50, 20, 10, 30], 4, 0.9),
            ],
        )
        text = "\n".join(compare_core_threads.report({"cases": {"small": c}}))
        self.assertIn(
            "直方体の並びの順 1 スレッド: 固有値計算 4 回、外側の反復 135 回", text
        )
        self.assertIn("色の順 4 スレッド: 固有値計算 4 回、外側の反復 110 回", text)
        self.assertIn("色の順の履歴はスレッドの数によらずバイトで同じ", text)


if __name__ == "__main__":
    unittest.main()
