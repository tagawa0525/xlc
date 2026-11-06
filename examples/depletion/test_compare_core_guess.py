"""compare_core_guess.py の要約（固有値計算の初期値だけを変えた core の燃焼の差と、固有値計算の数、反復の回数、時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_guess.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_guess


def core(k=1.0, boron=1000.0, burnup=1.0, power=0.5, solves=(1, 2)):
    """集合体 2 体・軸方向のノード 1 つの燃焼の履歴。点は 0 と 10 EFPD。2 体目のノードは核分裂がない"""
    return {
        "specific_power": 40.0,
        "states": {"hfp": {}},
        "axial_mesh": [0.0, 10.0],
        "assemblies": ["A", "R"],
        "heavy_metal": [[1.0], [0.0]],
        "points": [
            {
                "efpd": efpd,
                "days": efpd,
                "state": "hfp",
                "k": k,
                "boron": boron,
                "solves": n,
                "burnups": [[burnup * efpd], [0.0]],
                "powers": [[power], [0.0]],
            }
            for efpd, n in zip((0.0, 10.0), solves)
        ],
    }


def log(iterations, midpoints=((1.0, 900.0, 1),), total=5.0):
    """固有値計算ごとの反復の回数 `iterations`（1 回 0.25 s）、中点の行 (k, ホウ素, 回数)、全体の時間の core のログ"""
    lines = [
        f"2026-11-06T00:00:00Z  INFO nem{{iterations={i} seconds=0.25}}: xlc_core::nem: close time.busy=250ms time.idle=1µs"
        for i in iterations
    ]
    lines += [
        f"2026-11-06T00:00:00Z  INFO depletion step{{efpd=10.0}}: xlc_core::depletion: depletion midpoint k={k} boron={b} solves={n}"
        for k, b, n in midpoints
    ]
    lines.append(f"2026-11-06T00:00:00Z  INFO core: core seconds={total}")
    return "\n".join(lines) + "\n"


class Compare(unittest.TestCase):
    def test_differences_are_the_default_minus_the_flat_start(self):
        flat = core()
        default = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4, burnup=1.0 + 1e-9)
        c = compare_core_guess.compare(
            flat,
            log([100, 100, 100, 100], midpoints=((1.0, 900.0, 1),)),
            default,
            log(
                [100, 40, 30, 20], midpoints=((1.0 - 1e-9, 900.0 + 2e-4, 1),), total=4.0
            ),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(later["burnup"], 1e-9, delta=1e-15)
        # 核分裂のないノードは除く
        self.assertEqual(later["power"], 0.0)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(
            c["runs"]["default"],
            {
                "solves": 4,
                "iterations": 190,
                "eigenvalue_seconds": 1.0,
                "seconds": 4.0,
            },
        )
        self.assertEqual(c["runs"]["flat"]["iterations"], 400)

    def test_solve_counts_must_match_the_nem_spans(self):
        # 点 1 + 2 と中点 1 の 4 回に、span が 3 つしかない
        with self.assertRaises(SystemExit):
            compare_core_guess.compare(core(), log([1, 1, 1]), core(), log([1, 1, 1]))

    def test_histories_of_other_points_are_errors(self):
        other = core()
        other["points"][1]["days"] = 11.0
        with self.assertRaises(SystemExit):
            compare_core_guess.compare(core(), log([1] * 4), other, log([1] * 4))

    def test_power_in_one_history_only_is_an_error(self):
        other = core()
        other["points"][1]["powers"][1] = [0.1]
        with self.assertRaises(SystemExit):
            compare_core_guess.compare(core(), log([1] * 4), other, log([1] * 4))


if __name__ == "__main__":
    unittest.main()
