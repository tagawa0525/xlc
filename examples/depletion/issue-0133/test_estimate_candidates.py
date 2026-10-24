"""estimate_candidates.py（#0133 の候補の k への寄与の見積もり）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133 -p 'test_estimate_candidates.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import estimate_candidates as ec

NG = 69


def scores(capture: float, fission: float, nu: float) -> dict:
    return {
        "(n,gamma)": {"mean": [capture] + [0.0] * (NG - 1)},
        "fission": {"mean": [fission] + [0.0] * (NG - 1)},
        "nu-fission": {"mean": [nu] + [0.0] * (NG - 1)},
    }


# 生成 P = 1、吸収 A = 1、k = 1。U-238 の捕獲 0.2、Pu-239 の捕獲 0.1
RATES = {
    "k": 1.0,
    "rates": {"U-238": scores(0.2, 0.0, 0.0), "Pu-239": scores(0.1, 0.0, 1.0)},
    "materials": {"fuel": {"absorption": [1.0] + [0.0] * (NG - 1)}},
}


class Estimate(unittest.TestCase):
    def test_bounds_take_each_part_at_its_own_extreme(self):
        # U-238 の捕獲 +1〜+2% は k を −200〜−400 pcm、Pu-239 の捕獲 −2〜−1% は k を +100〜+200 pcm 動かす。和の範囲は
        # −300〜0 pcm（それぞれの端を同じ向きにそろえた −400 + 200 と −200 + 100 では狭すぎる）
        candidates = [
            (
                "c",
                [
                    ("U-238", 1, 1, "capture", 0.01, 0.02),
                    ("Pu-239", 1, 1, "capture", -0.02, -0.01),
                ],
            )
        ]
        ((_, low, high, _),) = ec.estimate(RATES, candidates)
        self.assertAlmostEqual(low, -300.0, places=6)
        self.assertAlmostEqual(high, 0.0, places=6)

    def test_a_nuclide_not_in_the_composition_is_reported(self):
        candidates = [
            (
                "c",
                [
                    ("U-238", 1, 1, "capture", 0.01, 0.01),
                    ("Pm-149", 1, 1, "capture", 0.01, 0.01),
                ],
            )
        ]
        ((_, low, _, missing),) = ec.estimate(RATES, candidates)
        self.assertEqual(missing, ["Pm-149"])
        self.assertAlmostEqual(low, -200.0, places=6)


if __name__ == "__main__":
    unittest.main()
