"""compare_core_partial_rods.py の要約（挿入の割合を変えた部分挿入のノードのミクロ燃焼と、棒のない・棒の入った燃焼の比較）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_partial_rods.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_partial_rods

NUCLIDES = ["U-235", "Xe-135"]
# 下のノード（棒のないノード）の数密度。どの燃焼でも同じにする
LOWER = [7.0e-4, 1.0e-8]


def core(points, banks=None):
    """集合体 1 体・軸方向のノード 2 つのミクロ燃焼の core の履歴。`points` は (k, 上のノードの数密度) で、燃焼度 0、1、2 の点"""
    return {
        "format": "xlc-core-history",
        "specific_power": 40.0,
        "state": {
            "boron": 1300.0,
            "b10": 0.199,
            "moderator": 565.0,
            "fuel": 565.0,
            "banks": banks or {},
        },
        "axial_mesh": [0.0, 10.75, 21.5],
        "assemblies": ["A01"],
        "heavy_metal": [[4500.0, 4500.0]],
        "nuclides": list(NUCLIDES),
        "points": [
            {
                "efpd": 25.0 * b,
                "burnup": float(b),
                "k": k,
                "densities": [[list(LOWER), n]],
            }
            for b, (k, n) in enumerate(points)
        ],
    }


def inputs():
    """棒のない燃焼、棒の入った燃焼、割合 0、0.5、1 の燃焼"""
    out = core(
        [(1.2, [7.0e-4, 0.0]), (1.15, [6.5e-4, 1.0e-8]), (1.1, [6.0e-4, 2.0e-8])]
    )
    rodded = core(
        [(0.9, [7.0e-4, 0.0]), (0.88, [6.8e-4, 0.5e-8]), (0.86, [6.6e-4, 1.0e-8])]
    )
    partial = [
        (
            0.0,
            core(
                [
                    (1.2, [7.0e-4, 0.0]),
                    (1.1501, [6.5e-4, 1.0e-8]),
                    (1.1, [6.0e-4 * (1 + 1e-6), 2.0e-8]),
                ],
                {"D": 172},
            ),
        ),
        (
            0.5,
            core(
                [
                    (1.0, [7.0e-4, 0.0]),
                    (1.0, [6.6e-4, 0.8e-8]),
                    (0.98, [6.3e-4, 1.6e-8]),
                ],
                {"D": 86},
            ),
        ),
        (
            1.0,
            core(
                [
                    (0.9, [7.0e-4, 0.0]),
                    (0.88, [6.8e-4, 0.5e-8]),
                    (0.86, [6.6e-4, 1.0e-8]),
                ],
                {"D": 0},
            ),
        ),
    ]
    return out, rodded, partial


class Summary(unittest.TestCase):
    def test_limits_compare_the_end_fractions_with_the_uniform_nodes(self):
        s = compare_core_partial_rods.summarize(*inputs())
        self.assertEqual(s["fractions"], [0.0, 0.5, 1.0])
        self.assertEqual(s["banks"], [{"D": 172}, {"D": 86}, {"D": 0}])
        p = s["points"][1]
        self.assertEqual(p["k"], [1.1501, 1.0, 0.88])
        self.assertAlmostEqual(p["zero_pcm"], (1 / 1.15 - 1 / 1.1501) * 1e5, places=9)
        self.assertEqual(p["one_pcm"], 0.0)
        u235 = s["densities"]["U-235"]
        self.assertAlmostEqual(u235["zero_limit"]["relative"], 1e-6, places=15)
        self.assertEqual(u235["zero_limit"]["burnup"], 2.0)
        self.assertEqual(u235["zero_limit"]["node"], 1)
        self.assertEqual(u235["one_limit"]["relative"], 0.0)
        # 最後の点の割合 1 の割合 0 に対する相対差は、軸方向のノードごと
        lower, upper = u235["last_change"]
        self.assertEqual(lower, 0.0)
        self.assertAlmostEqual(upper, 6.6 / (6.0 * (1 + 1e-6)) - 1, places=14)
        # 基準が 0 の点（燃焼度 0 の Xe-135）は比べない
        xe = s["densities"]["Xe-135"]
        self.assertEqual(xe["zero_limit"]["relative"], 0.0)
        self.assertAlmostEqual(xe["last_change"][1], -0.5, places=14)

    def test_monotonicity_counts_equal_neighbors_on_both_sides(self):
        s = compare_core_partial_rods.summarize(*inputs())
        # 燃焼度 0 の k は 1.2、1.0、0.9 で減り、数密度はすべて等しい
        self.assertEqual([p["k_monotonic"] for p in s["points"]], [True, True, True])
        self.assertEqual(s["densities"]["U-235"]["not_monotonic"], [])
        out, rodded, partial = inputs()
        # 割合 0.5 の点 1 の Xe-135 を両端より大きくすると、単調でない
        partial[1][1]["points"][1]["densities"][0][1][1] = 2.0e-8
        partial[1][1]["points"][2]["k"] = 1.2
        s = compare_core_partial_rods.summarize(out, rodded, partial)
        self.assertEqual(
            s["densities"]["Xe-135"]["not_monotonic"], [{"burnup": 1.0, "node": 1}]
        )
        self.assertEqual([p["k_monotonic"] for p in s["points"]], [True, True, False])

    def test_fractions_must_run_from_zero_to_one(self):
        out, rodded, partial = inputs()
        for bad in (partial[1:], partial[:-1], [partial[1], partial[0], partial[2]]):
            with self.assertRaises(SystemExit):
                compare_core_partial_rods.summarize(out, rodded, bad)

    def test_histories_must_differ_only_in_the_banks(self):
        for change in (
            lambda c: c["state"].update(boron=1000.0),
            lambda c: c.update(specific_power=38.0),
            lambda c: c.update(nuclides=["U-235", "Xe-136"]),
            lambda c: c["points"][1].update(burnup=1.5),
            lambda c: c.update(axial_mesh=[0.0, 10.0, 21.5]),
            lambda c: c.pop("nuclides"),
            lambda c: c.update(assemblies=["A01", "A02"]),
        ):
            out, rodded, partial = inputs()
            change(partial[1][1])
            with self.assertRaises(SystemExit):
                compare_core_partial_rods.summarize(out, rodded, partial)


if __name__ == "__main__":
    unittest.main()
