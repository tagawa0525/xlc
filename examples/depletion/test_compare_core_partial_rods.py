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


def core(points, banks=None, shift=0.0):
    """集合体 1 体・軸方向のノード 2 つのミクロ燃焼の core の履歴。`points` は (k, 上のノードの数密度) で、燃焼度 0、1、2 の点。
    ノードの燃焼度は、下が集合体の (1 + `shift`) 倍、上が (1 − `shift`) 倍"""
    return {
        "format": "xlc-core-history",
        "specific_power": 40.0,
        "states": {
            "nominal": {
                "boron": 1300.0,
                "b10": 0.199,
                "moderator": 565.0,
                "fuel": 565.0,
                "banks": banks or {},
            }
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
                "burnups": [[b * (1 + shift), b * (1 - shift)]],
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
                0.1,
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
                0.2,
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

    def test_node_burnups_of_every_fraction_are_kept_with_their_monotonicity(self):
        s = compare_core_partial_rods.summarize(*inputs())
        last = s["points"][2]
        self.assertEqual(last["node_burnups"], [[2.0, 2.0], [2.2, 1.8], [2.4, 1.6]])
        self.assertTrue(all(p["node_burnups_monotonic"] for p in s["points"]))
        out, rodded, partial = inputs()
        partial[1][1]["points"][2]["burnups"] = [[2.5, 1.5]]
        s = compare_core_partial_rods.summarize(out, rodded, partial)
        self.assertEqual(
            [p["node_burnups_monotonic"] for p in s["points"]], [True, True, False]
        )

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

    def test_change_per_node_burnup_is_checked_for_monotonicity_apart(self):
        # 数密度の燃焼度 0 からの変化をノードの燃焼度で割った値（燃焼度 0 の点は除く）
        s = compare_core_partial_rods.summarize(*inputs())
        for d in s["densities"].values():
            self.assertEqual(d["per_burnup_not_monotonic"], [])
        out, rodded, partial = inputs()
        # 割合 0.5 の点 2 の U-235 を、数密度は両端の間のまま、燃焼度あたりの変化は割合 0 より大きく減るようにする
        partial[1][1]["points"][2]["densities"][0][1][0] = 6.05e-4
        s = compare_core_partial_rods.summarize(out, rodded, partial)
        u235 = s["densities"]["U-235"]
        self.assertEqual(u235["not_monotonic"], [])
        self.assertEqual(u235["per_burnup_not_monotonic"], [{"burnup": 2.0, "node": 1}])

    def test_burnups_of_the_points_agree_to_rounding(self):
        # 集合体の燃焼度はノードの燃焼度の重金属の質量の重みの平均なので、燃焼によって丸めが違う
        out, rodded, partial = inputs()
        rodded["points"][1]["burnup"] = 1.0 + 2.2e-16
        s = compare_core_partial_rods.summarize(out, rodded, partial)
        self.assertEqual([p["burnup"] for p in s["points"]], [0.0, 1.0, 2.0])

    def test_banks_must_move_one_way_with_the_fraction(self):
        # 履歴と割合の組を取り違えると、バンクの位置が割合について単調でなくなる
        out, rodded, partial = inputs()
        swapped = [(0.0, partial[0][1]), (0.5, partial[2][1]), (1.0, partial[1][1])]
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, swapped)
        # 割合を逆に付けると、ステップが割合について増える（棒は上から入り、ステップが少ないほど深い）
        out, rodded, partial = inputs()
        reversed_ = [(f, c) for (f, _), (_, c) in zip(partial, reversed(partial))]
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, reversed_)
        # 順を保った取り違え: 割合 0.5 に、ステップが 0 と 1 の間の 1/4 の燃焼を付ける
        out, rodded, partial = inputs()
        partial[1][1]["states"]["nominal"]["banks"] = {"D": 129}
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, partial)
        # どの割合の燃焼もバンクを持たなければ、取り違えを確かめられない
        out, rodded, partial = inputs()
        for _, c in partial:
            c["states"]["nominal"]["banks"] = {}
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, partial)
        # 棒のない燃焼と棒の入った燃焼を取り違えると、棒の入った燃焼の最初の点の k が大きい
        out, rodded, partial = inputs()
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(rodded, out, partial)
        # 比べる先の燃焼はバンクを持たない
        out, rodded, partial = inputs()
        rodded["states"]["nominal"]["banks"] = {"D": 0}
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, partial)

    def test_no_fractions_is_an_error(self):
        out, rodded, _ = inputs()
        with self.assertRaises(SystemExit):
            compare_core_partial_rods.summarize(out, rodded, [])

    def test_the_report_of_a_summary_without_comparable_densities_is_drawn(self):
        out, rodded, partial = inputs()
        for core in [out, rodded] + [c for _, c in partial]:
            for p in core["points"]:
                p["densities"] = [[[0.0, 0.0], [0.0, 0.0]]]
        s = compare_core_partial_rods.summarize(out, rodded, partial)
        s["inputs"] = []
        lines = compare_core_partial_rods.report(s)
        self.assertTrue(any("相対差の絶対値の最大: –" in line for line in lines))

    def test_fractions_must_run_from_zero_to_one(self):
        out, rodded, partial = inputs()
        for bad in (partial[1:], partial[:-1], [partial[1], partial[0], partial[2]]):
            with self.assertRaises(SystemExit):
                compare_core_partial_rods.summarize(out, rodded, bad)

    def test_histories_must_differ_only_in_the_banks(self):
        for change in (
            lambda c: c["states"]["nominal"].update(boron=1000.0),
            lambda c: c.update(specific_power=38.0),
            lambda c: c.update(nuclides=["U-235", "Xe-136"]),
            lambda c: c["points"][1].update(burnup=1.5),
            lambda c: c.update(axial_mesh=[0.0, 10.0, 21.5]),
            lambda c: c.pop("nuclides"),
            lambda c: c.update(assemblies=["A01", "A02"]),
            lambda c: c.update(heavy_metal=[[4500.0, 4400.0]]),
            lambda c: c["points"][1].update(burnups=[[2.0, 0.0]]),
            lambda c: c["points"][1].update(densities=[[[], [6.6e-4, 0.8e-8]]]),
        ):
            out, rodded, partial = inputs()
            change(partial[1][1])
            with self.assertRaises(SystemExit):
                compare_core_partial_rods.summarize(out, rodded, partial)


if __name__ == "__main__":
    unittest.main()
