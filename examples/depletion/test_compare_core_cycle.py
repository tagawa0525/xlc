"""compare_core_cycle.py の要約（サイクルのステップの core の燃焼の、燃焼度の解析値との差、同じ EFPD の参照の燃焼との k の差、
I-135、Xe-135、Sm-149 の数密度）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_cycle.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_cycle

NUCLIDES = ["U-235", "I-135", "Xe-135", "Sm-149"]


def core(points, nodes=1, power=40.0, states=("nominal",)):
    """集合体 1 体・軸方向のノード `nodes` のミクロ燃焼の core の履歴。`points` は (EFPD, 日数, 出力の割合, 状態, 燃焼度, k, Xe-135)"""
    return {
        "format": "xlc-core-history",
        "specific_power": power,
        "states": {
            s: {
                "boron": 1300.0,
                "b10": 0.199,
                "moderator": 565.0,
                "fuel": 565.0,
                "banks": {},
            }
            for s in states
        },
        "axial_mesh": [21.5 * z / nodes for z in range(nodes + 1)],
        "assemblies": ["A01"],
        "heavy_metal": [[9000.0 / nodes] * nodes],
        "nuclides": list(NUCLIDES),
        "points": [
            {
                "efpd": e,
                "days": d,
                **({} if p is None else {"power": p}),
                "state": s,
                "burnup": b,
                "k": k,
                "boron": 1300.0,
                "solves": 1,
                "assembly_burnups": [b],
                "burnups": [[b] * nodes],
                "powers": [[1.0 / nodes] * nodes],
                "densities": [
                    [[7e-4, 1e-9 * (z + 1), x * (z + 1), 3e-8] for z in range(nodes)]
                ],
            }
            for e, d, p, s, b, k, x in points
        ],
    }


CYCLE = [
    (0.0, 0.0, None, "nominal", 0.0, 1.20, 0.0),
    (10.0, 20.0, 0.5, "nominal", 0.4 * (1 + 1e-13), 1.18, 2e-9),
    (10.0, 30.0, 0.0, "nominal", 0.4 * (1 + 1e-13), 1.21, 1e-12),
    (20.0, 40.0, 1.0, "nominal", 0.8, 1.17, 3e-9),
]

REFERENCE = [
    (0.0, 0.0, None, "nominal", 0.0, 1.20, 0.0),
    (10.0, 10.0, 1.0, "nominal", 0.4, 1.17, 3e-9),
    (20.0, 20.0, 1.0, "nominal", 0.8, 1.16, 3e-9),
]


def summarize(cycle=None, references=None):
    return compare_core_cycle.summarize(
        core(CYCLE) if cycle is None else cycle,
        {"full": core(REFERENCE)} if references is None else references,
    )


class Summary(unittest.TestCase):
    def test_burnups_are_compared_with_the_specific_power_times_efpd(self):
        s = summarize()
        p = s["points"]
        self.assertEqual([q["analytic_burnup"] for q in p], [0.0, 0.4, 0.4, 0.8])
        self.assertEqual(p[0]["burnup_relative"], 0.0)
        self.assertAlmostEqual(p[1]["burnup_relative"], 1e-13, delta=1e-15)
        self.assertEqual([q["power"] for q in p], [None, 0.5, 0.0, 1.0])
        self.assertEqual([q["days"] for q in p], [0.0, 20.0, 30.0, 40.0])
        self.assertEqual(s["nuclides"], ["I-135", "Xe-135", "Sm-149"])
        self.assertEqual(
            p[1]["densities"], {"I-135": [1e-9], "Xe-135": [2e-9], "Sm-149": [3e-8]}
        )

    def test_references_pair_with_the_points_at_power_of_the_same_efpd(self):
        # 崩壊だけのステップの点は、同じ EFPD の参照の点と比べない
        p = summarize()["points"]
        full = [q["references"]["full"] for q in p]
        self.assertIsNone(full[2])
        self.assertAlmostEqual(full[1]["rho_pcm"], (1 / 1.17 - 1 / 1.18) * 1e5)
        self.assertEqual(full[1]["k"], 1.17)
        self.assertEqual(full[1]["densities"]["Xe-135"], [3e-9])
        self.assertEqual(full[1]["node_burnups"], [0.4])
        self.assertEqual(full[0]["rho_pcm"], 0.0)

    def test_points_missing_from_a_reference_have_no_pair(self):
        reference = core(REFERENCE[:2])
        p = summarize(references={"full": reference})["points"]
        self.assertIsNone(p[3]["references"]["full"])

    def test_nodes_keep_their_densities(self):
        cycle = core(CYCLE, nodes=2)
        s = summarize(cycle, {"full": core(REFERENCE, nodes=2)})
        self.assertEqual(s["points"][3]["densities"]["Xe-135"], [3e-9, 6e-9])
        self.assertEqual(s["points"][3]["node_burnups"], [0.8, 0.8])

    def test_cores_of_other_shapes_or_powers_are_errors(self):
        two = core(CYCLE)
        two["assemblies"].append("A02")
        for cycle, references in [
            (two, None),
            (None, {"full": core(REFERENCE, power=38.41)}),
            (None, {"full": core(REFERENCE, nodes=2)}),
        ]:
            with self.assertRaises(SystemExit):
                summarize(cycle, references)

    def test_cores_must_track_the_watched_nuclides(self):
        cycle = core(CYCLE)
        cycle["nuclides"] = ["U-235", "I-135", "Xe-135", "Sm-150"]
        with self.assertRaises(SystemExit):
            summarize(cycle)
        cycle = core(CYCLE)
        del cycle["nuclides"]
        with self.assertRaises(SystemExit):
            summarize(cycle)

    def test_reports_give_the_largest_burnup_difference_and_the_reactivity_ranges(self):
        summary = {"inputs": [{"path": "a", "md5": "0" * 32}], **summarize()}
        text = "\n".join(compare_core_cycle.report(summary))
        self.assertIn("1.0e-13", text)
        self.assertIn("full", text)


if __name__ == "__main__":
    unittest.main()
