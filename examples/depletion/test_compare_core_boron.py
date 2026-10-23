"""compare_core_boron.py の要約（core の燃焼の点ごとの臨界ホウ素と、同じ燃焼度の一様な状態の探索）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_boron.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_boron

SEARCH = {"boron": [0.0, 1300.0], "k": 1.19, "tolerance": 1e-6, "max_solves": 20}


def history() -> dict:
    """集合体 1 体の探索の core の履歴（燃焼度 0、0.1、0.5 の点）"""
    return {
        "format": "xlc-core-history",
        "specific_power": 38.41,
        "state": {
            "search": dict(SEARCH),
            "b10": 0.199,
            "moderator": 565.0,
            "fuel": 565.0,
            "banks": {},
        },
        "axial_mesh": [0.0, 21.5],
        "assemblies": ["A01"],
        "heavy_metal": [[27049.1]],
        "points": [
            {
                "efpd": e,
                "burnup": b,
                "k": k,
                "boron": ppm,
                "solves": n,
                "assembly_burnups": [b],
                "burnups": [[b]],
                "powers": [[1.0]],
            }
            for e, b, k, ppm, n in [
                (0.0, 0.0, 1.19 + 4e-7, 1234.567891, 4),
                (2.6, 0.1, 1.19 - 2e-7, 851.234567, 5),
                (13.0, 0.5, 1.19, 812.000004, 3),
            ]
        ],
    }


# 同じ燃焼度を状態に与えた探索のケース（名前の順に解く）と、その標準出力
UNIFORM_CASE = """
[state.p0]
b10 = 0.199
moderator = 565
fuel = 565
burnup = 0.0
search = { boron = [0.0, 1300.0], k = 1.19, tolerance = 1e-6, max-solves = 20 }

[state.p1]
b10 = 0.199
moderator = 565
fuel = 565
burnup = 0.1
search = { boron = [0.0, 1300.0], k = 1.19, tolerance = 1e-6, max-solves = 20 }

[state.p2]
b10 = 0.199
moderator = 565
fuel = 565
burnup = 0.5
search = { boron = [0.0, 1300.0], k = 1.19, tolerance = 1e-6, max-solves = 20 }
"""

UNIFORM_OUT = """state p0
boron 1234.5679 ppm (searched in 4 solves)
k-eff 1.19000040
leakage per fission neutron 0.0000
state p1
boron 851.2346 ppm (searched in 5 solves)
k-eff 1.18999980
leakage per fission neutron 0.0000
state p2
boron 812.0000 ppm (searched in 3 solves)
k-eff 1.19000000
leakage per fission neutron 0.0000
"""


def summarize(core=None, case=UNIFORM_CASE, out=UNIFORM_OUT):
    return compare_core_boron.summarize(history() if core is None else core, case, out)


class Summary(unittest.TestCase):
    def test_points_have_the_boron_and_the_distance_from_the_target(self):
        s = summarize()
        self.assertEqual(s["search"], SEARCH)
        p = s["points"]
        self.assertEqual([q["burnup"] for q in p], [0.0, 0.1, 0.5])
        self.assertEqual([q["boron"] for q in p], [1234.567891, 851.234567, 812.000004])
        self.assertEqual([q["solves"] for q in p], [4, 5, 3])
        self.assertAlmostEqual(p[0]["k_minus_target"], 4e-7, places=15)
        self.assertAlmostEqual(s["max_abs_k_minus_target"], 4e-7, places=15)
        self.assertTrue(s["within_tolerance"])

    def test_the_boron_must_fall_with_burnup(self):
        self.assertTrue(summarize()["boron_falls"])
        core = history()
        core["points"][2]["boron"] = 900.0
        self.assertFalse(summarize(core)["boron_falls"])

    def test_points_pair_with_the_uniform_state_searches_by_burnup(self):
        s = summarize()
        p = s["points"]
        self.assertEqual([q["uniform_state"] for q in p], ["p0", "p1", "p2"])
        self.assertEqual([q["uniform_boron"] for q in p], [1234.5679, 851.2346, 812.0])
        self.assertEqual([q["uniform_solves"] for q in p], [4, 5, 3])
        self.assertEqual([q["uniform_k"] for q in p], [1.1900004, 1.1899998, 1.19])
        # 一様な状態の出力は 4 桁なので、点のホウ素を同じ桁に丸めて比べる
        self.assertEqual([q["boron_rounded_difference"] for q in p], [0.0, 0.0, 0.0])
        self.assertTrue(s["uniform_agrees"])

    def test_a_differing_uniform_search_is_reported(self):
        s = summarize(out=UNIFORM_OUT.replace("851.2346", "851.2400"))
        self.assertAlmostEqual(
            s["points"][1]["boron_rounded_difference"], -0.0054, places=9
        )
        self.assertFalse(s["uniform_agrees"])
        s = summarize(out=UNIFORM_OUT.replace("in 5 solves", "in 6 solves"))
        self.assertFalse(s["uniform_agrees"])

    def test_points_without_a_uniform_state_are_errors(self):
        with self.assertRaises(SystemExit):
            summarize(case=UNIFORM_CASE.replace("burnup = 0.5", "burnup = 0.6"))

    def test_the_history_must_search_on_one_node(self):
        core = history()
        del core["state"]["search"]
        core["state"]["boron"] = 1300.0
        with self.assertRaises(SystemExit):
            summarize(core)
        core = history()
        core["assemblies"].append("A02")
        with self.assertRaises(SystemExit):
            summarize(core)

    def test_the_uniform_states_must_search_with_the_same_settings(self):
        with self.assertRaises(SystemExit):
            summarize(case=UNIFORM_CASE.replace("k = 1.19", "k = 1.0", 1))

    def test_the_report_has_the_checks(self):
        s = summarize()
        s["inputs"] = [{"path": "core.msgpack", "md5": "0" * 32}]
        text = "\n".join(compare_core_boron.report(s))
        self.assertIn("|k − 目標|", text)
        self.assertIn("単調", text)


if __name__ == "__main__":
    unittest.main()
