"""compare_core_depletion.py の要約（core の燃焼の燃焼点ごとの k と、格子計算の燃焼の履歴の k の差）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_depletion.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_depletion


def inputs():
    """格子計算の履歴（燃焼度 0、1、2、3）、テーブル（燃焼度 0、2）、集合体 1 体の core の履歴（0、1、2 の点）"""
    lattice = {
        "format": "xlc-depletion-history",
        "points": [
            {"burnup": b, "k_inf": k}
            for b, k in [(0.0, 1.2), (1.0, 1.15), (2.0, 1.1), (3.0, 1.05)]
        ],
    }
    # 燃焼度 2 点 × ホウ素 2 点 × 減速材 1 点 × 燃料 1 点。公称（1300 ppm）の分岐の k∞ は 1.2 + 2e-8 と 1.1 − 1e-8
    table = {
        "format": "xlc-table",
        "split": 1,
        "axes": [
            {"name": "burnup", "points": [0.0, 2.0]},
            {"name": "boron", "points": [1000.0, 1300.0], "b10": 0.199},
            {"name": "moderator", "points": [565.0], "density": [0.743]},
            {"name": "fuel", "points": [565.0]},
        ],
        "states": [{"k_inf": k} for k in [1.3, 1.2 + 2e-8, 1.2, 1.1 - 1e-8]],
    }
    core = {
        "format": "xlc-core-history",
        "specific_power": 40.0,
        "state": {
            "boron": 1300.0,
            "b10": 0.199,
            "moderator": 565.0,
            "fuel": 565.0,
            "banks": {},
        },
        "axial_mesh": [0.0, 21.5],
        "assemblies": ["A01"],
        "heavy_metal": [[9000.0]],
        "points": [
            {
                "efpd": e,
                "burnup": b,
                "k": k,
                "assembly_burnups": [b],
                "burnups": [[b]],
                "powers": [[1.0]],
            }
            # 点の上は履歴と 1e-8 の差、点の間（1.0）は内挿の誤差
            for e, b, k in [
                (0.0, 0.0, 1.2 + 1e-8),
                (25.0, 1.0, 1.152),
                (50.0, 2.0, 1.1),
            ]
        ],
    }
    return lattice, table, core


class Summary(unittest.TestCase):
    def test_points_pair_with_the_history_by_burnup(self):
        s = compare_core_depletion.summarize(*inputs())
        self.assertEqual(s["table_burnups"], [0.0, 2.0])
        self.assertEqual([p["burnup"] for p in s["points"]], [0.0, 1.0, 2.0])
        self.assertEqual([p["on_table"] for p in s["points"]], [True, False, True])
        self.assertEqual([p["history_k"] for p in s["points"]], [1.2, 1.15, 1.1])
        self.assertEqual(s["points"][1]["efpd"], 25.0)

    def test_reactivity_differences_are_core_minus_history(self):
        s = compare_core_depletion.summarize(*inputs())
        p = s["points"]
        # Δρ = (1/k_履歴 − 1/k_core) × 1e5 pcm
        self.assertAlmostEqual(p[1]["rho_pcm"], (1 / 1.15 - 1 / 1.152) * 1e5, places=9)
        self.assertAlmostEqual(
            p[0]["rho_pcm"], (1 / 1.2 - 1 / (1.2 + 1e-8)) * 1e5, places=12
        )
        self.assertEqual(p[2]["rho_pcm"], 0.0)
        self.assertAlmostEqual(
            s["on_table_max_abs_pcm"], abs(p[0]["rho_pcm"]), places=15
        )
        self.assertAlmostEqual(
            s["between_max_abs_pcm"], abs(p[1]["rho_pcm"]), places=12
        )

    def test_table_points_compare_with_the_branch_lattice_k(self):
        s = compare_core_depletion.summarize(*inputs())
        p = s["points"]
        # 分岐の格子計算の k∞ は、core の状態のテーブルの状態の k∞
        self.assertEqual([q["branch_k"] for q in p], [1.2 + 2e-8, None, 1.1 - 1e-8])
        self.assertAlmostEqual(
            p[0]["branch_rho_pcm"],
            (1 / (1.2 + 2e-8) - 1 / (1.2 + 1e-8)) * 1e5,
            places=12,
        )
        self.assertIsNone(p[1]["branch_rho_pcm"])
        self.assertAlmostEqual(
            s["on_table_branch_max_abs_pcm"],
            max(abs(p[0]["branch_rho_pcm"]), abs(p[2]["branch_rho_pcm"])),
            places=15,
        )

    def test_a_state_off_the_table_points_is_an_error(self):
        lattice, table, core = inputs()
        core["state"]["boron"] = 1200.0
        with self.assertRaises(SystemExit):
            compare_core_depletion.summarize(lattice, table, core)

    def test_histories_that_search_the_boron_are_errors(self):
        # 探索した履歴は点ごとにホウ素が違い、状態の 1 つのホウ素でテーブルを引けない（compare_core_boron.py で比べる）
        lattice, table, core = inputs()
        del core["state"]["boron"]
        core["state"]["search"] = {
            "boron": [0.0, 1300.0],
            "k": 1.0,
            "tolerance": 1e-6,
            "max_solves": 20,
        }
        with self.assertRaises(SystemExit):
            compare_core_depletion.summarize(lattice, table, core)

    def test_states_without_a_boron_or_a_search_are_errors(self):
        lattice, table, core = inputs()
        del core["state"]["boron"]
        with self.assertRaisesRegex(SystemExit, "neither"):
            compare_core_depletion.summarize(lattice, table, core)

    def test_burnups_off_the_history_points_are_errors(self):
        lattice, table, core = inputs()
        core["points"][1]["burnup"] = 1.1
        with self.assertRaises(SystemExit):
            compare_core_depletion.summarize(lattice, table, core)

    def test_the_core_must_be_one_node(self):
        lattice, table, core = inputs()
        core["axial_mesh"] = [0.0, 10.0, 21.5]
        with self.assertRaises(SystemExit):
            compare_core_depletion.summarize(lattice, table, core)
        lattice, table, core = inputs()
        core["assemblies"].append("A02")
        with self.assertRaises(SystemExit):
            compare_core_depletion.summarize(lattice, table, core)

    def test_the_report_has_the_differences(self):
        s = compare_core_depletion.summarize(*inputs())
        s["inputs"] = [{"path": "core.msgpack", "md5": "0" * 32}]
        text = "\n".join(compare_core_depletion.report(s))
        self.assertIn("テーブルの点", text)
        self.assertIn("点の間", text)


if __name__ == "__main__":
    unittest.main()
