"""compare_core_micro.py の要約（core のミクロ燃焼の数密度のテーブルとの相対差と、k の履歴との差）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_micro.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_micro

STATE = {"boron": 1300.0, "b10": 0.199, "moderator": 565.0, "fuel": 565.0, "banks": {}}


def core(points, nuclides=None):
    """集合体 1 体・ノード 1 つの core の履歴。`points` は (EFPD, 燃焼度, k, 数密度)"""
    history = {
        "format": "xlc-core-history",
        "specific_power": 40.0,
        "state": dict(STATE),
        "axial_mesh": [0.0, 21.5],
        "assemblies": ["A01"],
        "heavy_metal": [[9000.0]],
        "points": [],
    }
    for e, b, k, n in points:
        point = {"efpd": e, "burnup": b, "k": k, "boron": 1300.0}
        if nuclides is not None:
            point["densities"] = [[n]]
        history["points"].append(point)
    if nuclides is not None:
        history["nuclides"] = nuclides
    return history


def inputs():
    """格子計算の履歴（燃焼度 0、1、2）、テーブル（燃焼度 0、2 の数密度）、ミクロ燃焼の core とミクロ燃焼しない core（0、1、2 の点）"""
    lattice = {
        "format": "xlc-depletion-history",
        "points": [
            {"burnup": b, "k_inf": k} for b, k in [(0.0, 1.2), (1.0, 1.15), (2.0, 1.1)]
        ],
    }
    table = {
        "format": "xlc-table",
        "axes": [
            {"name": "burnup", "points": [0.0, 2.0]},
            {"name": "boron", "points": [1300.0], "b10": 0.199},
        ],
        "micro": {
            "nuclides": ["U-235", "Xe-135"],
            "densities": [[7.0e-4, 0.0], [6.0e-4, 1.0e-8]],
        },
    }
    micro = core(
        [
            (0.0, 0.0, 1.2, [7.0e-4, 0.0]),
            (25.0, 1.0, 1.151, [6.5e-4, 0.9e-8]),
            (50.0, 2.0 + 1e-12, 1.1002, [6.006e-4, 0.99e-8]),
        ],
        ["U-235", "Xe-135"],
    )
    plain = core(
        [(0.0, 0.0, 1.2, None), (25.0, 1.0, 1.152, None), (50.0, 2.0, 1.1, None)]
    )
    return lattice, table, micro, plain


class Summary(unittest.TestCase):
    def test_densities_are_compared_with_the_table_at_its_burnup_points(self):
        s = compare_core_micro.summarize(*inputs())
        self.assertEqual([d["burnup"] for d in s["density_points"]], [0.0, 2.0])
        first, last = s["density_points"]
        # テーブルの数密度が 0 の核種は比べない
        self.assertEqual(first["relative"], {"U-235": 0.0, "Xe-135": None})
        self.assertAlmostEqual(
            last["relative"]["U-235"], 6.006e-4 / 6.0e-4 - 1, places=14
        )
        self.assertAlmostEqual(last["relative"]["Xe-135"], 0.99 - 1, places=14)
        # 核種ごとの絶対値の最大と、その燃焼度
        self.assertEqual(s["density_largest"]["Xe-135"]["burnup"], 2.0)
        self.assertAlmostEqual(
            s["density_largest"]["Xe-135"]["relative"], -0.01, places=14
        )
        self.assertEqual(s["density_largest"]["U-235"]["burnup"], 2.0)

    def test_k_differences_separate_the_correction_from_the_interpolation(self):
        p = compare_core_micro.summarize(*inputs())["points"]
        self.assertEqual([q["on_table"] for q in p], [True, False, True])
        # Δρ = (1/k_基準 − 1/k) × 1e5 pcm
        self.assertAlmostEqual(
            p[1]["micro_pcm"], (1 / 1.15 - 1 / 1.151) * 1e5, places=9
        )
        self.assertAlmostEqual(
            p[1]["plain_pcm"], (1 / 1.15 - 1 / 1.152) * 1e5, places=9
        )
        self.assertAlmostEqual(
            p[1]["correction_pcm"], (1 / 1.152 - 1 / 1.151) * 1e5, places=9
        )
        self.assertEqual(p[2]["plain_pcm"], 0.0)
        self.assertEqual(p[0]["correction_pcm"], 0.0)

    def test_the_core_histories_must_be_given_in_order_and_agree(self):
        lattice, table, micro, plain = inputs()
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, plain, micro)
        plain["state"]["boron"] = 1000.0
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, micro, plain)
        lattice, table, micro, plain = inputs()
        plain["points"][1]["burnup"] = 1.5
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, micro, plain)

    def test_burnups_off_the_history_and_other_nuclides_are_errors(self):
        lattice, table, micro, plain = inputs()
        micro["points"][1]["burnup"] = 1.1
        plain["points"][1]["burnup"] = 1.1
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, micro, plain)
        lattice, table, micro, plain = inputs()
        micro["nuclides"] = ["U-235", "Xe-136"]
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, micro, plain)

    def test_the_core_must_be_one_node(self):
        lattice, table, micro, plain = inputs()
        micro["axial_mesh"] = [0.0, 10.0, 21.5]
        with self.assertRaises(SystemExit):
            compare_core_micro.summarize(lattice, table, micro, plain)

    def test_the_report_has_the_differences(self):
        s = compare_core_micro.summarize(*inputs())
        s["inputs"] = [{"path": "core.msgpack", "md5": "0" * 32}]
        text = "\n".join(compare_core_micro.report(s))
        self.assertIn("ミクロ燃焼 − ミクロ燃焼なし", text)
        # 点の間だけの範囲（テーブルの点の値を混ぜない）
        self.assertIn("ミクロ燃焼 − 履歴（点の間）", text)
        self.assertIn("Xe-135 | -1.0000 | 2.0", text)


if __name__ == "__main__":
    unittest.main()
