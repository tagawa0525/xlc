"""compare_rodless.py: 制御棒の分岐を解いた格子の値 out の状態と、前の形の棒のない格子の状態の比べ方。

実行: python3 -m unittest discover -s examples/vera/p5 -p 'test_compare_rodless.py'
"""

import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_rodless

CONTINUOUS = [
    {"name": "boron", "points": [1100.0, 1400.0], "b10": 0.199},
    {"name": "moderator", "points": [560.0, 570.0], "density": [0.75207, 0.73294]},
]
ROD = {"name": "rod", "values": ["out", "aic", "b4c"]}


def state(k: float) -> dict:
    return {
        "k_inf": k,
        "total": [0.5, 1.2],
        "scatter": [[0.4, 0.02], [0.0, 1.0]],
        "pins": None,
    }


def rodless_table(ks: list[float]) -> dict:
    points = [[b, m] for b in range(2) for m in range(2)]
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "axes": CONTINUOUS,
        "split": 1,
        "points": points,
        "states": [state(k) for k in ks],
    }


def rod_table(out: list[float], first: str = "out") -> dict:
    """制御棒の軸を最後に持つテーブル。点は辞書式で、値 aic と b4c の状態は k を下げる"""
    values = [first, "aic", "b4c"]
    points = [[b, m, r] for b in range(2) for m in range(2) for r in range(3)]
    states = [state(out[2 * b + m] - 0.1 * r) for b, m, r in points]
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "axes": [*CONTINUOUS, {"name": "rod", "values": values}],
        "split": 1,
        "points": points,
        "states": states,
    }


def medium(scale: float) -> dict:
    return {
        "format": "xlc-medium",
        "bounds": [1e7, 1.0, 1e-5],
        "axes": CONTINUOUS,
        "states": [
            {
                "total": [scale * i, 2.0],
                "scatter": [[[0.1, 0.0], [0.0, 0.2]]],
                "fission": [0.0, 0.0],
                "nu_fission": [0.0, 0.0],
                "chi": [1.0, 0.0],
            }
            for i in range(4)
        ],
    }


KS = [1.01, 1.0, 0.98, 0.97]


class Compare(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.root = Path(self.dir.name)

    def tearDown(self):
        self.dir.cleanup()

    def write(self, name: str, data: dict) -> str:
        path = self.root / name
        path.write_bytes(msgpack.packb(data))
        return str(path)

    def test_out_states_equal_to_the_rodless_table_are_bit_equal(self):
        result = compare_rodless.compare(
            self.write("rod.msgpack", rod_table(KS)),
            self.write("rodless.msgpack", rodless_table(KS)),
        )
        self.assertEqual(result["kind"], "table")
        self.assertTrue(result["rod_axis"])
        self.assertEqual((result["states"], result["bit_equal"]), (4, 4))
        self.assertEqual(result["max_relative_difference"], 0.0)
        self.assertIsNone(result["field"])
        self.assertEqual(result["max_k_inf_pcm"], 0.0)

    def test_a_state_one_ulp_away_is_counted_with_its_field_and_k(self):
        ks = list(KS)
        ks[2] = math.nextafter(0.98, 1.0)
        result = compare_rodless.compare(
            self.write("rod.msgpack", rod_table(ks)),
            self.write("rodless.msgpack", rodless_table(KS)),
        )
        self.assertEqual((result["states"], result["bit_equal"]), (4, 3))
        self.assertEqual(result["field"], "k_inf")
        self.assertGreater(result["max_relative_difference"], 0.0)
        self.assertLess(result["max_relative_difference"], 1e-15)
        self.assertAlmostEqual(result["max_k_inf_pcm"], abs(ks[2] - 0.98) / 0.98 * 1e5)

    def test_a_node_table_without_the_rod_axis_is_compared_point_by_point(self):
        result = compare_rodless.compare(
            self.write("a.msgpack", rodless_table(KS)),
            self.write("b.msgpack", rodless_table(KS)),
        )
        self.assertFalse(result["rod_axis"])
        self.assertEqual((result["states"], result["bit_equal"]), (4, 4))

    def test_media_are_compared_state_by_state(self):
        result = compare_rodless.compare(
            self.write("a.msgpack", medium(1.0)), self.write("b.msgpack", medium(1.5))
        )
        self.assertEqual(result["kind"], "medium")
        self.assertEqual((result["states"], result["bit_equal"]), (4, 1))
        self.assertEqual(result["field"], "total")
        self.assertAlmostEqual(result["max_relative_difference"], 0.5 / 1.5)
        self.assertIsNone(result["max_k_inf_pcm"])

    def test_mismatched_shapes_are_errors(self):
        rodless = self.write("rodless.msgpack", rodless_table(KS))
        cases = {
            "first rod value": rod_table(KS, first="aic"),
            "axes": {
                **rod_table(KS),
                "axes": [
                    CONTINUOUS[0],
                    {**CONTINUOUS[1], "points": [550.0, 570.0]},
                    ROD,
                ],
            },
            "split": {**rod_table(KS), "split": 2},
            "point": {
                **rod_table(KS),
                "points": rod_table(KS)["points"][3:],
                "states": rod_table(KS)["states"][3:],
            },
        }
        for what, table in cases.items():
            with self.subTest(what), self.assertRaises(SystemExit):
                compare_rodless.compare(self.write("rod.msgpack", table), rodless)
        with self.assertRaises(SystemExit):
            compare_rodless.compare(self.write("m.msgpack", medium(1.0)), rodless)
        with self.assertRaises(SystemExit):
            compare_rodless.compare(
                self.write("a.msgpack", {**rodless_table(KS), "micro": {}}), rodless
            )

    def test_json_writes_nothing_to_stdout_and_report_names_each_output(self):
        rod = self.write("rod.msgpack", rod_table(KS))
        rodless = self.write("rodless.msgpack", rodless_table(KS))
        a, b = (
            self.write("a.msgpack", medium(1.0)),
            self.write("b.msgpack", medium(1.0)),
        )
        out = self.root / "compare.json"
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            compare_rodless.main(["--json", str(out), f"{rod}={rodless}", f"{a}={b}"])
        self.assertEqual(stdout.getvalue(), "")
        summary = json.loads(out.read_text())
        self.assertEqual(summary["format"], "xlc-compare-rodless")
        self.assertTrue(summary["all_bit_equal"])
        self.assertEqual([o["name"] for o in summary["outputs"]], [rod, a])
        self.assertEqual([i["path"] for i in summary["inputs"]], [rod, rodless, a, b])
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            compare_rodless.main(["--report", str(out)])
        text = stdout.getvalue()
        self.assertIn(f"{rod}: table with the rod axis, 4 of 4 states bit-equal", text)
        self.assertIn("all bit-equal: yes", text)

    def test_a_zero_of_the_other_sign_counts_as_a_difference_of_zero(self):
        rod, rodless = rod_table(KS), rodless_table(KS)
        rod["states"][0]["scatter"][1][0] = -0.0
        result = compare_rodless.compare(
            self.write("rod.msgpack", rod), self.write("rodless.msgpack", rodless)
        )
        self.assertEqual((result["states"], result["bit_equal"]), (4, 3))
        self.assertEqual(result["field"], "scatter")
        self.assertEqual(result["max_relative_difference"], 0.0)

    def test_a_bool_and_an_int_are_different_items(self):
        rod, rodless = rod_table(KS), rodless_table(KS)
        rod["states"][0]["pins"] = True
        rodless["states"][0]["pins"] = 1
        result = compare_rodless.compare(
            self.write("rod.msgpack", rod), self.write("rodless.msgpack", rodless)
        )
        self.assertEqual((result["bit_equal"], result["field"]), (3, "pins"))

    def test_non_finite_values_are_errors(self):
        rodless = self.write("rodless.msgpack", rodless_table(KS))
        for x in (float("nan"), float("inf")):
            rod = rod_table(KS)
            rod["states"][3]["total"][0] = x
            with self.subTest(x), self.assertRaises(SystemExit):
                compare_rodless.compare(self.write("rod.msgpack", rod), rodless)

    def test_states_not_matching_the_points_are_errors(self):
        rodless = self.write("rodless.msgpack", rodless_table(KS))
        short = rod_table(KS)
        short["states"] = short["states"][:-1]
        with self.assertRaises(SystemExit):
            compare_rodless.compare(self.write("rod.msgpack", short), rodless)
        with self.assertRaises(SystemExit):
            compare_rodless.compare(
                self.write("rod.msgpack", rod_table(KS)),
                self.write(
                    "r.msgpack", {**rodless_table(KS), "states": [state(1.0)] * 3}
                ),
            )

    def test_out_points_without_a_rodless_point_are_errors(self):
        rodless = rodless_table(KS)
        rodless["points"] = rodless["points"][1:]
        rodless["states"] = rodless["states"][1:]
        with self.assertRaises(SystemExit):
            compare_rodless.compare(
                self.write("rod.msgpack", rod_table(KS)),
                self.write("rodless.msgpack", rodless),
            )

    def test_a_pair_without_equals_is_an_error(self):
        with self.assertRaises(SystemExit):
            compare_rodless.main(["--json", str(self.root / "c.json"), "a.msgpack"])


if __name__ == "__main__":
    unittest.main()
