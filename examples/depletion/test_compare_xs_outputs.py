"""compare_xs_outputs.py（2 つの `xs nuclide` の出力を値で比べる。#0325）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_xs_outputs.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_xs_outputs as cxo


def output():
    """2 群、P0 と P1、σ0 が 1 つの核種 1 つの出力"""
    shielded = {name: [1.0, 2.0] for name in cxo.SHIELDED}
    temperature = {
        "temperature": 600.0,
        "chi": [0.75, 0.25],
        "scatter": [[[3.0, 1.0], [0.0, 2.0]], [[0.3, 0.1], [0.0, 0.2]]],
        "shielded": [shielded],
    }
    return {"nuclides": [{"name": "U-235", "temperatures": [temperature]}]}


class CompareTest(unittest.TestCase):
    def test_identical_outputs_are_the_same(self):
        self.assertEqual(cxo.compare(output(), output()), ["U-235 600 K same"])

    def test_scatter_difference_is_relative_to_the_p0_row_sum(self):
        b = output()
        b["nuclides"][0]["temperatures"][0]["scatter"][1][0][1] = 0.14
        self.assertEqual(
            cxo.compare(output(), b),
            ["U-235 600 K scatter P1: max |diff| / P0 row sum 1.000e-02 in 0 -> 1"],
        )

    def test_chi_and_shielded_differences_name_their_group(self):
        b = output()
        t = b["nuclides"][0]["temperatures"][0]
        t["chi"] = [0.7, 0.3]
        t["shielded"][0]["capture"][1] = 2.2
        lines = cxo.compare(output(), b)
        self.assertEqual(
            lines[0], "U-235 600 K chi: max |diff| 5.000e-02 in group 0, sum 1.000e-01"
        )
        self.assertEqual(
            lines[1],
            "U-235 600 K shielded capture: max rel diff 9.091e-02 at sigma0 0, group 1",
        )

    def test_different_group_counts_stop(self):
        # 散乱行列の行の数と自己遮蔽した量の長さが違う出力は、短い方に切り詰めて比べずに止める
        rows = output()
        rows["nuclides"][0]["temperatures"][0]["scatter"][0].pop()
        shielded = output()
        shielded["nuclides"][0]["temperatures"][0]["shielded"][0]["capture"].pop()
        for b in (rows, shielded):
            with self.assertRaises(ValueError):
                cxo.compare(output(), b)

    def test_different_group_bounds_or_sigma0_stop(self):
        # 群の数が同じでも、群境界か σ0 が違う出力は比べない
        for key, value in (("group_bounds", [2.0e7, 1.0, 1.0e-5]), ("sigma0", [1.0e5])):
            a, b = output(), output()
            a[key], b[key] = [2.0e7, 1.0e3, 1.0e-5], value
            if key == "sigma0":
                a[key] = [1.0e10]
            with self.assertRaises(ValueError):
                cxo.compare(a, b)

    def test_different_nuclides_stop(self):
        b = copy.deepcopy(output())
        b["nuclides"][0]["name"] = "U-238"
        with self.assertRaises(ValueError):
            cxo.compare(output(), b)


if __name__ == "__main__":
    unittest.main()
