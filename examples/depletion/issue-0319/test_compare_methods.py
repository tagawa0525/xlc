"""compare_methods.py（非分離共鳴の負のバックグラウンドの読み方ごとの、遮蔽した群定数の差）。

実行: python3 -m unittest discover -s examples/depletion/issue-0319 -p 'test_compare_methods.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_methods as cm


def library(capture, total, temperature=565.0):
    """3 群（10〜1e3〜1e5〜1e7 eV）、σ0 ごとの捕獲と全断面積だけを持つ 1 核種の xs のファイル"""
    return {
        "group_bounds": [1.0e7, 1.0e5, 1.0e3, 10.0],
        "sigma0": [1.0e10, 1.0],
        "nuclides": [
            {
                "name": "X",
                "temperatures": [
                    {
                        "temperature": temperature,
                        "shielded": [
                            {"capture": c, "total": t, "fission": [0.0] * 3}
                            for c, t in zip(capture, total)
                        ],
                    }
                ],
            }
        ],
    }


REFERENCE = library([[1.0, 2.0, 4.0], [1.0, 1.0, 2.0]], [[10.0] * 3, [8.0] * 3])


class DifferencesTest(unittest.TestCase):
    def test_takes_the_largest_relative_difference_of_the_groups_in_the_range(self):
        other = library([[1.0, 2.0, 4.0], [9.0, 1.1, 1.0]], [[10.0] * 3, [8.0] * 3])
        # 範囲 [2e3, 5e4] eV にかかるのは 2 群目（1e3〜1e5 eV）だけ
        (t,) = cm.differences(REFERENCE, other, 2.0e3, 5.0e4)
        self.assertEqual(t["temperature"], 565.0)
        capture = t["sigma0"][1]["capture"]
        self.assertEqual(capture["group"], 2)
        self.assertAlmostEqual(capture["relative"], 0.1)
        self.assertEqual((capture["reference"], capture["value"]), (1.0, 1.1))
        self.assertEqual(t["sigma0"][0]["capture"]["relative"], 0.0)

    def test_keeps_the_sign_of_a_negative_reference(self):
        other = library([[1.0, 2.0, 4.0], [1.0, 1.0, 2.0]], [[10.0] * 3, [8.0] * 3])
        negative = library([[1.0, 2.0, 4.0], [1.0, -0.5, 2.0]], [[10.0] * 3, [8.0] * 3])
        (t,) = cm.differences(negative, other, 2.0e3, 5.0e4)
        capture = t["sigma0"][1]["capture"]
        self.assertEqual((capture["reference"], capture["value"]), (-0.5, 1.0))
        self.assertAlmostEqual(capture["relative"], -3.0)

    def test_other_temperatures_or_groups_are_an_error(self):
        with self.assertRaises(SystemExit):
            cm.differences(
                REFERENCE, library([[1.0] * 3] * 2, [[1.0] * 3] * 2, 600.0), 2e3, 5e4
            )


if __name__ == "__main__":
    unittest.main()
