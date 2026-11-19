"""compare_shielding.py（足す核種の自己遮蔽した 1 群の断面積と無限希釈の差）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_shielding.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_shielding


def table(sigmas):
    """燃焼度 0 と 10 の点の、核種ごとの (n,γ) だけの冷却計算の断面積の表"""
    names = list(sigmas)
    return {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "points": [[0], [1]],
        "cooling": {
            "nuclides": names,
            "points": [
                [
                    {
                        "fission": 0.0,
                        "reactions": [{"mt": 102, "channels": [sigmas[n][k]]}],
                    }
                    for n in names
                ]
                for k in range(2)
            ],
        },
    }


def absorption(shares):
    return {
        "points": [
            {
                "burnup": 10.0,
                "k_inf": 0.5,
                "untransported_all": {
                    n: {"fraction": f, "captured_ratio": r}
                    for n, (f, r) in shares.items()
                },
            }
        ]
    }


class Compare(unittest.TestCase):
    def test_the_shielded_share_scales_the_estimate_by_the_cross_section_ratio(self):
        infinite = table({"A": [1.0, 100.0], "B": [1.0, 10.0]})
        shielded = table({"A": [1.0, 80.0], "B": [1.0, 10.0]})
        shares = absorption({"A": (1e-3, 0.5), "B": (1e-4, 1.0)})
        rows = compare_shielding.compare(infinite, shielded, shares, ["A", "B"])
        a = next(r for r in rows if r["nuclide"] == "A")
        self.assertAlmostEqual(a["points"][0]["ratio"], 0.8)
        # 1e5 × 1e-3 × 0.5 × (1 − 0.8) / 0.5 = 20 pcm 吸収が減る（反応度は +20）
        self.assertAlmostEqual(a["points"][0]["shielding_pcm"], 20.0)
        self.assertAlmostEqual(a["max_pcm"], 20.0)
        self.assertEqual([r["nuclide"] for r in rows], ["A", "B"])

    def test_nuclides_missing_from_either_table_are_errors(self):
        t = table({"A": [1.0, 1.0]})
        with self.assertRaises(SystemExit):
            compare_shielding.compare(t, t, absorption({"Z": (1e-3, 1.0)}), ["Z"])


if __name__ == "__main__":
    unittest.main()
