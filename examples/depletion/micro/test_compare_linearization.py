"""compare_linearization.py（ミクロ燃焼の線形の補正が、追う核種の組成の差の輸送の k∞ をどれだけ再現するか）。

実行: python3 -m unittest discover -s examples/depletion/micro -p 'test_compare_linearization.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_linearization as cl


def state(total, scatter, nu_fission, k_inf):
    return {
        "k_inf": k_inf,
        "total": total,
        "scatter": scatter,
        "nu_fission": nu_fission,
        "chi": [1.0, 0.0],
    }


def two_group_k(total, scatter, nu_fission):
    """χ = (1, 0)、上方散乱なしの 2 群の k∞"""
    a1 = total[0] - scatter[0][0] - scatter[0][1]
    a2 = total[1] - scatter[1][0] - scatter[1][1]
    s12 = scatter[0][1]
    return (nu_fission[0] + nu_fission[1] * s12 / a2) / (a1 + s12)


TOTAL = [0.5, 1.2]
SCATTER = [[0.45, 0.02], [0.0, 1.1]]
NU_FISSION = [0.006, 0.14]
# 追う核種 1 つの微視的断面積 (barn)
SIGMA = {
    "total": [10.0, 400.0],
    "scatter": [[9.0, 0.0], [0.0, 10.0]],
    "nu_fission": [1.0, 700.0],
    "fission": [0.4, 290.0],
    "reactions": [],
    "thermal_fission": 0.9,
}


def table(densities, shifted_k=None):
    """燃焼度 2 点、状態は点ごとに 1 つ。点 1 は点 0 と同じ少数群定数で、k∞ は `shifted_k`（なければ点 0 の値）"""
    k0 = two_group_k(TOTAL, SCATTER, NU_FISSION)
    states = [
        state(TOTAL, SCATTER, NU_FISSION, k0),
        state(TOTAL, SCATTER, NU_FISSION, k0 if shifted_k is None else shifted_k),
    ]
    out = {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "points": [[0], [1]],
        "states": states,
    }
    if densities is not None:
        out["micro"] = {
            "nuclides": ["U-235"],
            "densities": densities,
            "states": [[SIGMA], [SIGMA]],
        }
    return out


class KInf(unittest.TestCase):
    def test_two_group_k_is_the_infinite_medium_eigenvalue(self):
        got = cl.k_inf(state(TOTAL, SCATTER, NU_FISSION, 0.0))
        self.assertAlmostEqual(got, two_group_k(TOTAL, SCATTER, NU_FISSION), places=14)

    def test_the_linear_correction_adds_sigma_times_the_density_change(self):
        corrected = cl.corrected(
            state(TOTAL, SCATTER, NU_FISSION, 0.0), [SIGMA], [2.0e-4]
        )
        self.assertAlmostEqual(corrected["total"][1], 1.2 + 400.0 * 2.0e-4, places=14)
        self.assertAlmostEqual(
            corrected["nu_fission"][1], 0.14 + 700.0 * 2.0e-4, places=14
        )
        self.assertAlmostEqual(
            corrected["scatter"][1][1], 1.1 + 10.0 * 2.0e-4, places=14
        )
        self.assertEqual(corrected["chi"], [1.0, 0.0])


class Compare(unittest.TestCase):
    def test_rows_compare_the_linear_prediction_with_the_shifted_transport(self):
        densities = [[1.0e-4], [3.0e-4]]
        delta = 2.0e-4
        total = [t + s * delta for t, s in zip(TOTAL, SIGMA["total"])]
        scatter = [
            [x + s * delta for x, s in zip(row, srow)]
            for row, srow in zip(SCATTER, SIGMA["scatter"])
        ]
        nu_fission = [n + s * delta for n, s in zip(NU_FISSION, SIGMA["nu_fission"])]
        predicted = two_group_k(total, scatter, nu_fission)
        base = two_group_k(TOTAL, SCATTER, NU_FISSION)
        # 輸送の値は予測より 3 pcm 大きい反応度（1 − 1/k）のもの。誤差は ρ(予測) − ρ(輸送) で −3 pcm
        transport = 1.0 / (1.0 / predicted - 3.0e-5)
        summary = cl.compare(table(densities), table(None, transport))
        self.assertEqual(summary["format"], "xlc-compare-micro-linearization")
        self.assertEqual(summary["nuclides"], 1)
        (row,) = summary["rows"]
        self.assertEqual((row["from"], row["to"]), (0.0, 10.0))
        self.assertAlmostEqual(
            row["tracked_pcm"], (1 / base - 1 / transport) * 1e5, places=9
        )
        self.assertAlmostEqual(
            row["linear_pcm"], (1 / base - 1 / predicted) * 1e5, places=9
        )
        self.assertAlmostEqual(row["error_pcm"], -3.0, places=9)
        self.assertAlmostEqual(
            row["error_pcm"], row["linear_pcm"] - row["tracked_pcm"], places=9
        )
        self.assertAlmostEqual(row["density_change_percent"]["U-235"], 200.0, places=9)
        self.assertLess(summary["self_check_pcm"], 1e-9)

    def test_the_error_is_the_linear_minus_the_transport_effect_from_their_own_bases(
        self,
    ):
        # テーブルの k∞ が 2 群の k と少し違っても（収束の誤差）、誤差は 2 つの効きの差そのもの
        densities = [[1.0e-4], [3.0e-4]]
        nominal = table(densities)
        nominal["states"][0]["k_inf"] *= 1.0 + 1.0e-7
        shifted = table(None, 1.1)
        shifted["states"][0]["k_inf"] = nominal["states"][0]["k_inf"]
        row = cl.compare(nominal, shifted)["rows"][0]
        self.assertAlmostEqual(
            row["error_pcm"], row["linear_pcm"] - row["tracked_pcm"], places=10
        )

    def test_branch_axes_with_one_point_are_allowed(self):
        nominal, shifted = table([[1.0e-4], [3.0e-4]]), table(None)
        for t in (nominal, shifted):
            t["axes"].append({"name": "fuel", "points": [565.0]})
            t["points"] = [[0, 0], [1, 0]]
        self.assertEqual(len(cl.compare(nominal, shifted)["rows"]), 1)
        nominal["axes"][1]["points"] = [565.0, 600.0]
        with self.assertRaises(SystemExit):
            cl.compare(nominal, shifted)

    def test_tables_without_micro_data_or_with_other_burnups_are_errors(self):
        with self.assertRaises(SystemExit):
            cl.compare(table(None), table(None))
        other = table(None)
        other["axes"][0]["points"] = [0.0, 20.0]
        with self.assertRaises(SystemExit):
            cl.compare(table([[1.0e-4], [3.0e-4]]), other)

    def test_the_fresh_points_must_agree(self):
        # 燃焼度 0 の組成は 2 つの履歴で同じなので、k∞ も同じ
        shifted = table(None)
        shifted["states"][0]["k_inf"] *= 1.001
        with self.assertRaises(SystemExit):
            cl.compare(table([[1.0e-4], [3.0e-4]]), shifted)


if __name__ == "__main__":
    unittest.main()
