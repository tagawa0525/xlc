"""compare_rods_switch.py（制御棒を途中まで入れて抜く燃焼の、補正あり・なしの core と格子計算の参照の比較）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_rods_switch.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_rods_switch as c
from compare_residual import added, k_inf, rho


def few(s):
    """2 群の少数群定数（`s` で値を変える）"""
    x = 1.0 + s / 7.0
    return {
        "total": [0.55 * x, 1.3 * x],
        "absorption": [0.01 * x, 0.1 * x],
        "nu_fission": [0.007 * x, 0.14 * x],
        "fission": [0.0028 * x, 0.057 * x],
        "chi": [1.0, 0.0],
        "diffusion": [1.4 * x, 0.4 * x],
        "scatter": [[0.52 * x, 0.017 * x], [0.0015 * x, 1.2 * x]],
        "adf": [[1.0, 1.01]],
        "k_inf": 1.1 + s / 100.0,
    }


def residual(s):
    """形だけ合わせた残り（少数群定数の項目ごとの差）"""
    f = few(0.0)
    return {
        key: [[s * 1e-3 * v for v in row] for row in value]
        if isinstance(value[0], list)
        else [s * 1e-3 * v for v in value]
        for key, value in f.items()
        if key != "k_inf"
    }


def sigma(s):
    return {
        "total": [11.0 + s, 130.0 + s],
        "nu_fission": [1.1, 170.0 + s],
        "fission": [0.4, 70.0],
        "scatter": [[9.0, 0.1], [0.01, 8.0]],
    }


def table(l_density=1.0, l_rods=2.0):
    """燃焼度 0 と 20、制御棒の軸（out、aic）、追う核種 2 つの公称のテーブルと、2 つの変数の制御棒の値 out の補正"""
    axes = [
        {"name": "burnup", "points": [0.0, 20.0]},
        {"name": "boron", "points": [1300.0], "nominal": 1300.0},
        {"name": "moderator", "points": [565.0, 600.0], "nominal": 565.0},
        {"name": "fuel", "points": [565.0], "nominal": 565.0},
        {"name": "rod", "values": ["out", "aic"]},
    ]
    points = [[b, 0, m, 0, r] for b in range(2) for m in range(2) for r in range(2)]
    nominal = [[7.0e-4, 0.0], [5.0e-4, 2.0e-5]]

    def history(variable, densities, l):
        return {
            "variable": variable,
            "rods": [
                {
                    "rod": "out",
                    "burnup": [0.0, 20.0],
                    "untracked": [residual(0.0), residual(0.5)],
                    "linearization": [residual(0.0), residual(l)],
                    "densities": [nominal[0], densities],
                }
            ],
        }

    return {
        "axes": axes,
        "points": points,
        "states": [few(float(i)) for i in range(len(points))],
        "micro": {
            "nuclides": ["U-235", "Pu-240"],
            "densities": nominal,
            "states": [[sigma(float(i)), sigma(2.0 * i)] for i in range(len(points))],
        },
        "histories": [
            history("moderator-density", [5.2e-4, 2.6e-5], l_density),
            history("rods", [5.1e-4, 2.9e-5], l_rods),
        ],
    }


class StateIndex(unittest.TestCase):
    def test_the_state_is_the_nominal_state_of_the_rod_value(self):
        t = table()
        self.assertEqual(t["points"][c.state_index(t, 20.0, "aic")], [1, 0, 0, 0, 1])
        self.assertEqual(t["points"][c.state_index(t, 0.0, "out")], [0, 0, 0, 0, 0])


class Homogenized(unittest.TestCase):
    def history(self, area):
        return {
            "nuclides": ["U-235", "U-238", "Pu-240"],
            "regions": [{"volume": 2.0}, {"volume": 3.0}],
            "points": [
                {
                    "burnup": 0.0,
                    "k_inf": 1.2,
                    "densities": [[1e-3, 2e-2, 0.0], [2e-3, 2e-2, 0.0]],
                },
                {
                    "burnup": 20.0,
                    "k_inf": 1.0,
                    "densities": [[5e-4, 2e-2, 1e-5], [1e-3, 2e-2, 3e-5]],
                },
            ],
        }

    def test_densities_are_volume_averages_over_the_lattice_area(self):
        # 格子の面積 10 なら、N̄ = Σ V N / 10
        area = 10.0
        nominal = [(2e-3 + 6e-3) / area, (4e-2 + 6e-2) / area, 0.0]
        nuclides = ["U-235", "U-238", "Pu-240"]
        rows = c.homogenized(self.history(area), nuclides, nominal)
        self.assertAlmostEqual(
            rows[1]["densities"][2], (2e-5 + 9e-5) / area, delta=1e-18
        )
        self.assertEqual(rows[1]["k_inf"], 1.0)

    def test_another_initial_fuel_is_an_error(self):
        nuclides = ["U-235", "U-238", "Pu-240"]
        with self.assertRaises(SystemExit):
            c.homogenized(self.history(10.0), nuclides, [9e-4, 1e-2, 0.0])


class Nodes(unittest.TestCase):
    def test_cores_without_the_rod_history_are_errors(self):
        core = {
            "format": "xlc-core-history",
            "assemblies": ["A01"],
            "axial_mesh": [0.0, 21.5],
            "nuclides": ["U-235"],
            "history_variables": ["moderator-density"],
            "points": [],
        }
        with self.assertRaises(SystemExit):
            c.nodes(core, "core", ["U-235"])


class Split(unittest.TestCase):
    def test_a_node_at_the_expected_composition_takes_the_weights(self):
        # ノードの組成が見込み Σ w_j ΔN̄^j なら、ADR-0228 の a_j は w_j（Σ_j a_j P_j = P）
        (row,) = c.split_rows(table())
        for share in row["shares"]:
            self.assertAlmostEqual(share, c.SPLIT_WEIGHT, delta=1e-12)
        # 変数ごとに P / P_j で重ねると、どちらの変数にも w より大きい重みが掛かる
        for a in row["independent"]:
            self.assertGreater(a, c.SPLIT_WEIGHT)

    def test_without_linearizations_the_two_forms_are_the_same(self):
        (row,) = c.split_rows(table(0.0, 0.0))
        self.assertEqual(row["overcount_pcm"], 0.0)
        self.assertEqual(row["density_l_pcm"], 0.0)

    def test_the_overcount_adds_the_extra_weights_on_the_linearizations(self):
        t = table()
        (row,) = c.split_rows(t)
        # 重ねた形は、見込みの形より (P / P_j − w_j) L_j を余分に足す
        (none,) = c.split_rows(table(0.0, 0.0))
        self.assertNotEqual(row["overcount_pcm"], none["overcount_pcm"])
        w = c.SPLIT_WEIGHT
        extra_m, extra_r = (a - w for a in row["independent"])
        p = c.state_index(t, 20.0, "out")
        base = t["states"][p]
        lm = t["histories"][0]["rods"][0]["linearization"][1]
        lr = t["histories"][1]["rods"][0]["linearization"][1]
        # 線形の項の大きさで同じ向き（L の効きの符号が同じ）
        alone = rho(
            k_inf(added(added(base, c.scaled(lm, extra_m)), c.scaled(lr, extra_r)))
        ) - rho(k_inf(base))
        self.assertEqual(alone > 0, row["overcount_pcm"] > 0)


if __name__ == "__main__":
    unittest.main()
