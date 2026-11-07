"""kinetics_table.py の要約（燃焼度の点ごとの公称の状態の動特性と冷却計算のデータ。ADR-0238）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_kinetics_table.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import kinetics_table as kt


def state(scale: float) -> dict:
    """2 群の状態。上方散乱なしで、k∞ = (νΣf1 + νΣf2 Σs12 / Σr2) / Σr1"""
    return {
        "k_inf": (0.005 * scale + 0.1 * scale * 0.02 / 0.1) / 0.03,
        "total": [0.5, 1.0],
        "scatter": [[0.47, 0.02], [0.0, 0.9]],
        "chi": [1.0, 0.0],
        "nu_fission": [0.005 * scale, 0.1 * scale],
    }


def kinetics(scale: float) -> dict:
    """2 族。族 0 は崩壊定数 0.0125 /s、族 1 は群で違う核種の混ざり方"""
    return {
        "delayed_nu_fission": [[1.0e-5 * scale, 2.0e-4 * scale], [3.0e-5, 6.0e-4]],
        "delayed_precursors": [[8.0e-4 * scale, 1.6e-2 * scale], [1.0e-5, 4.0e-4]],
        "delayed_chi": [[1.0, 0.0], [1.0, 0.0]],
        "inverse_velocity": [5.0e-8, 2.0e-6],
    }


def table() -> dict:
    """燃焼度 2 点 × ホウ素 2 点（公称 1300）。公称でない状態は燃焼度 0 だけ"""
    points = [[0, 0], [0, 1], [1, 1]]
    scales = [1.0, 1.1, 0.9]
    cooling = lambda s: [
        {
            "fission": 40.0 * s,
            "thermal_fission": 0.9,
            "reactions": [{"mt": 102, "channels": [8.0 * s, 1.0 * s]}],
        },
        {
            "fission": 0.0,
            "thermal_fission": None,
            "reactions": [{"mt": 16, "channels": [0.01]}],
        },
    ]
    return {
        "axes": [
            {"name": "burnup", "points": [0.0, 10.0]},
            {"name": "boron", "points": [1000.0, 1300.0], "nominal": 1300.0},
        ],
        "points": points,
        "states": [state(s) for s in scales],
        "kinetics": {"states": [kinetics(s) for s in scales]},
        "cooling": {
            "nuclides": ["U-235", "Cs-133"],
            "points": [cooling(1.0), cooling(0.5)],
        },
    }


class Summary(unittest.TestCase):
    def test_takes_the_nominal_state_at_each_burnup(self):
        self.assertEqual(kt.nominal_points(table()), [1, 2])

    def test_a_missing_nominal_state_is_an_error(self):
        t = table()
        t["points"] = t["points"][:2]
        with self.assertRaises(ValueError):
            kt.nominal_points(t)

    def test_the_fundamental_mode_gives_the_lattice_k(self):
        data = kt.summary(table(), ["U-235", "Cs-133"])
        for row in data["points"]:
            self.assertAlmostEqual(row["k_minus_lattice_pcm"], 0.0, places=9)

    def test_beta_and_decay_constants_weight_the_families_by_the_flux(self):
        # φ = A⁻¹ χ: φ1 = 1/0.03、φ2 = 0.02 φ1 / 0.1。族 0 は群によらず ν_dΣf / (ν_dΣf/λ) = 1/80
        data = kt.summary(table(), ["U-235"])
        row = data["points"][0]
        phi = [1.0, 0.2]
        scale = 1.1
        production = 0.005 * scale * phi[0] + 0.1 * scale * phi[1]
        beta0 = (1.0e-5 * phi[0] + 2.0e-4 * phi[1]) * scale / production
        self.assertAlmostEqual(row["beta_pcm"][0], beta0 * 1e5, places=9)
        self.assertAlmostEqual(row["decay_constants"][0], 1.0 / 80.0, places=12)
        lam1 = (3.0e-5 * phi[0] + 6.0e-4 * phi[1]) / (1.0e-5 * phi[0] + 4.0e-4 * phi[1])
        self.assertAlmostEqual(row["decay_constants"][1], lam1, places=12)
        self.assertAlmostEqual(
            row["inverse_velocity"],
            (5.0e-8 * phi[0] + 2.0e-6 * phi[1]) / 1.2,
            places=18,
        )

    def test_cooling_keeps_the_picked_nuclides_and_sums_the_capture_channels(self):
        data = kt.summary(table(), ["U-235"])
        self.assertEqual(data["cooling_nuclides"], 2)
        second = data["points"][1]["cooling"]
        self.assertEqual(list(second), ["U-235"])
        self.assertAlmostEqual(second["U-235"]["capture"], 4.5)
        self.assertAlmostEqual(second["U-235"]["fission"], 20.0)

    def test_a_picked_nuclide_without_cooling_data_is_an_error(self):
        with self.assertRaises(ValueError):
            kt.summary(table(), ["U-235", "Pu-239"])

    def test_a_repeated_picked_nuclide_is_an_error(self):
        with self.assertRaises(ValueError):
            kt.summary(table(), ["U-235", "U-235"])

    def test_a_family_without_precursors_under_the_flux_is_an_error(self):
        t = table()
        for k in t["kinetics"]["states"]:
            k["delayed_precursors"][1] = [0.0, 0.0]
        with self.assertRaises(ValueError):
            kt.summary(t, ["U-235"])


if __name__ == "__main__":
    unittest.main()
