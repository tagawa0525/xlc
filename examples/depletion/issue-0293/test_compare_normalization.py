"""compare_normalization.py（#0293 の項目 3 の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0293 -p 'test_compare_normalization.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_normalization


def result(
    k,
    fission_power,
    u235,
    power=100.0,
    burnup=(0.0, 10.0, 20.0),
    normalization="energy-deposition",
):
    """openmc_depletion.py の extract の形（リング 2 つ、体積 1 と 3。U235 のほかの核種は 1）。"""
    others = {
        n: [1.0] * len(burnup) for n in ["Pu239", "Pu240", "Pu241", "Xe135", "Sm149"]
    }
    return {
        "burnup": list(burnup),
        "k": list(k),
        "k-std": [1.0e-4] * len(burnup),
        "rings": [
            {
                "name": "fuel-0",
                "volume": 1.0,
                "densities": {"U235": list(u235), **others},
            },
            {
                "name": "fuel-1",
                "volume": 3.0,
                "densities": {"U235": list(u235), **others},
            },
        ],
        "transportable": ["U235", "U238"],
        "power": power,
        "fission-power": list(fission_power),
        "fission-rates": {"U235": [1.0] * len(burnup), "Pu239": [0.0, 0.5, 1.0]},
        "fission-q-values": {"U235": 2.0e8, "Pu239": 2.0e8},
        "normalization": normalization,
    }


FISSION_Q = result(
    [1.2, 1.1, 1.0], [100.0, 100.0, 100.0], [1.0, 0.8, 0.6], normalization="fission-q"
)
DEPOSITION = result([1.2, 1.11, 1.01], [100 / 1.02] * 3, [1.0, 0.81, 0.62])


class Compare(unittest.TestCase):
    def test_the_extra_energy_is_the_set_power_over_the_fission_power(self):
        deposition = result(
            [1.2, 1.11, 1.01], [100 / 1.02, 100 / 1.04, 100 / 1.06], [1.0, 0.81, 0.62]
        )
        rows = compare_normalization.compare(FISSION_Q, deposition)
        self.assertAlmostEqual(rows[0]["delta"], 0.02)
        self.assertAlmostEqual(rows[2]["delta"], 0.06)

    def test_the_prediction_shifts_the_burnup_by_the_mean_extra_energy(self):
        deposition = result(
            [1.2, 1.11, 1.01], [100 / 1.02, 100 / 1.04, 100 / 1.06], [1.0, 0.81, 0.62]
        )
        rows = compare_normalization.compare(FISSION_Q, deposition)
        # δ̄ は 0〜20 MWd/kg の台形の平均: ((0.02 + 0.04) / 2 · 10 + (0.04 + 0.06) / 2 · 10) / 20 = 0.04
        self.assertAlmostEqual(rows[2]["mean_delta"], 0.04)
        self.assertAlmostEqual(rows[2]["burnup_shift"], 0.8)
        # −(dk/dB) B δ̄ / k: dk/dB は fission-q の 1 つ前の点との差分 −0.01
        self.assertAlmostEqual(rows[2]["predicted_pcm"], 0.01 * 20 * 0.04 / 1.0 * 1e5)
        self.assertIsNone(rows[0]["predicted_pcm"])

    def test_the_measured_difference_is_relative_to_the_fission_q_k(self):
        deposition = result([1.2, 1.11, 1.01], [100 / 1.02] * 3, [1.0, 0.81, 0.62])
        rows = compare_normalization.compare(FISSION_Q, deposition)
        self.assertAlmostEqual(rows[1]["measured_pcm"], (1.11 - 1.1) / 1.1 * 1e5)
        self.assertAlmostEqual(
            rows[1]["measured_std_pcm"], (2 * 1.0e-8) ** 0.5 / 1.1 * 1e5
        )

    def test_the_densities_are_volume_averaged_relative_differences(self):
        deposition = result([1.2, 1.11, 1.01], [100 / 1.02] * 3, [1.0, 0.88, 0.66])
        rows = compare_normalization.compare(FISSION_Q, deposition)
        self.assertAlmostEqual(rows[1]["densities"]["U235"], 0.88 / 0.8 - 1)
        self.assertAlmostEqual(rows[1]["densities"]["Pu239"], 0.0)

    def test_a_fission_q_run_off_its_set_power_is_an_error(self):
        broken = result(
            [1.2, 1.1, 1.0],
            [100.0, 99.0, 100.0],
            [1.0, 0.8, 0.6],
            normalization="fission-q",
        )
        with self.assertRaises(SystemExit):
            compare_normalization.compare(broken, DEPOSITION)

    def test_runs_not_normalized_as_named_are_an_error(self):
        with self.assertRaises(SystemExit):
            compare_normalization.compare(FISSION_Q, FISSION_Q)
        with self.assertRaises(SystemExit):
            compare_normalization.compare(DEPOSITION, DEPOSITION)

    def test_runs_on_different_burnup_points_are_an_error(self):
        other = result(
            [1.2, 1.1, 1.0], [100.0] * 3, [1.0, 0.8, 0.6], burnup=(0.0, 10.0, 30.0)
        )
        with self.assertRaises(SystemExit):
            compare_normalization.compare(FISSION_Q, other)

    def test_runs_transporting_different_nuclides_are_an_error(self):
        other = result([1.2, 1.1, 1.0], [100.0] * 3, [1.0, 0.8, 0.6])
        other["transportable"] = ["U235"]
        with self.assertRaises(SystemExit):
            compare_normalization.compare(FISSION_Q, other)


class Kappa(unittest.TestCase):
    """仕様書の核種ごとの定数 Kappa / Kappa0 を、fission-q の履歴の核分裂の率と、extract が書いたチェーンの Q で重み付けた割合。"""

    def test_the_kappa_ratio_is_weighted_by_the_chain_q_fission_power(self):
        kappa = {
            "U235": {"kappa": 1.04, "kappa0": 1.0},
            "Pu239": {"kappa": 1.08, "kappa0": 1.0},
        }
        rows = compare_normalization.kappa_delta(FISSION_Q, kappa)
        self.assertAlmostEqual(rows[0]["delta"], 0.04)
        # 点 2 は U235 と Pu239 が 1 : 1
        self.assertAlmostEqual(rows[2]["delta"], 0.06)
        self.assertAlmostEqual(rows[2]["uncovered"], 0.0)

    def test_fission_without_a_kappa_is_counted_as_uncovered(self):
        kappa = {"U235": {"kappa": 1.04, "kappa0": 1.0}}
        rows = compare_normalization.kappa_delta(FISSION_Q, kappa)
        self.assertAlmostEqual(rows[2]["delta"], 0.04)
        self.assertAlmostEqual(rows[2]["uncovered"], 0.5)


if __name__ == "__main__":
    unittest.main()
