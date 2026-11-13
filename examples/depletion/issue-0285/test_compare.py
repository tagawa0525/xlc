"""compare.py のピン出力、係数の比、リングの捕獲、ログの読み方（#0285、ADR-0250）。

実行: cd examples/depletion/issue-0285 && python3 -m unittest test_compare
"""

import sys
import tempfile
import unittest
from pathlib import Path
from typing import ClassVar

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare import coefficients, log, power, ring_worst, worst_ratio


def state(form_fast: list, form_thermal: list) -> dict:
    """1 行のセルの 2 群の状態。φ2/φ1 = 0.02 / (0.08 + 0.02) = 0.2、κΣf は [1, 5]"""
    return {
        "absorption": [0.01, 0.08],
        "scatter": [[0.5, 0.02], [0.02, 1.0]],
        "kappa_fission": [1.0, 5.0],
        "pins": {"form": [[form_fast], [form_thermal]]},
    }


class Power(unittest.TestCase):
    def test_cells_weigh_the_group_forms_by_the_group_fission_power(self):
        # 群の出力は 1 × 1 と 5 × 0.2 で等しく、セルの出力は 2 群の形状関数の平均
        cells = power(state([1.2, 0.8], [0.6, 1.4]))
        self.assertAlmostEqual(cells[0], 0.9)
        self.assertAlmostEqual(cells[1], 1.1)


def grid(k: list) -> list[dict]:
    """[減速材][燃料] の k∞ の 9 状態"""
    return [{"k_inf": k[m][f]} for m in range(3) for f in range(3)]


class Coefficients(unittest.TestCase):
    K: ClassVar = [[1.02, 1.01, 1.00], [1.03, 1.02, 1.01], [1.01, 1.00, 0.99]]

    def test_reactivities_are_taken_from_the_nominal_points(self):
        c = coefficients(grid(self.K))
        self.assertAlmostEqual(c["fuel"]["m1f0"], (1 / 1.02 - 1 / 1.03) * 1e5)
        self.assertAlmostEqual(c["moderator"]["m2f1"], (1 / 1.02 - 1 / 1.00) * 1e5)
        self.assertEqual(len(c["fuel"]), 6)
        self.assertEqual(len(c["moderator"]), 6)

    def test_the_worst_ratio_is_the_largest_relative_change(self):
        every = coefficients(grid(self.K))
        changed = [row[:] for row in self.K]
        changed[0][0] = 1.03
        worst = worst_ratio(coefficients(grid(changed)), every)
        self.assertEqual(worst["fuel"]["at"], "m0f0")
        self.assertEqual(worst["moderator"]["at"], "m0f0")
        # 1/k の差が 1/1.01 − 1/1.02 から 1/1.01 − 1/1.03 へ、ほぼ 2 倍
        self.assertAlmostEqual(
            worst["fuel"]["ratio_minus_one"],
            (1 / 1.01 - 1 / 1.03) / (1 / 1.01 - 1 / 1.02) - 1.0,
        )


class Rings(unittest.TestCase):
    def test_the_worst_ring_capture_names_its_state_and_material(self):
        every = [{("fuel-0", "U-238"): 0.1, ("fuel-4", "U-238"): 0.2}] * 2
        ratio = [
            {("fuel-0", "U-238"): 0.1, ("fuel-4", "U-238"): 0.2},
            {("fuel-0", "U-238"): 0.099, ("fuel-4", "U-238"): 0.204},
        ]
        worst = ring_worst(ratio, every)["U-238"]
        self.assertEqual((worst["state"], worst["material"]), (1, "fuel-4"))
        self.assertAlmostEqual(worst["relative"], 0.02)

    def test_the_log_gives_spans_states_and_ring_rates(self):
        text = (
            "t  INFO shielding: xlc_lattice::shielding: close time.busy=45.0s time.idle=1µs\n"
            't DEBUG ring_rates: ring reaction rates per neutron produced material="fuel-0" '
            'nuclide="U-238" capture=0.05 fission=0.0\n'
            "t  INFO lattice: lattice state=0 iterations=3 tracks=1 seconds=50.0 total_seconds=51\n"
            "t  INFO ratio-shielding: xlc_lattice::shielding: close time.busy=500ms time.idle=1µs\n"
            "t  INFO lattice: lattice state=1 iterations=3 tracks=1 seconds=6.0 total_seconds=57\n"
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "a.log"
            path.write_text(text)
            parsed = log(path)
        self.assertEqual(parsed["spans"]["shielding"], [45.0])
        self.assertEqual(parsed["spans"]["ratio-shielding"], [0.5])
        self.assertEqual(parsed["seconds"], [50.0, 6.0])
        self.assertEqual(parsed["rings"], [{("fuel-0", "U-238"): 0.05}, {}])


if __name__ == "__main__":
    unittest.main()
