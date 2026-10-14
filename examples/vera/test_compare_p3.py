"""compare_p3.py の統計誤差のツール（OpenMC の *-std との比）。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_p3.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_p3


class OffsetStd(unittest.TestCase):
    def test_propagates_independent_layer_errors(self):
        """AO = (U − L) / (U + L)。下半分の層だけに誤差 0.1 があるとき ∂AO/∂p = −2U/(U+L)² = −0.5。"""
        bounds, axial = [0.0, 1.0, 2.0], [1.0, 1.0]
        self.assertAlmostEqual(
            compare_p3.offset_std(bounds, axial, [0.1, 0.0]), 0.05, places=12
        )
        self.assertAlmostEqual(
            compare_p3.offset_std(bounds, axial, [0.1, 0.1]),
            0.05 * math.sqrt(2),
            places=12,
        )


class SigmaMultiple(unittest.TestCase):
    def test_returns_the_largest_ratio_and_its_index(self):
        multiple, at = compare_p3.sigma_multiple([1.0, -6.0, 2.0], [1.0, 2.0, 1.0])
        self.assertAlmostEqual(multiple, 3.0, places=12)
        self.assertEqual(at, 1)


class PinSigma(unittest.TestCase):
    def test_median_relative_std_and_largest_multiple_over_fuel_cells(self):
        ours = [[1.0, 0.0], [2.2, 2.0]]
        ref = [[1.0, 0.0], [2.0, 2.0]]
        std = [[0.1, 0.0], [0.1, 0.1]]
        median, multiple, at = compare_p3.pin_sigma(ours, ref, std)
        self.assertAlmostEqual(median, 0.05, places=12)
        self.assertAlmostEqual(multiple, 2.0, places=12)
        self.assertEqual(at, (1, 0))


if __name__ == "__main__":
    unittest.main()
