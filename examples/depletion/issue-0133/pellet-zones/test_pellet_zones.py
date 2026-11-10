"""pellet_zones.py（#0275 のペレットを径方向の単位に分けた燃焼の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133/pellet-zones -p 'test_pellet_zones.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pellet_zones as pz


class RangeMean(unittest.TestCase):
    def test_whole_units_inside_the_range_average_by_volume(self):
        self.assertAlmostEqual(pz.range_mean([1, 2, 3, 4], [1, 1, 1, 1], 0.0, 0.5), 1.5)
        self.assertAlmostEqual(pz.range_mean([1, 2, 3, 4], [1, 1, 1, 1], 0.5, 1.0), 3.5)

    def test_a_unit_cut_by_the_range_counts_its_share_of_the_volume(self):
        # 5 単位の 3 つ目（体積の 0.4〜0.6）は、半分ずつ内と外に入る
        values = [1, 2, 3, 4, 5]
        self.assertAlmostEqual(pz.range_mean(values, [1] * 5, 0.0, 0.5), 1.8)
        self.assertAlmostEqual(pz.range_mean(values, [1] * 5, 0.5, 1.0), 4.2)

    def test_unequal_volumes_set_the_fractions(self):
        self.assertAlmostEqual(pz.range_mean([1, 3], [1, 3], 0.0, 0.5), 2.0)


class Pu239Halves(unittest.TestCase):
    def test_one_unit_has_no_distribution(self):
        self.assertIsNone(pz.pu239_halves_ratio([2.0], [1.0]))

    def test_the_same_profile_in_more_units_gives_the_same_ratio(self):
        coarse = pz.pu239_halves_ratio([1, 2, 3, 4], [1] * 4)
        fine = pz.pu239_halves_ratio([1, 1, 2, 2, 3, 3, 4, 4], [1] * 8)
        self.assertAlmostEqual(coarse, 3.5 / 1.5)
        self.assertAlmostEqual(fine, coarse)

    def test_no_pu239_in_the_inner_half_has_no_ratio(self):
        self.assertIsNone(pz.pu239_halves_ratio([0.0, 0.0, 1.0, 1.0], [1] * 4))


class Signed(unittest.TestCase):
    def test_pcm_keeps_two_decimals_with_the_record_minus_sign(self):
        self.assertEqual(pz.signed(0.37, 2), "+0.37")
        self.assertEqual(pz.signed(-20.456, 2), "−20.46")

    def test_values_rounding_to_zero_have_no_sign(self):
        self.assertEqual(pz.signed(-0.001, 2), "0.00")


if __name__ == "__main__":
    unittest.main()
