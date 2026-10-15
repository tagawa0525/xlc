"""compare_aws.py の要約（区間に分けた AWS の燃焼と、分けない手元の燃焼の k と数密度の差）。

実行: python3 -m unittest discover -s examples/depletion/aws-check -p 'test_compare_aws.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_aws


def result(k, std, u235):
    """2 点の OpenMC の燃焼の結果（msgpack の形。1 リング）"""
    names = ["U235", "U238", "Pu239", "Xe135", "Sm149"]
    return {
        "format": "openmc-depletion",
        "burnup": [0.0, 1.0],
        "k": k,
        "k-std": std,
        "transportable": ["U235"],
        "rings": [
            {
                "name": "fuel",
                "volume": 2.0,
                "densities": {n: [1.0, 2.0] for n in names} | {"U235": u235},
            }
        ],
    }


class Summary(unittest.TestCase):
    def test_k_difference_is_in_sigma_of_the_combined_error(self):
        local = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws = result([1.17, 1.1203], [0.0, 3e-4], [1.0, 2.0])
        s = compare_aws.summarize(local, aws)
        self.assertAlmostEqual(s["k"][1]["diff_pcm"], 30.0, places=6)
        # 合成の σ = sqrt(3e-4^2 + 4e-4^2) = 5e-4 = 50 pcm
        self.assertAlmostEqual(s["k"][1]["sigma_pcm"], 50.0, places=6)
        self.assertAlmostEqual(s["k"][1]["ratio"], 0.6, places=6)

    def test_identical_points_have_zero_ratio_without_dividing_by_zero(self):
        r = result([1.17, 1.12], [0.0, 0.0], [1.0, 2.0])
        s = compare_aws.summarize(r, r)
        self.assertEqual([k["ratio"] for k in s["k"]], [0.0, 0.0])

    def test_density_difference_is_relative_to_the_local_run(self):
        local = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.002])
        s = compare_aws.summarize(local, aws)
        self.assertTrue(math.isclose(s["density"]["U235"][1], 1e-3, rel_tol=1e-9))
        self.assertEqual(s["density"]["U235"][0], 0.0)

    def test_density_is_none_when_the_local_run_has_zero(self):
        local = result([1.17, 1.12], [3e-4, 4e-4], [0.0, 2.0])
        s = compare_aws.summarize(local, local)
        self.assertIsNone(s["density"]["U235"][0])

    def test_different_burnup_points_stop(self):
        local = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws["burnup"] = [0.0, 2.0]
        with self.assertRaises(SystemExit):
            compare_aws.summarize(local, aws)

    def test_different_transported_nuclides_stop(self):
        local = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws = result([1.17, 1.12], [3e-4, 4e-4], [1.0, 2.0])
        aws["transportable"] = ["U238"]
        with self.assertRaises(SystemExit):
            compare_aws.summarize(local, aws)


if __name__ == "__main__":
    unittest.main()
