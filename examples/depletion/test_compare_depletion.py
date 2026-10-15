"""compare_depletion.py の要約（燃焼点ごとの k の差と、核種の数密度の相対差）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_depletion.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_depletion


def histories(rings: int = 1):
    """2 点の燃焼の、xlc と OpenMC の履歴（体積 1 と 3 のリング）"""
    names = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149", "Am-242m"]
    volumes = [1.0, 3.0][:rings]
    xlc = {
        "format": "xlc-depletion-history",
        "nuclides": names,
        "transportable": ["U-235", "Am-242m"],
        "regions": [
            {"material": f"fuel-{r}", "volume": v} for r, v in enumerate(volumes)
        ],
        "seconds": {
            "shielding": 9.0,
            "moc": 0.5,
            "rates": 0.01,
            "matrix": 0.02,
            "cram": 0.03,
        },
        "solves": 3,
        "points": [
            {
                "burnup": b,
                "k_inf": k,
                # [リング][核種]: U-235 だけが値を持ち、リング r は 1 + r 倍
                "densities": [[d * (1 + r), 0, 0, 0, 0, 0] for r in range(rings)],
            }
            for b, k, d in [(0.0, 1.17, 2.0), (1.0, 1.12, 1.0)]
        ],
    }
    omc = {
        "format": "openmc-depletion",
        "burnup": [0.0, 1.0],
        "k": [1.171, 1.1215],
        "k-std": [0.0002, 0.0002],
        "transportable": ["Am242_m1", "U235"],
        "transport-seconds": [100.0, 110.0, 120.0],
        "total-seconds": 331.0,
        "cram-seconds": [0.5, 0.25],
        "rings": [
            {
                "name": f"fuel-{r}",
                "volume": v,
                "densities": {
                    "U235": [2.0 * (1 + r), 0.8 * (1 + r)],
                    "U238": [0.0, 0.0],
                    "Pu239": [0.0, 0.0],
                    "Xe135": [0.0, 0.0],
                    "Sm149": [0.0, 0.0],
                },
            }
            for r, v in enumerate(volumes)
        ],
    }
    return xlc, omc


class Summary(unittest.TestCase):
    def test_times_of_both_codes(self):
        s = compare_depletion.summarize(*histories())
        self.assertEqual(s["seconds"]["xlc"]["shielding"], 9.0)
        self.assertEqual(s["seconds"]["xlc"]["solves"], 3)
        o = s["seconds"]["openmc"]
        self.assertEqual(
            (o["transport"], o["fastest"], o["slowest"]), (330.0, 100.0, 120.0)
        )
        self.assertEqual((o["cram"], o["total"], o["solves"]), (0.75, 331.0, 3))

    def test_k_differences_and_their_change_from_the_first_point(self):
        s = compare_depletion.summarize(*histories())
        self.assertEqual(s["burnup"], [0.0, 1.0])
        diffs = [p["diff_pcm"] for p in s["k"]]
        self.assertAlmostEqual(diffs[0], -100.0, places=9)
        self.assertAlmostEqual(diffs[1], -150.0, places=9)
        self.assertAlmostEqual(s["k"][1]["change_pcm"], -50.0, places=9)
        self.assertAlmostEqual(s["k"][1]["openmc_std_pcm"], 20.0, places=9)

    def test_pin_densities_are_volume_averages_over_the_rings(self):
        s = compare_depletion.summarize(*histories(rings=2))
        # xlc: (1·1 + 3·2) / 4 = 1.75、OpenMC: (1·0.8 + 3·1.6) / 4 = 1.4
        self.assertAlmostEqual(s["pin"]["U-235"][1], 1.75 / 1.4 - 1.0, places=12)
        self.assertAlmostEqual(s["pin"]["U-235"][0], 0.0, places=12)
        # OpenMC が 0 の量は比べない
        self.assertIsNone(s["pin"]["Xe-135"][1])
        self.assertAlmostEqual(s["rings"][1]["U-235"], 2.0 / 1.6 - 1.0, places=12)

    def test_the_transported_nuclides_must_agree(self):
        xlc, omc = histories()
        omc["transportable"] = ["U235"]
        with self.assertRaises(SystemExit):
            compare_depletion.summarize(xlc, omc)

    def test_several_rings_pair_by_name(self):
        xlc, omc = histories(rings=2)
        xlc["regions"].reverse()
        for r, v in zip(xlc["regions"], [1.0, 3.0]):
            r["volume"] = v
        with self.assertRaises(SystemExit):
            compare_depletion.summarize(xlc, omc)

    def test_the_rings_must_have_the_same_volumes(self):
        xlc, omc = histories(rings=2)
        omc["rings"][1]["volume"] = 3.1
        with self.assertRaises(SystemExit):
            compare_depletion.summarize(xlc, omc)


if __name__ == "__main__":
    unittest.main()
