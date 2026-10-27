"""ufg_free_gas.py（ADR-0203 の UFG の自由ガスの散乱核の測定の要約）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133 -p 'test_ufg_free_gas.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ufg_free_gas as ug

ABSORPTION_LOG = """k 0.825540
U-238:27-27: k 0.825879 (+33.9 pcm)
*:27-27: k 0.825700 (+16.0 pcm)
U-238:14-26: k 0.825823 (+28.3 pcm)
*:14-26: k 0.825800 (+26.0 pcm)
"""

FINE_LOG = """4.0-9.877 eV: 195 bins
  4.018-5.012 eV: xlc 0.27630 OpenMC 0.27928 比 -0.0107
  5.012-6.281 eV: xlc 0.21730 OpenMC 0.21941 比 -0.0096
  6.281-7.87 eV: xlc 0.18403 OpenMC 0.18235 比 +0.0092
  7.87-9.863 eV: xlc 0.32237 OpenMC 0.31896 比 +0.0107
  U-236: 区間の断面積を xlc の中性子束で重みづけた群の捕獲 / OpenMC − 1 = -0.0170
"""


class Absorption(unittest.TestCase):
    def test_contributions_are_the_swapped_changes_with_the_sign_flipped(self):
        k, swaps = ug.parse_absorption(ABSORPTION_LOG)
        self.assertEqual(k, 0.825540)
        self.assertEqual(swaps["U-238:27-27"], -33.9)
        self.assertEqual(swaps["*:14-26"], -26.0)

    def test_a_log_without_the_unswapped_k_is_an_error(self):
        with self.assertRaises(ug.SummaryError):
            ug.parse_absorption("U-238:27-27: k 0.825879 (+33.9 pcm)\n")


class Fine(unittest.TestCase):
    def test_bins_give_the_bounds_and_the_ratio(self):
        bins = ug.parse_fine(FINE_LOG)
        self.assertEqual(len(bins), 4)
        self.assertEqual(bins[0], {"low": 4.018, "high": 5.012, "ratio": -0.0107})
        self.assertEqual(bins[3]["ratio"], 0.0107)

    def test_a_log_without_bins_is_an_error(self):
        with self.assertRaises(ug.SummaryError):
            ug.parse_fine("4.0-9.877 eV: 195 bins\n")


SPECS = ["U-238:27-27", "*:27-27", "U-238:14-26", "*:14-26"]

TRANSFER_LOG = """material water
groups 27..=27; 1 pins (1 solved); SPH factors recomputed from the reference: max |Δμ| 3.6e-12
OpenMC k 1.171998
xlc                      k 1.172193 (1/A 1.172193), k − OpenMC   +19.5 pcm; max |out-scatter xs / OpenMC − 1| 3.52e-2
swapped, SPH kept        k 1.172100 (1/A 1.172100), k − OpenMC   +10.2 pcm; max |out-scatter xs / OpenMC − 1| 1.22e-15
swapped, SPH recomputed  k 1.172101 (1/A 1.172101), k − OpenMC   +10.3 pcm; max |out-scatter xs / OpenMC − 1| 6.66e-16
material fuel
groups 14..=26; 1 pins (1 solved); SPH factors recomputed from the reference: max |Δμ| 3.6e-12
OpenMC k 1.171998
xlc                      k 1.172193 (1/A 1.172193), k − OpenMC   +19.5 pcm; max |out-scatter xs / OpenMC − 1| 1.77e-1
swapped, SPH kept        k 1.172243 (1/A 1.172243), k − OpenMC   +24.5 pcm; max |out-scatter xs / OpenMC − 1| 3.63e-7
swapped, SPH recomputed  k 1.172240 (1/A 1.172240), k − OpenMC   +24.2 pcm; max |out-scatter xs / OpenMC − 1| 2.22e-8
"""


# 移行のログと同じ置き換えない k の吸収のログ
SAME_K = ABSORPTION_LOG.replace("k 0.825540", "k 1.172193", 1)


class Transfer(unittest.TestCase):
    def test_contributions_are_the_changes_with_the_sph_kept_and_the_sign_flipped(self):
        logs = {"static": (SAME_K, None, TRANSFER_LOG)}
        state = ug.summarize_state("1C 0", "1C-0", [0.825], ["static"], logs, SPECS)
        transfers = state["results"]["static"]["transfers"]
        self.assertEqual(list(transfers), ["water:27-27", "fuel:14-26"])
        self.assertAlmostEqual(transfers["water:27-27"], (1.172193 - 1.172100) * 1e5)
        self.assertAlmostEqual(transfers["fuel:14-26"], (1.172193 - 1.172243) * 1e5)

    def test_kernels_without_a_transfer_log_have_none(self):
        logs = {"static": (ABSORPTION_LOG, None, None)}
        state = ug.summarize_state("1C 0", "1C-0", [0.825], ["static"], logs, SPECS)
        self.assertIsNone(state["results"]["static"]["transfers"])

    def test_a_transfer_log_of_another_unswapped_k_is_an_error(self):
        # 移行の置き換えは、同じ散乱核の吸収のログと同じ置き換えない k から測ったもの
        other = TRANSFER_LOG.replace("k 1.172193", "k 1.172000")
        logs = {"static": (SAME_K, None, other)}
        with self.assertRaises(ug.SummaryError):
            ug.summarize_state("1C 0", "1C-0", [0.825], ["static"], logs, SPECS)


class State(unittest.TestCase):
    def test_kernels_without_logs_are_left_out(self):
        logs = {"static": (ABSORPTION_LOG, None, None), "exact": None}
        state = ug.summarize_state(
            "2A 40", "2A-40", [0.825], ["static", "exact"], logs, SPECS
        )
        self.assertEqual(list(state["results"]), ["static"])

    def test_differences_use_the_mean_openmc_k_and_the_static_k(self):
        logs = {
            "static": (ABSORPTION_LOG, None, None),
            "exact": (
                ABSORPTION_LOG.replace("k 0.825540", "k 0.825835", 1),
                FINE_LOG,
                None,
            ),
        }
        state = ug.summarize_state(
            "2A 40", "2A-40", [0.82480, 0.82500], ["static", "exact"], logs, SPECS
        )
        static, exact = state["results"]["static"], state["results"]["exact"]
        self.assertAlmostEqual(static["diff_pcm"], (0.825540 - 0.82490) * 1e5)
        self.assertAlmostEqual(exact["change_pcm"], (0.825835 - 0.825540) * 1e5)
        self.assertIsNone(static["tilt"])
        self.assertEqual(len(exact["tilt"]), 4)

    def test_the_static_mode_must_come_first(self):
        logs = {"exact": (ABSORPTION_LOG, None, None)}
        with self.assertRaises(ug.SummaryError):
            ug.summarize_state("2A 40", "2A-40", [0.825], ["exact"], logs, SPECS)

    def test_modes_must_swap_the_same_list(self):
        logs = {
            "static": (ABSORPTION_LOG, None, None),
            "exact": ("k 0.825835\nU-238:27-27: k 0.826 (+16.5 pcm)\n", None, None),
        }
        with self.assertRaises(ug.SummaryError):
            ug.summarize_state(
                "2A 40", "2A-40", [0.825], ["static", "exact"], logs, SPECS
            )

    def test_the_state_keeps_its_declared_swaps(self):
        logs = {"static": (ABSORPTION_LOG, None, None)}
        state = ug.summarize_state("2A 40", "2A-40", [0.825], ["static"], logs, SPECS)
        self.assertEqual(state["swaps"], SPECS)

    def test_logs_must_swap_the_declared_list(self):
        # 途中で切れたログや、宣言と違う一覧で流したログは止める
        logs = {"static": (ABSORPTION_LOG, None, None)}
        with self.assertRaises(ug.SummaryError):
            ug.summarize_state(
                "2A 40", "2A-40", [0.825], ["static"], logs, [*SPECS, "Pu-239:27-27"]
            )

    def test_swap_lists_skip_comments_and_blank_lines(self):
        text = "# 説明\nU-238:27-27\n\n*:27-27\n"
        self.assertEqual(ug.parse_swaps(text), ["U-238:27-27", "*:27-27"])


if __name__ == "__main__":
    unittest.main()
