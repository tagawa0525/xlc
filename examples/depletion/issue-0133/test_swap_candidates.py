"""swap_candidates.py（#0133 の燃焼で増える k の差の候補ごとの、置き換えで測る k への寄与）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133 -p 'test_swap_candidates.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import swap_candidates as sc

FRESH = ["U-234", "U-235", "U-236", "U-238", "O-16"]
BURNT = [*FRESH, "Pu-239", "Pu-241", "Xe-135", "Sm-149", "Rh-103", "Am-241"]

ABSORPTION_LOG = """k 1.086000
U-238:1-13: k 1.085500 (-50.0 pcm)
U-238:27-27: k 1.086300 (+30.0 pcm)
"""

TRANSFER_LOG = """material water
groups 1..=13; 1 pins (1 solved); SPH factors recomputed from the reference: max |Δμ| 1.0e-9; after the swap: max |μ'/μ − 1| 1.00e-3
OpenMC k 1.085000
xlc                      k 1.086000 (1/A 1.086000), k − OpenMC  +100.0 pcm; max |out-scatter xs / OpenMC − 1| 5.00e-2
swapped, SPH kept        k 1.084800 (1/A 1.084800), k − OpenMC   -20.0 pcm; max |out-scatter xs / OpenMC − 1| 1.00e-12
swapped, SPH recomputed  k 1.084700 (1/A 1.084700), k − OpenMC   -30.0 pcm; max |out-scatter xs / OpenMC − 1| 1.00e-12
"""


class Lists(unittest.TestCase):
    def test_fresh_fuel_has_no_plutonium_or_fission_product_lines(self):
        keys = [key for key, _ in sc.absorption_items(FRESH)]
        self.assertNotIn("thermal-pu", keys)
        self.assertNotIn("fp", keys)
        self.assertNotIn("xe-sm", keys)
        self.assertIn("ufg27-u238", keys)

    def test_other_nuclides_leave_out_u238_and_oxygen(self):
        specs = dict(sc.absorption_items(BURNT))
        self.assertEqual(
            specs["ufg27-other"],
            ",".join(
                f"{n}:27-27"
                for n in [
                    "U-234",
                    "U-235",
                    "U-236",
                    "Pu-239",
                    "Pu-241",
                    "Xe-135",
                    "Sm-149",
                    "Rh-103",
                    "Am-241",
                ]
            ),
        )

    def test_fission_products_are_the_nuclides_lighter_than_actinides(self):
        specs = dict(sc.absorption_items(BURNT))
        self.assertEqual(specs["fp"], "Xe-135:14-69,Sm-149:14-69,Rh-103:14-69")

    def test_fast_lines_hold_only_actinides(self):
        specs = dict(sc.absorption_items(BURNT))
        self.assertEqual(specs["fast-u238"], "U-238:1-13")
        self.assertEqual(
            specs["fast-other"],
            "U-234:1-13,U-235:1-13,U-236:1-13,Pu-239:1-13,Pu-241:1-13,Am-241:1-13",
        )


class Logs(unittest.TestCase):
    def test_absorption_log_gives_the_change_of_k_per_line(self):
        k0, changes = sc.parse_absorption(ABSORPTION_LOG)
        self.assertEqual(k0, 1.086)
        self.assertAlmostEqual(changes["U-238:1-13"], -50.0, places=6)
        self.assertAlmostEqual(changes["U-238:27-27"], 30.0, places=6)

    def test_transfer_log_uses_the_kept_sph_factors(self):
        found = sc.transfer_runs(TRANSFER_LOG)
        self.assertEqual(found[("water", 1, 13)]["openmc"], 1.085)
        self.assertAlmostEqual(found[("water", 1, 13)]["change"], -120.0, places=6)


ITEMS = [("fast-u238", "U-238:1-13"), ("ufg27-u238", "U-238:27-27")]
TRANSFERS = [("transfer-water-fast", "water", 1, 13)]


class Malformed(unittest.TestCase):
    def test_a_repeated_absorption_line_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.parse_absorption(ABSORPTION_LOG + "U-238:1-13: k 1.085500 (-50.0 pcm)\n")

    def test_a_repeated_transfer_block_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.transfer_runs(TRANSFER_LOG + TRANSFER_LOG)

    def test_an_unreadable_groups_line_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.transfer_runs(TRANSFER_LOG.replace("groups 1..=13;", "groups 1 to 13;"))

    def test_a_block_without_its_own_xlc_k_is_an_error(self):
        second = "\n".join(
            line
            for line in TRANSFER_LOG.splitlines()
            if not line.startswith(("material", "xlc "))
        ).replace("groups 1..=13", "groups 14..=27")
        with self.assertRaises(sc.SummaryError):
            sc.transfer_runs(TRANSFER_LOG + second + "\n")


class Summary(unittest.TestCase):
    def test_contribution_is_minus_the_change_averaged_over_runs(self):
        state = sc.summarize_state(
            "1C 5 MWd/kg",
            ITEMS,
            TRANSFERS,
            [
                {
                    "k": 1.0850,
                    "k_std": 0.0002,
                    "absorption": ABSORPTION_LOG,
                    "transfer": TRANSFER_LOG,
                    "fast": True,
                },
                {
                    "k": 1.0852,
                    "k_std": 0.0002,
                    "absorption": ABSORPTION_LOG.replace("-50.0", "-70.0").replace(
                        "1.085500", "1.085300"
                    ),
                    "transfer": None,
                    "fast": True,
                },
            ],
        )
        self.assertAlmostEqual(state["diff_pcm"], 1e5 * (1.086 - 1.0851), places=6)
        self.assertAlmostEqual(state["items"]["fast-u238"]["pcm"], 60.0, places=6)
        self.assertEqual(state["items"]["fast-u238"]["runs"], [50.0, 70.0])
        self.assertAlmostEqual(
            state["items"]["transfer-water-fast"]["pcm"], 120.0, places=6
        )

    def test_runs_without_the_group_above_10_mev_leave_out_the_fast_lines(self):
        state = sc.summarize_state(
            "2A 40 MWd/kg",
            ITEMS,
            TRANSFERS,
            [
                {
                    "k": 1.0850,
                    "k_std": 0.0002,
                    "absorption": ABSORPTION_LOG,
                    "transfer": None,
                    "fast": False,
                },
                {
                    "k": 1.0850,
                    "k_std": 0.0002,
                    "absorption": ABSORPTION_LOG.replace("-50.0", "-70.0").replace(
                        "1.085500", "1.085300"
                    ),
                    "transfer": None,
                    "fast": True,
                },
            ],
        )
        self.assertEqual(state["items"]["fast-u238"]["runs"], [70.0])
        self.assertEqual(state["items"]["ufg27-u238"]["runs"], [-30.0, -30.0])

    def test_a_line_missing_from_a_log_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.summarize_state(
                "1C 0 MWd/kg",
                ITEMS,
                TRANSFERS,
                [
                    {
                        "k": 1.0,
                        "k_std": 0.0002,
                        "absorption": "k 1.000000\n",
                        "transfer": None,
                        "fast": True,
                    }
                ],
            )

    def test_rest_pairs_each_run_with_its_own_openmc_k(self):
        items = [*ITEMS, ("all", "*:14-69")]
        log = ABSORPTION_LOG + "*:14-69: k 1.085800 (-20.0 pcm)\n"
        state = sc.summarize_state(
            "1C 5 MWd/kg",
            items,
            TRANSFERS,
            [
                {
                    "k": 1.0852,
                    "k_std": 0.0002,
                    "absorption": log,
                    "transfer": None,
                    "fast": True,
                },
                {
                    "k": 1.0850,
                    "k_std": 0.0002,
                    "absorption": log,
                    "transfer": TRANSFER_LOG,
                    "fast": True,
                },
            ],
        )
        # 移行のある実行だけ: (1.086 − 1.085) × 1e5 − 20 − 120
        self.assertAlmostEqual(state["items"]["rest"]["pcm"], -40.0, places=6)
        self.assertEqual(len(state["items"]["rest"]["runs"]), 1)

    def test_a_transfer_log_against_another_openmc_run_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.summarize_state(
                "1C 5 MWd/kg",
                ITEMS,
                TRANSFERS,
                [
                    {
                        "k": 1.0852,
                        "k_std": 0.0002,
                        "absorption": ABSORPTION_LOG,
                        "transfer": TRANSFER_LOG,
                        "fast": True,
                    }
                ],
            )

    def test_a_transfer_log_may_differ_from_the_absorption_k_by_one_printed_digit(self):
        state = sc.summarize_state(
            "1C 5 MWd/kg",
            ITEMS,
            TRANSFERS,
            [
                {
                    "k": 1.085,
                    "k_std": 0.0002,
                    "absorption": ABSORPTION_LOG,
                    "transfer": TRANSFER_LOG.replace(
                        "k 1.086000 (1/A", "k 1.086001 (1/A"
                    ),
                    "fast": True,
                }
            ],
        )
        self.assertEqual(len(state["items"]["transfer-water-fast"]["runs"]), 1)

    def test_a_transfer_log_of_another_xlc_solution_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            sc.summarize_state(
                "1C 5 MWd/kg",
                ITEMS,
                TRANSFERS,
                [
                    {
                        "k": 1.085,
                        "k_std": 0.0002,
                        "absorption": ABSORPTION_LOG,
                        "transfer": TRANSFER_LOG.replace(
                            "k 1.086000 (1/A", "k 1.086002 (1/A"
                        ),
                        "fast": True,
                    }
                ],
            )

    def test_runs_must_share_the_xlc_k(self):
        with self.assertRaises(sc.SummaryError):
            sc.summarize_state(
                "1C 0 MWd/kg",
                ITEMS,
                TRANSFERS,
                [
                    {
                        "k": 1.0,
                        "k_std": 0.0002,
                        "absorption": ABSORPTION_LOG,
                        "transfer": None,
                        "fast": True,
                    },
                    {
                        "k": 1.0,
                        "k_std": 0.0002,
                        "absorption": ABSORPTION_LOG.replace(
                            "k 1.086000", "k 1.087000"
                        ),
                        "transfer": None,
                        "fast": True,
                    },
                ],
            )


if __name__ == "__main__":
    unittest.main()
