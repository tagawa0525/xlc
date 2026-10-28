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

MATERIALS_LOG = """k 1.086000
@water=moderator:1-13: k 1.085900 (-10.0 pcm)
@water=moderator:46-69: k 1.086200 (+20.0 pcm)
"""

OTHER_LOG = """k 1.086000
O-16:1-13: k 1.085990 (-1.0 pcm)
*!1-13: k 1.085980 (-2.0 pcm)
*^1-13: k 1.086030 (+3.0 pcm)
"""

JOINT_SPEC = (
    "*:14-69,O-16:1-13,*!1-13,*^1-13,"
    "@water=moderator:1-13,@water=moderator:46-69,~water:1-13"
)

JOINT_LOG = f"""k 1.086000
*:14-69: k 1.085800 (-20.0 pcm)
O-16:1-13: k 1.085990 (-1.0 pcm)
*!1-13: k 1.085980 (-2.0 pcm)
*^1-13: k 1.086030 (+3.0 pcm)
@water=moderator:1-13: k 1.085900 (-10.0 pcm)
@water=moderator:46-69: k 1.086200 (+20.0 pcm)
~water:1-13: k 1.084800 (-120.0 pcm)
{JOINT_SPEC}: k 1.084400 (-160.0 pcm)
"""

CHI_LOG = """k 1.086000
chi: k 1.086050 (+5.0 pcm)
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

    def test_material_lines_split_each_material_into_the_four_ranges(self):
        specs = dict(sc.material_items())
        self.assertEqual(specs["material-water-fast"], "@water=moderator:1-13")
        self.assertEqual(specs["material-water-ufg"], "@water=moderator:14-27")
        self.assertEqual(specs["material-water-epithermal"], "@water=moderator:28-45")
        self.assertEqual(
            specs["material-zircaloy-thermal"], "@zircaloy=zircaloy-4:46-69"
        )
        self.assertEqual(len(specs), 8)

    def test_other_lines_split_the_fuel_and_its_non_actinides(self):
        specs = dict(sc.other_items(BURNT))
        products = ["O-16", "Xe-135", "Sm-149", "Rh-103"]
        self.assertEqual(
            specs["fast-products"], ",".join(f"{n}:1-13" for n in products)
        )
        self.assertEqual(specs["other-fast"], "*!1-13")
        self.assertEqual(
            specs["other-fast-products"], ",".join(f"{n}!1-13" for n in products)
        )
        self.assertEqual(specs["other-low"], "*!14-69")
        self.assertEqual(
            specs["other-low-products"], ",".join(f"{n}!14-69" for n in products)
        )
        self.assertEqual(specs["excess-fast"], "*^1-13")
        self.assertEqual(
            specs["excess-fast-products"], ",".join(f"{n}^1-13" for n in products)
        )

    def test_other_and_excess_lines_touch_group_one_by_their_groups(self):
        self.assertTrue(sc.touches_group_one("*!1-13"))
        self.assertTrue(sc.touches_group_one("O-16^1-13"))
        self.assertFalse(sc.touches_group_one("*!14-69"))

    def test_the_joint_line_swaps_every_separate_line_at_once(self):
        lines = sc.joint_lines(
            [("fast-u238", "U-238:1-13"), ("all", "*:14-69")],
            OTHERS,
            [
                ("material-water-fast", "@water=moderator:1-13"),
                ("material-water-thermal", "@water=moderator:46-69"),
            ],
            [("transfer-water-fast", "water", 1, 13)],
        )
        self.assertEqual(
            lines,
            [
                "*:14-69",
                "O-16:1-13",
                "*!1-13",
                "*^1-13",
                "@water=moderator:1-13",
                "@water=moderator:46-69",
                "~water:1-13",
                JOINT_SPEC,
            ],
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


ITEMS = [
    ("fast-u238", "U-238:1-13"),
    ("ufg27-u238", "U-238:27-27"),
    ("all", "*:14-69"),
]
MATERIALS = [
    ("material-water-fast", "@water=moderator:1-13"),
    ("material-water-thermal", "@water=moderator:46-69"),
]
TRANSFERS = [("transfer-water-fast", "water", 1, 13)]
# 同時の置き換えに入れない行（other-fast-products）を含む
OTHERS = [
    ("fast-products", "O-16:1-13"),
    ("other-fast", "*!1-13"),
    ("other-fast-products", "O-16!1-13"),
    ("excess-fast", "*^1-13"),
]
ALL_LOG = ABSORPTION_LOG + "*:14-69: k 1.085800 (-20.0 pcm)\n"


def run(
    k=1.085,
    absorption=ALL_LOG,
    materials=None,
    transfer=None,
    fast=True,
    other=OTHER_LOG + "O-16!1-13: k 1.085990 (-1.0 pcm)\n",
):
    return {
        "k": k,
        "k_std": 0.0002,
        "absorption": absorption,
        "materials": materials,
        "transfer": transfer,
        "other": other,
        "fast": fast,
    }


# 最後の実行が移行と燃料の吸収、最初の実行が物質の吸収を同時の置き換えに渡す
def summarize(runs, joint_log=JOINT_LOG, chi_log=CHI_LOG, materials_run=0):
    return sc.summarize_state(
        "1C 5 MWd/kg",
        ITEMS,
        OTHERS,
        MATERIALS,
        TRANSFERS,
        runs,
        {"run": len(runs) - 1, "materials_run": materials_run, "log": joint_log},
        {"log": chi_log},
    )


def two_runs():
    return [
        run(k=1.0852, materials=MATERIALS_LOG),
        run(transfer=TRANSFER_LOG),
    ]


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
        runs = [
            run(k=1.0852, materials=MATERIALS_LOG),
            run(
                absorption=ALL_LOG.replace("-50.0", "-70.0").replace(
                    "1.085500", "1.085300"
                ),
                transfer=TRANSFER_LOG,
            ),
        ]
        state = summarize(runs)
        self.assertAlmostEqual(state["diff_pcm"], 1e5 * (1.086 - 1.0851), places=6)
        self.assertAlmostEqual(state["items"]["fast-u238"]["pcm"], 60.0, places=6)
        self.assertEqual(state["items"]["fast-u238"]["runs"], [50.0, 70.0])
        self.assertAlmostEqual(
            state["items"]["transfer-water-fast"]["pcm"], 120.0, places=6
        )

    def test_material_contributions_average_the_runs_with_material_logs(self):
        runs = [
            run(k=1.0852, materials=MATERIALS_LOG),
            run(
                k=1.0851,
                materials=MATERIALS_LOG.replace("-10.0", "-14.0").replace(
                    "1.085900", "1.085860"
                ),
            ),
            run(transfer=TRANSFER_LOG),
        ]
        state = summarize(runs)
        self.assertEqual(state["items"]["material-water-fast"]["runs"], [10.0, 14.0])
        self.assertAlmostEqual(
            state["items"]["material-water-thermal"]["pcm"], -20.0, places=6
        )

    def test_runs_without_the_group_above_10_mev_leave_out_the_fast_lines(self):
        runs = [
            run(materials=MATERIALS_LOG, fast=False),
            run(
                absorption=ALL_LOG.replace("-50.0", "-70.0").replace(
                    "1.085500", "1.085300"
                ),
                materials=MATERIALS_LOG,
                transfer=TRANSFER_LOG,
            ),
        ]
        state = summarize(runs, materials_run=1)
        self.assertEqual(state["items"]["fast-u238"]["runs"], [70.0])
        self.assertEqual(state["items"]["ufg27-u238"]["runs"], [-30.0, -30.0])
        self.assertEqual(state["items"]["material-water-fast"]["runs"], [10.0])
        self.assertEqual(
            state["items"]["material-water-thermal"]["runs"], [-20.0, -20.0]
        )

    def test_other_contributions_average_the_runs_with_other_logs(self):
        runs = [
            run(
                k=1.0852,
                materials=MATERIALS_LOG,
                other=OTHER_LOG.replace("-2.0", "-6.0").replace("1.085980", "1.085940")
                + "O-16!1-13: k 1.085990 (-1.0 pcm)\n",
            ),
            run(other=None),
            run(transfer=TRANSFER_LOG),
        ]
        joint = JOINT_LOG.replace(
            "*!1-13: k 1.085980 (-2.0 pcm)", "*!1-13: k 1.085940 (-6.0 pcm)"
        )
        state = summarize(runs, joint_log=joint)
        self.assertEqual(state["items"]["other-fast"]["runs"], [6.0, 2.0])
        self.assertEqual(state["items"]["fast-products"]["runs"], [1.0, 1.0])
        self.assertEqual(state["items"]["excess-fast"]["runs"], [-3.0, -3.0])
        self.assertEqual(state["items"]["other-fast-products"]["runs"], [1.0, 1.0])
        # 同時 160 − 別々の和（吸収の全部 20 + ほかの吸収と過剰 1 + 6 − 3 + 物質 10 − 20 + 移行 120）
        self.assertAlmostEqual(state["items"]["interaction"]["pcm"], 26.0, places=6)

    def test_other_lines_of_runs_without_the_group_above_10_mev_are_left_out(self):
        runs = [
            run(fast=False),
            run(materials=MATERIALS_LOG, transfer=TRANSFER_LOG),
        ]
        state = summarize(runs, materials_run=1)
        self.assertEqual(state["items"]["other-fast"]["runs"], [2.0])

    def test_a_joint_other_line_off_the_materials_run_is_an_error(self):
        joint = JOINT_LOG.replace(
            "*^1-13: k 1.086030 (+3.0 pcm)", "*^1-13: k 1.086050 (+5.0 pcm)"
        )
        with self.assertRaises(sc.SummaryError):
            summarize(two_runs(), joint_log=joint)

    def test_the_joint_swap_needs_other_lines_on_the_materials_run(self):
        runs = [
            run(k=1.0852, materials=MATERIALS_LOG, other=None),
            run(transfer=TRANSFER_LOG),
        ]
        with self.assertRaises(sc.SummaryError):
            summarize(runs)

    def test_a_line_missing_from_a_log_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            summarize([run(absorption="k 1.086000\n", transfer=TRANSFER_LOG)])

    def test_interaction_is_the_joint_swap_less_the_separate_ones(self):
        state = summarize(two_runs())
        # 同時 160 − 別々の和（吸収の全部 20 + 物質 10 − 20 + 移行 120）
        self.assertAlmostEqual(state["items"]["interaction"]["pcm"], 30.0, places=6)
        self.assertAlmostEqual(state["items"]["chi"]["pcm"], -5.0, places=6)

    def test_a_transfer_log_outside_the_joint_run_is_an_error(self):
        # 表の移行の行は同時の置き換えの実行だけから取るので、ほかの実行の移行のログは表に入らない
        with self.assertRaisesRegex(sc.SummaryError, "other than the joint one"):
            summarize(
                [
                    run(materials=MATERIALS_LOG, transfer=TRANSFER_LOG),
                    run(transfer=TRANSFER_LOG),
                ]
            )

    def test_transfer_rows_take_the_lines_of_the_joint_log(self):
        # swap_transfer の変化は小数 6 桁の 2 つの k の差なので、同時の置き換えのログの行（同じ置き換えの 1 つの k の差）を表に使う
        joint = JOINT_LOG.replace(
            "~water:1-13: k 1.084800 (-120.0 pcm)",
            "~water:1-13: k 1.084799 (-120.1 pcm)",
        ).replace("(-160.0 pcm)", "(-160.1 pcm)")
        state = summarize(two_runs(), joint_log=joint)
        self.assertAlmostEqual(
            state["items"]["transfer-water-fast"]["pcm"], 120.1, places=6
        )
        self.assertAlmostEqual(state["items"]["interaction"]["pcm"], 30.0, places=6)

    def test_rest_is_the_difference_less_the_joint_swap_and_chi(self):
        state = summarize(two_runs())
        # 同時の置き換えの実行の差 (1.086 − 1.085) × 1e5 − 160 − (−5)
        self.assertAlmostEqual(state["items"]["rest"]["pcm"], -55.0, places=6)
        self.assertEqual(len(state["items"]["rest"]["runs"]), 1)

    def test_a_joint_log_that_disagrees_with_the_separate_logs_is_an_error(self):
        for wrong in (
            JOINT_LOG.replace(
                "~water:1-13: k 1.084800 (-120.0 pcm)",
                "~water:1-13: k 1.084700 (-130.0 pcm)",
            ),
            JOINT_LOG.replace(
                "*:14-69: k 1.085800 (-20.0 pcm)", "*:14-69: k 1.085700 (-30.0 pcm)"
            ),
            JOINT_LOG.replace(
                "@water=moderator:1-13: k 1.085900 (-10.0 pcm)",
                "@water=moderator:1-13: k 1.085800 (-20.0 pcm)",
            ),
            JOINT_LOG.replace("k 1.086000\n", "k 1.086100\n", 1),
        ):
            with (
                self.subTest(wrong=wrong.splitlines()[0]),
                self.assertRaises(sc.SummaryError),
            ):
                summarize(two_runs(), joint_log=wrong)

    def test_a_joint_log_without_the_joint_line_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            summarize(two_runs(), joint_log=JOINT_LOG.rsplit(JOINT_SPEC, 1)[0])

    def test_a_chi_log_of_another_xlc_solution_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            summarize(two_runs(), chi_log=CHI_LOG.replace("k 1.086000", "k 1.087000"))

    def test_a_transfer_log_against_another_openmc_run_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            summarize([run(k=1.0852, materials=MATERIALS_LOG, transfer=TRANSFER_LOG)])

    def test_a_transfer_log_may_differ_from_the_absorption_k_by_one_printed_digit(self):
        state = summarize(
            [
                run(k=1.0852, materials=MATERIALS_LOG),
                run(
                    transfer=TRANSFER_LOG.replace("k 1.086000 (1/A", "k 1.086001 (1/A")
                ),
            ]
        )
        self.assertEqual(len(state["items"]["transfer-water-fast"]["runs"]), 1)

    def test_a_transfer_log_of_another_xlc_solution_is_an_error(self):
        with self.assertRaises(sc.SummaryError):
            summarize(
                [
                    run(k=1.0852, materials=MATERIALS_LOG),
                    run(
                        transfer=TRANSFER_LOG.replace(
                            "k 1.086000 (1/A", "k 1.086002 (1/A"
                        )
                    ),
                ]
            )

    def test_runs_must_share_the_xlc_k(self):
        with self.assertRaises(sc.SummaryError):
            summarize(
                [
                    run(k=1.0852, materials=MATERIALS_LOG),
                    run(
                        absorption=ALL_LOG.replace("k 1.086000", "k 1.087000"),
                        transfer=TRANSFER_LOG,
                    ),
                ]
            )


if __name__ == "__main__":
    unittest.main()
