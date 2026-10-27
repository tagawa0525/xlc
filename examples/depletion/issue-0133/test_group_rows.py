"""group_rows.py（#0152 の群の移行の行を核と群の中の重みに分けた測定の要約）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133 -p 'test_group_rows.py'
"""

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import group_rows as gr

ROWS_LOG = """kernel 2A-40 material water 27-27
group 27 out/own: xlc 5.932341e-01 3.967022e-01 state 5.935841e-01 3.982400e-01 kernel-openmc 5.887341e-01 4.032593e-01 kernel-xlc 5.938165e-01 3.982808e-01 kernel-lethargy 5.945169e-01 3.976478e-01
kernel MOX1 material water 27-27
group 27 out/own: xlc 5.932341e-01 3.967022e-01 state 5.935841e-01 3.982400e-01 kernel-openmc 5.906456e-01 4.014076e-01 kernel-xlc 5.938789e-01 3.982091e-01 kernel-lethargy 5.946445e-01 3.975824e-01
"""

SPECTRUM_LOG = """material water
7.879-9.877 eV: xlc 0.28915 OpenMC 0.28619 比 +0.0104
4-5.014 eV: xlc 0.25439 OpenMC 0.25327 比 +0.0044
material fuel
7.879-9.877 eV: xlc 0.31304 OpenMC 0.31151 比 +0.0049
4-5.014 eV: xlc 0.27385 OpenMC 0.27451 比 -0.0024
"""


def transfer_log(base: float, swapped: float, material="water", groups="27..=27"):
    first, last = groups.split("..=")
    return f"""material {material} {first}-{last}
groups {groups}; 1 pins (1 solved); SPH factors recomputed from the reference: max |Δμ| 1.5e-12; after the swap: max |μ'/μ − 1| 6.81e-5
OpenMC k 1.171998
xlc                      k {base:.6f} (1/A {base:.6f}), k − OpenMC   +31.7 pcm; max |out-scatter xs / OpenMC − 1| 5.90e-4
swapped, SPH kept        k {swapped:.6f} (1/A {swapped:.6f}), k − OpenMC   +35.6 pcm; max |out-scatter xs / OpenMC − 1| 3.33e-16
swapped, SPH recomputed  k {swapped:.6f} (1/A {swapped:.6f}), k − OpenMC   +35.6 pcm; max |out-scatter xs / OpenMC − 1| 8.88e-16
"""


class Rows(unittest.TestCase):
    def test_rows_give_the_out_and_own_transfer_of_each_weight_by_kernel(self):
        rows = gr.parse_rows(ROWS_LOG)
        row = rows["2A-40"]["water:27"]
        self.assertEqual(row["xlc"], [0.5932341, 0.3967022])
        self.assertEqual(row["kernel-xlc"], [0.5938165, 0.3982808])
        self.assertEqual(rows["MOX1"]["water:27"]["kernel-openmc"][0], 0.5906456)

    def test_a_group_line_before_its_kernel_is_an_error(self):
        with self.assertRaises(gr.SummaryError):
            gr.parse_rows(ROWS_LOG.split("\n", 1)[1])


class Spectrum(unittest.TestCase):
    def test_spectrum_gives_the_ratio_by_material_from_low_energy(self):
        spectrum = gr.parse_spectrum(SPECTRUM_LOG)
        self.assertEqual(
            spectrum["water"][0], {"low": 4.0, "high": 5.014, "ratio": 0.0044}
        )
        self.assertEqual(spectrum["fuel"][1]["ratio"], 0.0049)

    def test_an_unreadable_line_is_an_error(self):
        with self.assertRaises(gr.SummaryError):
            gr.parse_spectrum("material water\n4-5 eV: ratio\n")


class Summarize(unittest.TestCase):
    def test_contributions_are_the_swapped_changes_with_the_sign_flipped(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "logs").mkdir()
            (d / "summary.toml").write_text(
                """[[spec]]
spec = "water:27-27"
label = "減速材、27 群"

[[kernel]]
key = "2A-40"
label = "2A-40 の核"

[[state]]
key = "1C-0"
name = "1C 0"
kernels = ["2A-40"]
"""
            )
            (d / "logs/1C-0-transfer.txt").write_text(transfer_log(1.172316, 1.172354))
            (d / "logs/1C-0-2A-40-xlc-transfer.txt").write_text(
                transfer_log(1.172316, 1.172386)
            )
            (d / "logs/1C-0-2A-40-lethargy-transfer.txt").write_text(
                transfer_log(1.172316, 1.172473)
            )
            (d / "logs/1C-0-rows.txt").write_text(ROWS_LOG.split("kernel MOX1")[0])
            summary = gr.summarize(d / "summary.toml")
        json.dumps(summary)
        (state,) = summary["states"]
        swap = state["swaps"]["water:27-27"]
        self.assertAlmostEqual(swap["openmc"], -3.8)
        self.assertAlmostEqual(swap["kernels"]["2A-40"]["xlc"], -7.0)
        self.assertAlmostEqual(swap["kernels"]["2A-40"]["lethargy"], -15.7)
        self.assertEqual(state["k"], 1.172316)
        self.assertEqual(state["openmc_k"], 1.171998)
        self.assertEqual(state["spectrum"], None)
        self.assertEqual(state["rows"]["2A-40"]["water:27"]["state"][0], 0.5935841)

    def test_logs_from_another_unswapped_k_are_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "logs").mkdir()
            (d / "summary.toml").write_text(
                """[[spec]]
spec = "water:27-27"
label = "減速材"

[[kernel]]
key = "2A-40"
label = "核"

[[state]]
key = "1C-0"
name = "1C 0"
kernels = ["2A-40"]
"""
            )
            (d / "logs/1C-0-transfer.txt").write_text(transfer_log(1.172316, 1.172354))
            (d / "logs/1C-0-2A-40-xlc-transfer.txt").write_text(
                transfer_log(1.172000, 1.172386)
            )
            (d / "logs/1C-0-2A-40-lethargy-transfer.txt").write_text(
                transfer_log(1.172316, 1.172473)
            )
            (d / "logs/1C-0-rows.txt").write_text(ROWS_LOG)
            with self.assertRaises(gr.SummaryError):
                gr.summarize(d / "summary.toml")


if __name__ == "__main__":
    unittest.main()
