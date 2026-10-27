"""ufg_transfer.py（ADR-0207 の UFG の範囲の群の移行の形の測定の要約）。

実行: python3 -m unittest discover -s examples/depletion/issue-0133 -p 'test_ufg_transfer.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ufg_transfer as ut
from group_rows import SummaryError


def transfer_log(base: float, swapped: float, openmc: float = 0.824863) -> str:
    return f"""material water 27-27
groups 27..=27; 1 pins (1 solved); SPH factors recomputed from the reference: max |Δμ| 1.5e-12; after the swap: max |μ'/μ − 1| 6.81e-5
OpenMC k {openmc:.6f}
xlc                      k {base:.6f} (1/A {base:.6f}), k − OpenMC   +31.7 pcm; max |out-scatter xs / OpenMC − 1| 5.90e-4
swapped, SPH kept        k {swapped:.6f} (1/A {swapped:.6f}), k − OpenMC   +35.6 pcm; max |out-scatter xs / OpenMC − 1| 3.33e-16
swapped, SPH recomputed  k {swapped:.6f} (1/A {swapped:.6f}), k − OpenMC   +35.6 pcm; max |out-scatter xs / OpenMC − 1| 8.88e-16
"""


ROWS_LOG = """kernel 2A-40 material water 27-27
group 27 out/own: xlc 5.866043e-01 4.032975e-01 state 5.887940e-01 4.030164e-01 kernel-openmc 5.887341e-01 4.032593e-01 kernel-xlc 5.884444e-01 4.035452e-01 kernel-lethargy 5.945169e-01 3.976478e-01
"""

CONFIG = """[[spec]]
spec = "water:27-27"
label = "減速材、27 群"

[[kernel]]
key = "2A-40"
label = "2A-40 の核"

[[mode]]
key = "dilute"
label = "無限希釈"

[[mode]]
key = "ring"
label = "リング"

[[state]]
key = "2A-40"
name = "2A 40"
kernels = ["2A-40"]
rows = ["water:27-27"]

[[state]]
key = "1C-0"
name = "1C 0"
kernels = ["2A-40"]
"""


class Summarize(unittest.TestCase):
    def write(self, directory: Path, logs: dict[str, str]) -> Path:
        (directory / "logs").mkdir()
        for name, text in logs.items():
            (directory / "logs" / name).write_text(text)
        config = directory / "summary.toml"
        config.write_text(CONFIG)
        return config

    def logs(self) -> dict[str, str]:
        return {
            "2A-40-dilute-transfer.txt": transfer_log(0.825977, 0.825532),
            "2A-40-ring-transfer.txt": transfer_log(0.825628, 0.825500),
            "2A-40-dilute-rows.txt": ROWS_LOG,
            "2A-40-ring-rows.txt": ROWS_LOG,
            "1C-0-dilute-transfer.txt": transfer_log(1.172316, 1.172354, 1.171998),
            "1C-0-ring-transfer.txt": transfer_log(1.172000, 1.172100, 1.171998),
            "1C-0-dilute-rows.txt": "",
            "1C-0-ring-rows.txt": "",
        }

    def test_each_mode_gives_its_k_and_swaps_against_one_openmc_k(self):
        with tempfile.TemporaryDirectory() as d:
            summary = ut.summarize(self.write(Path(d), self.logs()))
        state = summary["states"][0]
        self.assertEqual(state["openmc_k"], 0.824863)
        self.assertEqual(state["modes"]["dilute"]["k"], 0.825977)
        self.assertAlmostEqual(
            state["modes"]["dilute"]["swaps"]["water:27-27"], 44.5, places=6
        )
        self.assertAlmostEqual(
            state["modes"]["ring"]["swaps"]["water:27-27"], 12.8, places=6
        )
        self.assertEqual(
            state["rows"]["ring"]["2A-40"]["water:27"]["xlc"][0], 0.5866043
        )
        # rows を比べない状態は核を持たない
        self.assertEqual(summary["states"][1]["kernels"], [])

    def test_modes_against_another_openmc_k_are_an_error(self):
        logs = self.logs()
        logs["2A-40-ring-transfer.txt"] = transfer_log(0.825628, 0.825500, 0.824000)
        with tempfile.TemporaryDirectory() as d, self.assertRaises(SummaryError):
            ut.summarize(self.write(Path(d), logs))

    def test_rows_missing_a_configured_group_are_an_error(self):
        logs = self.logs()
        logs["2A-40-ring-rows.txt"] = ""
        with tempfile.TemporaryDirectory() as d, self.assertRaises(SummaryError):
            ut.summarize(self.write(Path(d), logs))


if __name__ == "__main__":
    unittest.main()
