"""compare_p3.py の統計誤差のツール（OpenMC の *-std との比）。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_p3.py'
"""

import contextlib
import io
import json
import math
import subprocess
import sys
import tempfile
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


class AxialSummary(unittest.TestCase):
    def test_returns_the_numbers_and_prints_nothing(self):
        """差 [+0.1, −0.2]: RMS = √((0.01 + 0.04) / 2)、最大は 2 層目。AO = (上 − 下) / (上 + 下)"""
        bounds = [0.0, 1.0, 2.0]
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            stats = compare_p3.axial_summary(
                "ref", bounds, [1.1, 0.8], [1.0, 1.0], levels=True
            )
        self.assertEqual(out.getvalue(), "")
        self.assertAlmostEqual(stats["rms"], math.sqrt(0.025), places=12)
        self.assertAlmostEqual(
            stats["relative_rms"], math.sqrt((0.1**2 + 0.2**2) / 2), places=12
        )
        self.assertAlmostEqual(stats["max"], -0.2, places=12)
        self.assertEqual(stats["max_level"], 2)
        self.assertEqual(stats["max_bounds"], [1.0, 2.0])
        self.assertAlmostEqual(
            stats["ao_subject"], (0.8 - 1.1) / (0.8 + 1.1), places=12
        )
        self.assertAlmostEqual(stats["ao_reference"], 0.0, places=12)
        self.assertEqual(
            stats["levels"],
            [
                {"lo": 0.0, "hi": 1.0, "subject": 1.1, "reference": 1.0},
                {"lo": 1.0, "hi": 2.0, "subject": 0.8, "reference": 1.0},
            ],
        )

    def test_omits_the_level_list_when_not_asked(self):
        stats = compare_p3.axial_summary("ref", [0.0, 1.0], [1.0], [1.0], levels=False)
        self.assertNotIn("levels", stats)


class AxialReport(unittest.TestCase):
    def test_prints_the_header_and_one_line_per_level(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            compare_p3.axial_report("ref", [0.0, 1.0, 2.0], [1.1, 0.8], [1.0, 1.0])
        lines = out.getvalue().splitlines()
        self.assertTrue(lines[0].startswith("axial vs ref: rms 0.1581 "))
        self.assertEqual(
            lines[2],
            "  level  2    1.000-   2.000: xlc 0.80000 ref 1.00000 diff -0.20000",
        )
        self.assertEqual(len(lines), 3)


class Report(unittest.TestCase):
    """--json は JSON だけを書き、--report がその JSON から文章を描く（参照解が手元にあるときだけ）"""

    ROOT = Path(__file__).resolve().parents[2]
    P3 = ROOT / "examples/vera/p3"
    OPENMC = ROOT / "references/openmc/aws/m3/openmc-3A"
    SPEC = ROOT / "references/vera/spec.txt"

    def run_compare(
        self, tmp: str, *extra: str
    ) -> tuple[subprocess.CompletedProcess, dict]:
        out = Path(tmp) / "compare.json"
        run = subprocess.run(
            [
                sys.executable,
                str(Path(compare_p3.__file__)),
                "--json",
                str(out),
                "3A",
                str(self.SPEC),
                str(self.P3 / "core-3A.out"),
                str(self.P3 / "core-3A-pins.msgpack"),
                str(self.P3 / "core-3A-layers.msgpack"),
                *extra,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        return run, json.loads(out.read_text())

    needs = unittest.skipUnless(
        (OPENMC / "k.json").exists()
        and SPEC.exists()
        and (P3 / "core-3A-layers.msgpack").exists(),
        "needs the OpenMC reference and the core-3A outputs",
    )

    @needs
    def test_json_run_prints_nothing_and_report_describes_the_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            run, data = self.run_compare(tmp, str(self.OPENMC))
            self.assertEqual(run.stdout, "")
            report = subprocess.run(
                [
                    sys.executable,
                    str(Path(compare_p3.__file__)),
                    "--report",
                    str(Path(tmp) / "compare.json"),
                ],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        self.assertIn(f"{data['k']['xlc']:.6f}", report)
        self.assertIn(f"{data['k']['openmc']:.6f}", report)
        self.assertIn(f"{data['k']['diff_keno_pcm']:+.0f} pcm", report)
        self.assertEqual(len(data["layer_pins"]["levels"]), 49)
        self.assertEqual(report.count("\n  level "), 49 * 2 + 49)

    @needs
    def test_a_run_without_openmc_has_no_openmc_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, data = self.run_compare(tmp)
        self.assertNotIn("openmc", data["k"])
        self.assertNotIn("axial_vs_openmc", data)
        self.assertEqual(len(data["axial_vs_keno"]["levels"]), 49)


if __name__ == "__main__":
    unittest.main()
