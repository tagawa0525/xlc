"""compare_p3.py の統計誤差の道具（OpenMC の *-std との比）。

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


class AxialReport(unittest.TestCase):
    def test_returns_the_numbers_it_prints(self):
        """差 [+0.1, −0.2]: RMS = √((0.01 + 0.04) / 2)、最大は 2 層目。AO = (上 − 下) / (上 + 下)"""
        bounds = [0.0, 1.0, 2.0]
        with contextlib.redirect_stdout(io.StringIO()):
            stats = compare_p3.axial_report(
                "ref", bounds, [1.1, 0.8], [1.0, 1.0], levels=False
            )
        self.assertAlmostEqual(stats["rms"], math.sqrt(0.025), places=12)
        self.assertAlmostEqual(
            stats["relative_rms"], math.sqrt((0.1**2 + 0.2**2) / 2), places=12
        )
        self.assertAlmostEqual(stats["max"], -0.2, places=12)
        self.assertEqual(stats["max_level"], 2)
        self.assertAlmostEqual(
            stats["ao_subject"], (0.8 - 1.1) / (0.8 + 1.1), places=12
        )
        self.assertAlmostEqual(stats["ao_reference"], 0.0, places=12)


class JsonOutput(unittest.TestCase):
    """--json の出力は、標準出力に出す数と同じ値を持つ（参照解が手元にあるときだけ）"""

    ROOT = Path(__file__).resolve().parents[2]
    P3 = ROOT / "examples/vera/p3"
    OPENMC = ROOT / "references/openmc/aws/m3/openmc-3A"
    SPEC = ROOT / "references/vera/spec.txt"

    @unittest.skipUnless(
        (OPENMC / "k.json").exists()
        and SPEC.exists()
        and (P3 / "core-3A-layers.msgpack").exists(),
        "needs the OpenMC reference and the core-3A outputs",
    )
    def test_matches_the_printed_k_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
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
                    str(self.OPENMC),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            data = json.loads(out.read_text())
        self.assertEqual(data["problem"], "3A")
        self.assertIn(f"{data['k']['xlc']:.6f}", run.stdout)
        self.assertIn(f"{data['k']['openmc']:.6f}", run.stdout)
        self.assertIn(f"{data['k']['diff_keno_pcm']:+.0f} pcm", run.stdout)
        self.assertEqual(len(data["layer_pins"]["levels"]), 49)


if __name__ == "__main__":
    unittest.main()
