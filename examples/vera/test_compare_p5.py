"""compare_p5.py の要約（JSON）と、JSON から描く文章の報告。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_p5.py'
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

import msgpack
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_p5

N = 17


def core_out(path: Path, k: float, rates: np.ndarray) -> None:
    """core の標準出力（k-eff の行と、反射体の位置 `.` で囲んだ集合体の核分裂率）"""
    rows = [f"k-eff {k:.8f}", "(`.` is a reflector site, `-` is empty)"]
    padded = np.full((rates.shape[0] + 2, rates.shape[1] + 2), np.nan)
    padded[1:-1, 1:-1] = rates
    for row in padded:
        rows.append(" ".join("." if np.isnan(x) else f"{x:.5f}" for x in row))
    path.write_text("\n".join(rows) + "\n")


def pins_map(assemblies: np.ndarray, tilt: float) -> np.ndarray:
    """集合体 3 × 3 の 17 × 17 のピン（燃料のセルは一様に集合体の値。tilt で左右に傾ける）"""
    pins = np.kron(assemblies, np.ones((N, N)))
    pins *= 1 + tilt * np.linspace(-1, 1, pins.shape[1])[None, :]
    return pins


def write_rows(path: Path, rows: np.ndarray) -> None:
    path.write_bytes(msgpack.packb({"format": "x", "rows": rows.tolist()}))


# 付録 J の 1/8 の表の行ごとの集合体の数（compare_p5.OCTANT_ROWS）
def spec_text() -> str:
    """5A-2D の k、1/8 の集合体出力の表（値は行と列から決めた 5 桁）、表 M3-6 の最大のピン出力だけの仕様書"""
    lines = ["5A-2D 1.004085 0.000008", "", "5A-2D Assembly Powers (ENDF/B-VII.0)", ""]
    for row, count in enumerate(compare_p5.OCTANT_ROWS):
        lines.append(
            " ".join(f"{0.8 + 0.01 * row + 0.02 * c:.5f}" for c in range(count))
        )
    lines += [
        "",
        compare_p5.POWER_TABLE,
        "Maximum 1.2 1.3 1.4",
        "Maximum Uncertainty 0.1 0.1 0.1",
        "Maximum 1.5 1.6 1.7",
    ]
    return "\n".join(lines) + "\n"


class Fixture(unittest.TestCase):
    """3 × 3 の集合体の炉心。xlc は OpenMC より、中心を除いて 2% 高い"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        reference = np.array([[1.0, 1.0, 1.0], [1.0, 0.5, 1.0], [1.0, 1.0, 1.0]])
        ours = reference * np.where(reference == 0.5, 1.0, 1.02)
        self.reference, self.ours = reference, ours
        core_out(self.dir / "core.out", 1.0010, ours)
        write_rows(self.dir / "core-pins.msgpack", pins_map(ours, 0.0))
        results = self.dir / "results"
        results.mkdir()
        (results / "k.json").write_text(json.dumps({"k": 1.0000, "k-std": 2e-5}))
        write_rows(results / "assemblies.msgpack", reference)
        write_rows(results / "assemblies-std.msgpack", reference * 0.005)
        pins = pins_map(reference, 0.0)
        write_rows(results / "pins.msgpack", pins)
        write_rows(results / "pins-std.msgpack", pins * 0.002)
        (self.dir / "spec.txt").write_text(spec_text())

    def run_tool(self, *args):
        return subprocess.run(
            [sys.executable, str(HERE / "compare_p5.py"), *map(str, args)],
            capture_output=True,
            text=True,
            check=False,
        )


class AssemblySummary(unittest.TestCase):
    def test_returns_the_numbers_and_prints_nothing(self):
        """xlc [1.1, 0.9]、参照 [1, 1]（平均で規格化して同じ）: 相対差 +10% と −10%、RMS 10%"""
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            d = compare_p5.assembly_summary(
                np.array([[1.1, 0.9]]), np.array([[1.0, 1.0]])
            )
        self.assertEqual(out.getvalue(), "")
        self.assertAlmostEqual(d["relative"]["rms"], 0.1, places=12)
        self.assertAlmostEqual(d["relative"]["max"], 0.1, places=12)
        self.assertEqual(d["relative"]["at"], [0, 0])
        self.assertEqual(d["relative"]["n"], 2)

    def test_counts_the_differences_in_sigma_of_the_reference(self):
        """差 +0.1 を標準偏差 0.06 で割って約 1.7 倍、−0.1 を 0.04 で割って約 2.5 倍（2σ を超えるのは 1 個）"""
        d = compare_p5.assembly_summary(
            np.array([[1.1, 0.9]]), np.array([[1.0, 1.0]]), np.array([[0.06, 0.04]])
        )
        self.assertAlmostEqual(d["sigma"]["max"], 2.5, places=12)
        self.assertEqual(d["sigma"]["at"], [0, 1])
        self.assertEqual(d["sigma"]["beyond_2sigma"], 1)
        self.assertEqual(d["sigma"]["n"], 2)

    def test_names_the_position_of_the_largest_difference(self):
        d = compare_p5.assembly_summary(np.array([[1.1, 0.9]]), np.array([[1.0, 1.0]]))
        self.assertEqual(d["relative"]["name"], "R01")
        self.assertEqual(d["absolute"]["name"], "R01")

    def test_a_zero_standard_deviation_stops_instead_of_writing_inf(self):
        with self.assertRaises(SystemExit):
            compare_p5.assembly_summary(
                np.array([[1.1, 0.9]]), np.array([[1.0, 1.0]]), np.array([[0.0, 0.1]])
            )


class Names(unittest.TestCase):
    def test_assembly_names_use_the_columns_of_the_core_map(self):
        """列は西から R P N M L K J H G F E D C B A、行は北から 01。(1, 11) は D02"""
        self.assertEqual(compare_p5.assembly_name([1, 11]), "D02")
        self.assertEqual(compare_p5.assembly_name([7, 7]), "H08")


class OpenmcJson(Fixture):
    def test_json_has_the_numbers_and_prints_nothing(self):
        run = self.run_tool(
            "--json",
            self.dir / "out.json",
            "openmc",
            "5A",
            self.dir / "spec.txt",
            self.dir / "core.out",
            self.dir / "core-pins.msgpack",
            self.dir / "results",
        )
        self.assertEqual((run.returncode, run.stdout), (0, ""), run.stderr)
        s = json.loads((self.dir / "out.json").read_text())
        self.assertEqual(s["format"], "xlc-compare-p5")
        self.assertEqual(s["kind"], "openmc")
        k = s["k"]
        self.assertAlmostEqual(k["diff_pcm"], 1e5 * (1.0010 / 1.0 - 1), places=6)
        self.assertAlmostEqual(k["diff_std_pcm"], 2.0, places=9)
        # OpenMC − KENO-VI（絶対差 × 1e5）と、両方の統計誤差の合成
        self.assertAlmostEqual(k["keno"], 1.004085, places=12)
        self.assertAlmostEqual(k["diff_keno_pcm"], 1e5 * (1.0 - 1.004085), places=6)
        self.assertAlmostEqual(
            k["diff_keno_std_pcm"], 1e5 * math.hypot(2e-5, 8e-6), places=9
        )
        self.assertEqual(
            sorted(i["path"].rsplit("/", 1)[-1] for i in s["inputs"]),
            sorted(
                [
                    "spec.txt",
                    "core.out",
                    "core-pins.msgpack",
                    "k.json",
                    "assemblies.msgpack",
                    "assemblies-std.msgpack",
                    "pins.msgpack",
                    "pins-std.msgpack",
                ]
            ),
        )
        self.assertEqual(s["assemblies"]["relative"]["n"], 9)
        self.assertEqual(len(s["assemblies"]["quarter"]), 2)
        self.assertEqual(s["pins"]["all"]["n"], 9 * N * N)

    def test_report_draws_the_k_line_from_the_json(self):
        first = self.run_tool(
            "--json",
            self.dir / "out.json",
            "openmc",
            "5A",
            self.dir / "spec.txt",
            self.dir / "core.out",
            self.dir / "core-pins.msgpack",
            self.dir / "results",
        )
        self.assertEqual(first.returncode, 0, first.stderr)
        run = self.run_tool("--report", self.dir / "out.json")
        lines = run.stdout.splitlines()
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertTrue(lines[0].startswith("inputs: "))
        self.assertEqual(
            lines[1],
            "k: xlc 1.001000, OpenMC 1.000000 +- 0.000020, +100.0 pcm",
        )
        self.assertTrue(
            any(l.startswith("assemblies (relative): 9, rms") for l in lines)
        )
        self.assertTrue(
            math.isfinite(
                json.loads((self.dir / "out.json").read_text())["k"]["diff_pcm"]
            )
        )

    def openmc_args(self, problem="5A"):
        return [
            "--json",
            self.dir / "out.json",
            "openmc",
            problem,
            self.dir / "spec.txt",
            self.dir / "core.out",
            self.dir / "core-pins.msgpack",
            self.dir / "results",
        ]

    def test_an_unknown_problem_name_is_not_written(self):
        run = self.run_tool(*self.openmc_args("9Z"))
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse((self.dir / "out.json").exists())

    def test_a_zero_openmc_k_std_is_not_written_as_inf(self):
        (self.dir / "results" / "k.json").write_text(
            json.dumps({"k": 1.0, "k-std": 0.0})
        )
        run = self.run_tool(*self.openmc_args())
        self.assertNotEqual(run.returncode, 0)
        self.assertFalse((self.dir / "out.json").exists())


class KenoJson(Fixture):
    """KENO-VI との比較。参照の集合体出力は仕様書の 1/8 の表から広げたもの、xlc はそれを 1% 高くして 1 つだけ 3% 高くしたもの"""

    def test_json_has_the_k_difference_and_the_largest_assembly(self):
        spec = spec_text().splitlines()
        reference = compare_p5.keno_assemblies(spec, "5A")
        ours = reference * 1.01
        i, j = (int(x) for x in np.argwhere(~np.isnan(reference))[40])
        ours[i, j] = reference[i, j] * 1.03
        ours = np.where(np.isnan(reference), np.nan, ours)
        core_out(self.dir / "core.out", 1.004085 * 1.001, ours)
        write_rows(self.dir / "core-pins.msgpack", np.ones((17, 17)))
        run = self.run_tool(
            "--json",
            self.dir / "keno.json",
            "keno",
            "5A",
            self.dir / "spec.txt",
            self.dir / "core.out",
            self.dir / "core-pins.msgpack",
        )
        self.assertEqual((run.returncode, run.stdout), (0, ""), run.stderr)
        s = json.loads((self.dir / "keno.json").read_text())
        self.assertEqual((s["kind"], s["problem"]), ("keno", "5A"))
        self.assertAlmostEqual(s["k"]["diff_pcm"], 100.0, places=2)
        self.assertEqual(s["k"]["reference_std"], 8e-6)
        self.assertEqual(s["assemblies"]["relative"]["at"], [i, j])
        self.assertEqual(
            s["assemblies"]["relative"]["name"], f"{'RPNMLKJHGFEDCBA'[j]}{i + 1:02d}"
        )
        self.assertEqual(s["max_pin"]["reference"], 1.5)


if __name__ == "__main__":
    unittest.main()
