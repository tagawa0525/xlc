"""compare_keno.compare と compare_keno_lattice の KENO-VI の読み方。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_keno_lattice.py'
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_keno
import compare_p3

SPEC = Path(__file__).resolve().parents[2] / "references/vera/spec.txt"


@unittest.skipUnless(SPEC.exists(), "references/vera/spec.txt がない")
class KenoValues(unittest.TestCase):
    def setUp(self):
        self.spec = SPEC.read_text().split("\n")

    def test_k_of_2b_is_the_endf_b_vii_0_table_value(self):
        """仕様書の ENDF/B-VII.0 の表（付録）の 2B は 1.183360 ± 0.000024。"""
        self.assertEqual(compare_p3.keno_k(self.spec, "2B"), (1.183360, 0.000024))

    def test_k_of_1b_is_the_endf_b_vii_0_table_value(self):
        """Problem 1 は ENDF/B-VII.0 の表が ENDF/B-VI.8 の表より先にある。1B は 1.182149 ± 0.000068。"""
        self.assertEqual(compare_p3.keno_k(self.spec, "1B"), (1.182149, 0.000068))

    def test_pins_equal_to_the_keno_table_have_no_difference(self):
        """KENO-VI の 1/8 の表を 8 本の対称な位置に広げた分布は、差が 0。"""
        octant = compare_keno.keno_octant(self.spec, compare_keno.CAPTIONS["2B"])
        rows = [[0.0] * 17 for _ in range(17)]
        for (i, j), value in octant.items():
            for a, b in [
                (i, j),
                (j, i),
                (-i, j),
                (j, -i),
                (i, -j),
                (-j, i),
                (-i, -j),
                (-j, -i),
            ]:
                rows[8 + a][8 + b] = value
        n, rms, worst, _ = compare_keno.compare(self.spec, "2B", rows)
        self.assertGreater(n, 30)
        self.assertLess(
            rms, 1e-3
        )  # 表の値の平均で 1 に規格化し直すので、表の平均が 1 からずれた分だけ残る
        self.assertLess(abs(worst), 1e-3)

    def test_json_has_k_difference_and_no_stdout(self):
        """k 1.18 の 2B は、KENO-VI の 1.183360 に対して −336 pcm。要約は JSON だけに書き、標準出力には何も出さない。"""
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "lattice.out"
            out.write_text("k 1.18\n")
            summary = Path(d) / "compare.json"
            run = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve().parent / "compare_keno_lattice.py"),
                    "--json",
                    str(summary),
                    str(SPEC),
                    "2B",
                    str(out),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(run.stdout, "")
            data = json.loads(summary.read_text())
        self.assertEqual(data["format"], "xlc-compare-keno-lattice")
        self.assertAlmostEqual(data["k"]["diff_pcm"], -336.0, places=6)
        self.assertIsNone(data["pins"])


if __name__ == "__main__":
    unittest.main()
