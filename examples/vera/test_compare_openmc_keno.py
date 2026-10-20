"""compare_openmc_keno.py の k の差の計算（xlc − OpenMC、OpenMC − KENO-VI、統計誤差の合成）。

実行: python3 -m unittest discover -s examples/vera -p 'test_compare_openmc_keno.py'
"""

import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = HERE.parents[1] / "references/vera/spec.txt"


@unittest.skipUnless(SPEC.exists(), "references/vera/spec.txt がない")
class Differences(unittest.TestCase):
    def test_k_differences_and_errors(self):
        """1B の KENO-VI は 1.182149 ± 0.000068。xlc 1.184, OpenMC 1.1830 ± 0.00004 なら xlc − OpenMC は +100 ± 4 pcm、OpenMC − KENO-VI は +85.1 pcm で誤差は二乗和の平方根。"""
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "xlc.json").write_text(
                json.dumps({"k": {"xlc": 1.184}, "problem": "1B"})
            )
            (d / "results").mkdir()
            (d / "results/k.json").write_text(
                json.dumps({"k": 1.183, "k-std": 4e-5, "problem": "1B"})
            )
            out = d / "compare.json"
            run = subprocess.run(
                [
                    sys.executable,
                    str(HERE / "compare_openmc_keno.py"),
                    "--json",
                    str(out),
                    str(SPEC),
                    "1B",
                    str(d / "xlc.json"),
                    str(d / "results"),
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(run.stdout, "")
            k = json.loads(out.read_text())["k"]
        self.assertAlmostEqual(k["diff_xlc_openmc_pcm"], 100.0, places=6)
        self.assertAlmostEqual(k["diff_xlc_openmc_std_pcm"], 4.0, places=6)
        self.assertAlmostEqual(k["diff_openmc_keno_pcm"], 85.1, places=6)
        self.assertAlmostEqual(
            k["diff_openmc_keno_std_pcm"], math.hypot(4.0, 6.8), places=6
        )


if __name__ == "__main__":
    unittest.main()
