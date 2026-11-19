"""openmc_depletion.py の、ピンセルの燃焼の入力の読み方（ADR-0173）。

OpenMC のある環境で実行する（入力を読むだけで、輸送は解かない）:
nix develop .#references -c python3 -m unittest discover -s examples/depletion -p 'test_openmc_depletion.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import openmc_depletion

INPUT = """
problem = "1C"
rings = 1
cross-sections = "openmc-library/cross_sections.xml"
chain = "openmc-chain.xml"
power = 38.41
burnup = [0.1, 0.5]
particles = 1000
batches = 10
inactive = 5
seed = 1
threads = 1
normalization = "fission-q"
"""


def read_text(text: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        source = Path(tmp) / "openmc-1C.toml"
        source.write_text(text)
        return openmc_depletion.read(source)


class Normalization(unittest.TestCase):
    def test_the_chain_fission_q_reads(self):
        self.assertEqual(read_text(INPUT)["normalization"], "fission-q")

    def test_energy_deposition_can_be_chosen(self):
        run = read_text(INPUT.replace('"fission-q"', '"energy-deposition"'))
        self.assertEqual(run["normalization"], "energy-deposition")

    def test_an_input_without_normalization_is_an_error(self):
        with self.assertRaises(SystemExit):
            read_text(INPUT.replace('normalization = "fission-q"\n', ""))

    def test_an_unknown_normalization_is_an_error(self):
        with self.assertRaises(SystemExit):
            read_text(INPUT.replace('"fission-q"', '"source-rate"'))


if __name__ == "__main__":
    unittest.main()
