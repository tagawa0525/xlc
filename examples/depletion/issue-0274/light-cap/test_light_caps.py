"""light_caps.py の比較（分離共鳴を持たない核種の自由ガスの上限を変えたピンセル。ADR-0246）。

実行: cd examples/depletion/issue-0274/light-cap && python3 -m unittest test_light_caps
"""

import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import light_caps as lc

SCATTER = [[0.4, 0.02], [0.0, 1.2]]


def write(directory: Path, name: str, k: float, absorption: list[float]) -> str:
    path = directory / name
    table = {"k_inf": k, "absorption": absorption, "scatter": SCATTER}
    path.write_bytes(msgpack.packb({"states": [table]}))
    return str(path)


class LightCapsTest(unittest.TestCase):
    def test_reactivity_and_absorption_are_relative_to_the_reference(self):
        rho, thermal, absorption = lc.compare(
            {"k_inf": 1.25, "absorption": [0.01, 0.1], "scatter": SCATTER},
            {
                "k_inf": 1.2502,
                "absorption": [0.010001, 0.1],
                "scatter": [[0.4, 0.0202], [0.0, 1.3]],
            },
        )
        self.assertAlmostEqual(rho, (1 / 1.2502 - 1 / 1.25) * 1e5)
        # 熱群への移行は熱群より上の群からの和で、熱群の中の散乱を含まない
        self.assertAlmostEqual(thermal, 1000.0, places=6)
        self.assertAlmostEqual(absorption[0], 10.0, places=6)
        self.assertEqual(absorption[1], 0.0)

    def test_the_multiple_of_kt_is_read_from_the_name(self):
        self.assertEqual(lc.multiple("a/l-350.msgpack"), 350.0)
        self.assertEqual(lc.multiple("l-dilute-400.msgpack"), 400.0)
        self.assertEqual(lc.multiple("l-inf.msgpack"), float("inf"))
        self.assertIsNone(lc.multiple("pin-1C-0-light-inf.msgpack"))

    def test_rows_are_sorted_by_the_multiple_with_the_cap_energies(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            reference = write(d, "l-inf.msgpack", 1.25, [0.01])
            paths = [
                write(d, "l-400.msgpack", 1.2503, [0.0099]),
                write(d, "l-300.msgpack", 1.2502, [0.0101]),
                write(d, "other.msgpack", 1.25, [0.01]),
            ]
            rows = lc.rows(reference, paths)
        self.assertEqual(
            [r.split()[0] for r in rows],
            ["l-300.msgpack", "l-400.msgpack", "other.msgpack"],
        )
        # 300 kT は 900 K で 23.27 eV、580 K で 14.99 eV
        self.assertEqual(rows[0].split()[1:3], ["23.27", "14.99"])
        self.assertEqual(rows[2].split()[1:3], ["+0.00", "+0.00"])


if __name__ == "__main__":
    unittest.main()
