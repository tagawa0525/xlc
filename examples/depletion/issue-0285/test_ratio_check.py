"""ratio_check.py のログの読み方と重みの平均（#0285、ADR-0250）。

実行: cd examples/depletion/issue-0285 && python3 -m unittest test_ratio_check
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ratio_check import blocks


def line(ring: int, group: int, equivalence: float, rate: float) -> str:
    return (
        "t DEBUG ratio_check: ratio check of the ring capture cell=0 "
        f'ring={ring} nuclide="U-238" group={group} equivalence={equivalence} state=0.98 '
        f"mean=0.97 dancoff=0.99 rate={rate}"
    )


class Blocks(unittest.TestCase):
    def test_ratios_are_rate_weighted_per_ring_and_group_and_split_by_reference(self):
        close = " INFO ratio-shielding: xlc_lattice::shielding: close time.busy=1s time.idle=1µs"
        text = "\n".join(
            [
                line(0, 13, 0.99, 1.0),
                line(0, 14, 0.96, 3.0),
                line(1, 13, 0.98, 2.0),
                " INFO ratio-shielding:lumps: xlc_lattice::shielding: close time.busy=1s",
                close,
                line(0, 13, 1.01, 1.0),
                close,
            ]
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "a.log"
            path.write_text(text + "\n")
            found = blocks(path)
        self.assertEqual(len(found), 2)
        first = found[0]
        # リング 0: (0.99 × 1 + 0.96 × 3) / 4 − 1
        self.assertAlmostEqual(first["rings"][(0, "U-238")]["equivalence"], -0.0325)
        # 群 13: (0.99 × 1 + 0.98 × 2) / 3 − 1
        self.assertAlmostEqual(first["groups"][13]["equivalence"], 2.95 / 3 - 1)
        self.assertAlmostEqual(first["rings"][(1, "U-238")]["state"], -0.02)
        self.assertAlmostEqual(found[1]["rings"][(0, "U-238")]["equivalence"], 0.01)


if __name__ == "__main__":
    unittest.main()
