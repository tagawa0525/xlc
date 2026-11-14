"""transfer_check.py のログの読み方と反応率の差の和（#0285、ADR-0250）。

実行: cd examples/depletion/issue-0285 && python3 -m unittest test_transfer_check
"""

import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_check import blocks


def line(
    cell: int,
    ring: int,
    group: int,
    weight: float,
    state: float,
    elastic: float,
    static: str,
) -> str:
    return (
        "t DEBUG ratio_check: ratio check of the ring transfer "
        f'cell={cell} ring={ring} nuclide="H-1" group={group} weight={weight} '
        f"state_out={state} state_up=0.1 state_below=0.0 "
        f"elastic_out={elastic} elastic_up=0.1 elastic_below=0.0 "
        f"static_out={static} static_up={static} static_below={static}"
    )


CLOSE = (
    " INFO ratio-shielding: xlc_lattice::shielding: close time.busy=1s time.idle=1µs"
)


class Blocks(unittest.TestCase):
    def test_differences_are_weighted_and_summed_over_pins_and_split_by_reference(self):
        text = "\n".join(
            [
                line(0, 7, 20, 2.0, 1.0, 1.5, "1.25"),
                line(1, 7, 20, 1.0, 2.0, 1.0, "2.0"),
                line(0, 7, 21, 4.0, 1.0, 1.0, "1.0"),
                CLOSE,
                line(0, 7, 20, 1.0, 1.0, 0.5, "1.0"),
                CLOSE,
            ]
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "a.log"
            path.write_text(text + "\n")
            found = blocks(path)
        self.assertEqual(len(found), 2)
        group = found[0]["groups"][(7, "H-1", 20)]
        # 2 × (1.5 − 1) + 1 × (1 − 2)
        self.assertAlmostEqual(group["elastic_out"], 0.0)
        # 2 × (1.25 − 1) + 1 × (2 − 2)
        self.assertAlmostEqual(group["static_out"], 0.5)
        self.assertAlmostEqual(group["state_out"], 4.0)
        ring = found[0]["rings"][(7, "H-1")]
        self.assertAlmostEqual(ring["state_out"], 8.0)
        self.assertAlmostEqual(ring["static_out"], 0.5)
        self.assertAlmostEqual(found[1]["groups"][(7, "H-1", 20)]["elastic_out"], -0.5)

    def test_rows_without_static_counts_give_nan(self):
        # 静止した標的の行を数えなかった実行では static の値は NaN で、和も NaN（0 と取り違えない）
        text = "\n".join([line(0, 7, 20, 1.0, 1.0, 1.0, "NaN"), CLOSE])
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "a.log"
            path.write_text(text + "\n")
            found = blocks(path)
        self.assertTrue(math.isnan(found[0]["groups"][(7, "H-1", 20)]["static_out"]))


if __name__ == "__main__":
    unittest.main()
