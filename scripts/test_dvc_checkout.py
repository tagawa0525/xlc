"""dvc-checkout.py の検査（分けるのは dvc checkout の出力の文字列を受け、dvc に依らない）。

実行: python3 -m unittest discover -s scripts -p 'test_dvc_checkout.py'
"""

import importlib.util
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "dvc_checkout", Path(__file__).resolve().parent / "dvc-checkout.py"
)
chk = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(chk)

ROOT = Path("/repo")


def output(no_hash, failed):
    """dvc 3.67 の dvc checkout の出力の形（警告は絶対パス、失敗の一覧はリポジトリの直下からのパス）"""
    lines = [
        f"WARNING: No file hash info found for '{ROOT / p}'. It won't be created."
        for p in no_hash
    ]
    lines.append("A       examples/a/out.msgpack")
    if failed:
        lines.append("ERROR: Checkout failed for following targets:")
        lines.extend(failed)
        lines.append("Is your cache up to date?")
        lines.append("<https://error.dvc.org/missing-files>")
    return "\n".join(lines) + "\n"


class Classify(unittest.TestCase):
    def test_outputs_without_a_hash_are_unrun_not_missing(self):
        out = output(
            ["examples/d/openmc-mox1", "examples/v/openmc-1B"],
            ["examples/v/openmc-1B", "examples/d/openmc-mox1"],
        )
        self.assertEqual(
            chk.classify(out, ROOT),
            (["examples/d/openmc-mox1", "examples/v/openmc-1B"], []),
        )

    def test_an_output_with_a_hash_that_did_not_restore_is_missing(self):
        out = output(
            ["examples/d/openmc-mox1"],
            ["examples/d/openmc-mox1", "examples/a/lattice.msgpack"],
        )
        self.assertEqual(
            chk.classify(out, ROOT),
            (["examples/d/openmc-mox1"], ["examples/a/lattice.msgpack"]),
        )

    def test_a_failure_without_the_list_of_targets_is_an_error(self):
        with self.assertRaises(ValueError):
            chk.classify("ERROR: unexpected error - lock is held\n", ROOT)

    def test_a_warning_outside_the_root_is_an_error(self):
        out = (
            "WARNING: No file hash info found for '/elsewhere/x'. It won't be created.\n"
            "ERROR: Checkout failed for following targets:\nx\nIs your cache up to date?\n"
        )
        with self.assertRaises(ValueError):
            chk.classify(out, ROOT)


if __name__ == "__main__":
    unittest.main()
