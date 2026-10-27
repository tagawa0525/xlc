"""compare-dvc-lock.py の検査（lock の比較は文字列を受け、git に依らない）。

実行: python3 -m unittest discover -s scripts -p 'test_compare_dvc_lock.py'
"""

import importlib.util
import textwrap
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "compare_dvc_lock", Path(__file__).resolve().parent / "compare-dvc-lock.py"
)
cmp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cmp)

LOCK = """\
schema: '2.0'
stages:
  a:
    cmd: run a
    deps:
    - path: in.txt
      hash: md5
      md5: {dep}
      size: 1
    params:
      p.toml:
        x: 1
    outs:
    - path: out.txt
      hash: md5
      md5: {out}
      size: 2
"""


def lock(dep="d1", out="o1"):
    return textwrap.dedent(LOCK).format(dep=dep, out=out)


class CompareLocks(unittest.TestCase):
    def test_same_stages_with_changed_deps_hash_has_no_difference(self):
        self.assertEqual(cmp.compare_locks("f", lock(dep="d1"), lock(dep="d2")), [])

    def test_changed_cmd_is_reported_with_stage_and_field(self):
        head = lock().replace("run a", "run b")
        self.assertEqual(
            cmp.compare_locks("f", lock(), head), ["f: stage a: cmd differs"]
        )

    def test_changed_outs_hash_is_reported(self):
        self.assertEqual(
            cmp.compare_locks("f", lock(out="o1"), lock(out="o2")),
            ["f: stage a: outs differs"],
        )

    def test_changed_params_is_reported(self):
        head = lock().replace("x: 1", "x: 2")
        self.assertEqual(
            cmp.compare_locks("f", lock(), head), ["f: stage a: params differs"]
        )

    def test_added_and_removed_stages_are_reported(self):
        head = lock() + "  b:\n    cmd: run b\n"
        self.assertEqual(cmp.compare_locks("f", lock(), head), ["f: stage b: added"])
        self.assertEqual(cmp.compare_locks("f", head, lock()), ["f: stage b: removed"])

    def test_missing_lock_is_an_empty_stage_set(self):
        self.assertEqual(cmp.compare_locks("f", None, lock()), ["f: stage a: added"])

    def test_stage_without_params_equals_stage_without_params(self):
        text = "schema: '2.0'\nstages:\n  s:\n    cmd: c\n"
        self.assertEqual(cmp.compare_locks("f", text, text), [])


if __name__ == "__main__":
    unittest.main()
