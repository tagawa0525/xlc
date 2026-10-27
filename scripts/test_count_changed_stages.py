"""count-changed-stages.py の検査（数えは dvc status --json の dict を受け、dvc に依らない）。

実行: python3 -m unittest discover -s scripts -p 'test_count_changed_stages.py'
"""

import importlib.util
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "count_changed_stages", Path(__file__).resolve().parent / "count-changed-stages.py"
)
cnt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cnt)

STATUS = {
    "examples/a/dvc.yaml:s1": [{"changed deps": {"x": "modified"}}],
    "examples/a/dvc.yaml:s2": [
        {"changed deps": {"x": "modified"}},
        {"changed outs": {"o": "modified"}},
    ],
    "examples/a/dvc.yaml:s3": [{"changed outs": {"o": "deleted"}}],
    "examples/b/dvc.yaml:s1": [{"changed deps": {"y": "modified"}}],
    "build": [{"changed deps": {"src": "modified"}}],
}


class Count(unittest.TestCase):
    def test_counts_stages_per_yaml_and_kind(self):
        self.assertEqual(
            cnt.count_changed(STATUS),
            {
                "examples/a/dvc.yaml": {"changed deps": 2, "changed outs": 2},
                "examples/b/dvc.yaml": {"changed deps": 1, "changed outs": 0},
                "dvc.yaml": {"changed deps": 1, "changed outs": 0},
            },
        )

    def test_empty_status_has_no_yaml(self):
        self.assertEqual(cnt.count_changed({}), {})

    def test_unknown_change_kind_is_an_error(self):
        with self.assertRaises(ValueError):
            cnt.count_changed({"dvc.yaml:s": [{"not in cache": {"o": "x"}}]})


if __name__ == "__main__":
    unittest.main()
