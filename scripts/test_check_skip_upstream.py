"""check-skip-upstream.py の検査（判定は dvc status --with-deps --json の dict を受け、dvc に依らない）。

実行: python3 -m unittest discover -s scripts -p 'test_check_skip_upstream.py'
"""

import importlib.util
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "check_skip_upstream", Path(__file__).resolve().parent / "check-skip-upstream.py"
)
skip = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(skip)

TARGET = "examples/vera/dvc.yaml:openmc-1B-model"


class Blockers(unittest.TestCase):
    def test_upstream_with_only_changed_deps_does_not_block(self):
        status = {
            "build": [{"changed deps": {"crates": "modified"}}],
            "examples/vera/dvc.yaml:xs": [{"changed deps": {"bin/xs": "modified"}}],
        }
        self.assertEqual(skip.blockers(status, TARGET), [])

    def test_the_target_itself_is_not_counted(self):
        status = {
            TARGET: [
                {"changed deps": {"x": "modified"}},
                {"changed outs": {"o": "deleted"}},
                "changed command",
            ]
        }
        self.assertEqual(skip.blockers(status, TARGET), [])

    def test_always_changed_item_does_not_block_but_its_changed_outs_does(self):
        status = {
            "examples/vera/dvc.yaml:references-shell": ["always changed"],
        }
        self.assertEqual(skip.blockers(status, TARGET), [])
        status["examples/vera/dvc.yaml:references-shell"].append(
            {"changed outs": {"references-shell.txt": "modified"}}
        )
        self.assertEqual(
            skip.blockers(status, TARGET),
            [("examples/vera/dvc.yaml:references-shell", ["changed outs"])],
        )

    def test_upstream_with_other_kinds_blocks(self):
        status = {
            "examples/vera/dvc.yaml:xs": [
                {"changed deps": {"bin/xs": "modified"}},
                {"changed outs": {"vera.msgpack": "modified"}},
            ],
            "examples/vera/dvc.yaml:xs-range": ["changed command", "changed checksum"],
        }
        self.assertEqual(
            skip.blockers(status, TARGET),
            [
                ("examples/vera/dvc.yaml:xs", ["changed outs"]),
                (
                    "examples/vera/dvc.yaml:xs-range",
                    ["changed checksum", "changed command"],
                ),
            ],
        )

    def test_root_stage_target_matches_the_key_without_the_yaml(self):
        status = {"build": [{"changed outs": {"bin/xs": "modified"}}]}
        self.assertEqual(skip.blockers(status, "dvc.yaml:build"), [])

    def test_unknown_kind_is_an_error(self):
        status = {"examples/vera/dvc.yaml:xs": [{"not in cache": {}}]}
        with self.assertRaises(ValueError):
            skip.blockers(status, TARGET)


if __name__ == "__main__":
    unittest.main()
