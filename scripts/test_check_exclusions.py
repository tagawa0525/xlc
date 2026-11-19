"""check-exclusions.py の検査（ファイルの中身と issue の状態を dict で受け、git に依らない）。

実行: python3 -m unittest discover -s scripts -p 'test_check_exclusions.py'
"""

import importlib.util
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "check_exclusions", Path(__file__).resolve().parent / "check-exclusions.py"
)
exc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exc)

STATUSES = {"0001": "open", "0002": "closed"}

# 検査は行頭の skip のデコレータを拾うので、このファイルが当たらないように組み立てる
SKIP = "@" + "unittest.skip"
IGNORE = "#" + "[ignore"


def problems(files):
    return exc.violations(exc.find_markers(files), STATUSES)


class TextlintIgnore(unittest.TestCase):
    def test_entries_under_a_comment_with_an_open_issue_pass(self):
        files = {".textlintignore": "# 後で直す（#0001）\na.toml\nb.toml\n"}
        self.assertEqual(problems(files), [])

    def test_entries_without_an_issue_are_reported_once_per_entry(self):
        files = {".textlintignore": "# 理由だけ\na.toml\nb.toml\n"}
        found = problems(files)
        self.assertEqual(len(found), 2)
        self.assertIn(".textlintignore:2", found[0])

    def test_a_blank_line_ends_the_comment_of_a_block(self):
        files = {".textlintignore": "# 後で直す（#0001）\na.toml\n\nb.toml\n"}
        found = problems(files)
        self.assertEqual(len(found), 1)
        self.assertIn(".textlintignore:4", found[0])

    def test_entries_whose_only_issue_is_closed_are_reported(self):
        files = {".textlintignore": "# 直るまで（#0002）\na.toml\n"}
        found = problems(files)
        self.assertEqual(len(found), 1)
        self.assertIn("#0002", found[0])


class SkippedTests(unittest.TestCase):
    def test_rust_ignore_needs_an_open_issue_or_an_adr(self):
        files = {
            "crates/a/src/tests.rs": "\n".join(
                [
                    IGNORE + ' = "slow (#0001)"]',
                    IGNORE + ' = "needs a DVC output (ADR-0130)"]',
                    IGNORE + ' = "survey of #0002 / #0001"]',
                    IGNORE + "]",
                    IGNORE + ' = "until fixed (#0002)"]',
                ]
            )
        }
        found = problems(files)
        self.assertEqual(len(found), 2)
        self.assertIn("crates/a/src/tests.rs:4", found[0])
        self.assertIn("crates/a/src/tests.rs:5", found[1])

    def test_unconditional_python_skips_need_a_reference(self):
        files = {
            "examples/test_a.py": "\n".join(
                [
                    SKIP + '("後で")',
                    SKIP + '("後で（#0001）")',
                    SKIP + 'Unless(SPEC.exists(), "参照がない")',
                ]
            )
        }
        found = problems(files)
        self.assertEqual(len(found), 1)
        self.assertIn("examples/test_a.py:1", found[0])


class InputComments(unittest.TestCase):
    def test_an_unverified_premise_needs_an_open_issue_in_its_comment_block(self):
        files = {
            "examples/a/xs.toml": "# 核種の組は文献の一覧（本文は\n# 確かめていない。ADR-0169）\nx = 1\n",
            "examples/b/xs.toml": "# 仮の値（#0001 で決める）\nx = 1\n",
            "examples/c/xs.toml": "# 暫定（#0002）\nx = 1\n",
        }
        found = problems(files)
        self.assertEqual(len(found), 2)
        self.assertIn("examples/a/xs.toml:2", found[0])
        self.assertIn("examples/c/xs.toml:1", found[1])

    def test_an_exclusion_that_cites_only_closed_issues_is_reported(self):
        files = {
            "examples/a/xs.toml": "# 表を作る（#0001）\n# Pm-147 は NaN になるので外す（#0002）\nx = 1\n",
            "examples/b/xs.toml": "# 不連続因子を外した定数\n# 時間のために外した（#0001、#0002 の後）\nx = 1\n",
        }
        found = problems(files)
        self.assertEqual(len(found), 1)
        self.assertIn("examples/a/xs.toml:2", found[0])
        self.assertIn("#0002", found[0])

    def test_trial_inputs_of_an_issue_are_not_checked(self):
        files = {
            "examples/depletion/issue-0002/pin.toml": "# 時間のために外した（#0002）\n",
            "examples/depletion/trials-x/pin.toml": "# 仮の値\n",
        }
        self.assertEqual(problems(files), [])

    def test_values_are_not_comments(self):
        files = {"examples/a/xs.toml": 'note = "仮の値"\n'}
        self.assertEqual(problems(files), [])


class Statuses(unittest.TestCase):
    def test_an_unknown_issue_is_an_error(self):
        files = {".textlintignore": "# 後で（#0999）\na.toml\n"}
        with self.assertRaises(ValueError):
            problems(files)

    def test_listing_counts_markers_by_kind(self):
        files = {
            ".textlintignore": "# 後で直す（#0001）\na.toml\nb.toml\n",
            "examples/b/xs.toml": "# 仮の値（#0001 で決める）\nx = 1\n",
        }
        counts = exc.count_by_kind(exc.find_markers(files))
        self.assertEqual(counts, {"textlintignore": 2, "provisional": 1})

    def test_listing_counts_open_issues_and_review_leftovers_by_name(self):
        names = {
            "0001": "0001-frame-gap.md",
            "0002": "0002-core-review-leftovers.md",
            "0003": "0003-adr-review-points.md",
            "0004": "0004-rod-branch-leftovers.md",
        }
        statuses = {"0001": "open", "0002": "open", "0003": "closed", "0004": "open"}
        self.assertEqual(
            exc.count_open_issues(statuses, names),
            {"open-issues": 3, "open-review-leftover-issues": 2},
        )


if __name__ == "__main__":
    unittest.main()
