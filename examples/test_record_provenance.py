"""record_provenance.py の検査（dvc.lock の読み方、md5 を入れたコミット、JSON の inputs との照合）。

実行: python3 -m unittest discover -s examples -p 'test_record_provenance.py'
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import record_provenance as rp

LOCK = """schema: '2.0'
stages:
  a:
    cmd: echo a
    deps:
    - path: in.txt
      hash: md5
      md5: {a_in}
      size: 1
    outs:
    - path: a.out
      hash: md5
      md5: {a_out}
      size: 2
  b:
    cmd: echo b
    deps:
    - path: a.out
      hash: md5
      md5: {a_out}
      size: 2
    outs:
    - path: b.json
      hash: md5
      md5: {b_out}
      size: 3
    - path: dir
      hash: md5
      md5: {d}.dir
      size: 4
      nfiles: 1
"""
H = {k: k * 32 for k in "0123456789abcdef"}


def lock(a_in="1", a_out="2", b_out="3", d="4") -> str:
    return LOCK.format(a_in=H[a_in], a_out=H[a_out], b_out=H[b_out], d=H[d])


class Parse(unittest.TestCase):
    def test_reads_deps_and_outs_of_one_stage(self):
        stage = rp.parse_stage(lock(), "b")
        self.assertEqual(stage["deps"], [("a.out", H["2"])])
        self.assertEqual(stage["outs"], [("b.json", H["3"]), ("dir", H["4"] + ".dir")])
        self.assertEqual(rp.parse_stage(lock(), "a")["deps"], [("in.txt", H["1"])])

    def test_unknown_stage_is_an_error(self):
        with self.assertRaises(rp.ProvenanceError):
            rp.parse_stage(lock(), "c")


class Git(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.git("init", "-q", "-b", "main")

    def git(self, *args: str) -> str:
        env = {
            **os.environ,
            "GIT_AUTHOR_NAME": "t",
            "GIT_AUTHOR_EMAIL": "t@example.com",
            "GIT_COMMITTER_NAME": "t",
            "GIT_COMMITTER_EMAIL": "t@example.com",
        }
        return subprocess.run(
            ["git", "-C", str(self.dir), *args],
            check=True,
            capture_output=True,
            text=True,
            env=env,
        ).stdout.strip()

    def commit(self, text: str, date: str) -> str:
        (self.dir / "dvc.lock").write_text(text)
        self.git("add", "dvc.lock")
        env_date = f"{date}T12:00:00"
        subprocess.run(
            ["git", "-C", str(self.dir), "commit", "-q", "-m", "m"],
            check=True,
            capture_output=True,
            env={
                **os.environ,
                "GIT_AUTHOR_NAME": "t",
                "GIT_AUTHOR_EMAIL": "t@example.com",
                "GIT_COMMITTER_NAME": "t",
                "GIT_COMMITTER_EMAIL": "t@example.com",
                "GIT_AUTHOR_DATE": env_date,
                "GIT_COMMITTER_DATE": env_date,
            },
        )
        return self.git("rev-parse", "--short=8", "HEAD")

    def test_names_the_commit_that_first_locked_the_current_md5(self):
        first = self.commit(lock(), "2026-10-01")
        # b だけが再実行した relock（a の出力は同じ md5）と、a の deps の変更なしの再書き込み
        self.commit(lock(b_out="5"), "2026-10-02")
        third = self.commit(lock(b_out="5"), "2026-10-03")
        lock_path = self.dir / "dvc.lock"
        a = rp.locked_since(lock_path, "a", "a.out")
        self.assertEqual(a, (first, "2026-10-01"))
        b = rp.locked_since(lock_path, "b", "b.json")
        self.assertEqual(b[1], "2026-10-02")
        self.assertNotEqual(b[0], third)

    def test_a_md5_that_returns_to_an_old_value_counts_from_its_return(self):
        self.commit(lock(), "2026-10-01")
        self.commit(lock(a_out="6"), "2026-10-02")
        back = self.commit(lock(), "2026-10-03")
        self.assertEqual(
            rp.locked_since(self.dir / "dvc.lock", "a", "a.out"),
            (back, "2026-10-03"),
        )

    def test_uncommitted_and_untracked_are_said_not_hidden(self):
        self.commit(lock(), "2026-10-01")
        (self.dir / "dvc.lock").write_text(lock(a_out="7"))
        self.assertEqual(
            rp.locked_since(self.dir / "dvc.lock", "a", "a.out")[0], "uncommitted"
        )
        lone = self.dir / "sub"
        lone.mkdir()
        (lone / "dvc.lock").write_text(lock())
        self.assertEqual(rp.locked_since(lone / "dvc.lock", "a", "a.out")[0], "unknown")


class Inputs(unittest.TestCase):
    def test_json_inputs_must_match_the_locked_deps(self):
        stage = rp.parse_stage(lock(), "b")
        ok = {"inputs": [{"path": "a.out", "md5": H["2"]}]}
        self.assertEqual(rp.input_problems(ok, stage), [])
        bad = {"inputs": [{"path": "a.out", "md5": H["9"]}]}
        self.assertEqual(len(rp.input_problems(bad, stage)), 1)
        missing = {"inputs": [{"path": "zzz", "md5": H["9"]}]}
        self.assertEqual(len(rp.input_problems(missing, stage)), 1)


if __name__ == "__main__":
    unittest.main()
