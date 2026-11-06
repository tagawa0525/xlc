"""record_provenance.py の検査（dvc.lock と .dvc の読み方、md5 を入れたコミット、JSON の inputs との照合）。

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

    def test_reads_a_path_that_dvc_folded_onto_the_next_line(self):
        """dvc は長いパスを `- path:` の次の行に折る。折られた dep を読み落とすと、JSON の inputs が dep でないと誤る"""
        lock = (
            "schema: '2.0'\nstages:\n  a:\n    cmd: echo a\n    deps:\n"
            "    - path: \n        ../very/long/path/to/the/file.msgpack\n"
            "      hash: md5\n      md5: " + "1" * 32 + "\n      size: 1\n"
        )
        self.assertEqual(
            rp.parse_stage(lock, "a")["deps"],
            [("../very/long/path/to/the/file.msgpack", "1" * 32)],
        )

    def test_a_folded_path_without_its_continuation_is_an_error(self):
        """続きの行が無いのに md5 を読み進めると、その dep を黙って捨てる"""
        lock = (
            "schema: '2.0'\nstages:\n  a:\n    cmd: echo a\n    deps:\n"
            "    - path: \n      hash: md5\n      md5: "
            + "1" * 32
            + "\n      size: 1\n"
        )
        with self.assertRaises(rp.ProvenanceError):
            rp.parse_stage(lock, "a")

    def test_a_folded_path_may_be_indented_by_any_amount(self):
        lock = (
            "schema: '2.0'\nstages:\n  a:\n    cmd: echo a\n    deps:\n"
            "    - path: \n          ../x/y.msgpack\n"
            "      hash: md5\n      md5: " + "2" * 32 + "\n      size: 1\n"
        )
        self.assertEqual(
            rp.parse_stage(lock, "a")["deps"], [("../x/y.msgpack", "2" * 32)]
        )

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

    def commit(self, text: str, date: str, name: str = "dvc.lock") -> str:
        (self.dir / name).write_text(text)
        self.git("add", name)
        env_date = f"{date}T12:00:00"
        subprocess.run(
            [
                "git",
                "-c",
                "core.hooksPath=/dev/null",  # 利用者の commit-msg hook を持ち込まない
                "-C",
                str(self.dir),
                "commit",
                "-q",
                "-m",
                "chore: m",
            ],
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
        # b だけが再実行した relock（a の出力は同じ md5）と、a の入力だけが変わった relock（出力の md5 は同じ）
        self.commit(lock(b_out="5"), "2026-10-02")
        third = self.commit(lock(a_in="8", b_out="5"), "2026-10-03")
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


class Show(Git):
    """JSON を持つステージの出力の検査（JSON が無い、md5 が違う、inputs が無い）"""

    def stage_with_json(self, body: str | None, md5: str) -> None:
        (self.dir / "dvc.yaml").write_text("stages: {}\n")
        lock = (
            "schema: '2.0'\nstages:\n  c:\n    cmd: x\n    deps:\n    - path: in\n"
            f"      hash: md5\n      md5: {H['1']}\n      size: 1\n    outs:\n"
            f"    - path: c.json\n      hash: md5\n      md5: {md5}\n      size: 1\n"
        )
        self.commit(lock, "2026-10-01")
        if body is not None:
            (self.dir / "c.json").write_text(body)

    def run_show(self) -> tuple[int, str]:
        import contextlib
        import io

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = rp.show(self.dir / "dvc.yaml", "c")
        return status, out.getvalue()

    def test_missing_json_is_reported(self):
        self.stage_with_json(None, H["3"])
        status, out = self.run_show()
        self.assertEqual(status, 1)
        self.assertIn("MISSING", out)

    def test_json_differing_from_the_lock_is_reported(self):
        self.stage_with_json("{}", H["3"])
        status, out = self.run_show()
        self.assertEqual(status, 1)
        self.assertIn("MISMATCH", out)

    def test_json_without_inputs_is_reported(self):
        import hashlib

        body = '{"format": "xlc-compare-p3"}'
        self.stage_with_json(body, hashlib.md5(body.encode()).hexdigest())
        status, out = self.run_show()
        self.assertEqual(status, 1)
        self.assertIn("has no inputs", out)

    def test_a_changed_dep_makes_the_stage_uncommitted(self):
        self.stage_with_json(None, H["3"])
        lock = self.dir / "dvc.lock"
        lock.write_text(lock.read_text().replace(H["1"], H["8"]))
        self.assertEqual(rp.locked_since(lock, "c", "c.json")[0], "uncommitted")


DVC = """outs:
- md5: {md5}
  size: {size}
  hash: md5
  path: run.log
"""


def dvc(md5: str = "1", size: int = 10) -> str:
    return DVC.format(md5=H[md5], size=size)


class DvcFile(Git):
    """`dvc add` で追う生成物の .dvc（ADR-0222）"""

    def test_reads_outs_in_any_key_order(self):
        text = (
            "outs:\n- path: a.log\n  md5: " + H["1"] + "\n  size: 1\n"
            "- hash: md5\n  md5: " + H["2"] + ".dir\n  nfiles: 2\n  path: d\n"
        )
        self.assertEqual(
            rp.parse_dvc(text), [("a.log", H["1"]), ("d", H["2"] + ".dir")]
        )

    def test_an_out_without_md5_is_an_error(self):
        with self.assertRaises(rp.ProvenanceError):
            rp.parse_dvc("outs:\n- path: a.log\n  size: 1\n")

    def test_a_file_without_outs_is_an_error(self):
        with self.assertRaises(rp.ProvenanceError):
            rp.parse_dvc("deps: []\n")

    def test_lists_each_md5_with_the_commit_that_put_it_in(self):
        first = self.commit(dvc("1"), "2026-10-01", "run.log.dvc")
        second = self.commit(dvc("2"), "2026-10-02", "run.log.dvc")
        # md5 の同じ書き直し（大きさだけが違う）は、md5 を入れたコミットを変えない
        self.commit(dvc("2", size=11), "2026-10-03", "run.log.dvc")
        fourth = self.commit(dvc("1"), "2026-10-04", "run.log.dvc")
        history = rp.dvc_history(self.dir / "run.log.dvc")
        self.assertEqual(history["now"], [("run.log", H["1"])])
        self.assertIsNone(history["state"])
        self.assertEqual(
            history["history"]["run.log"],
            [
                (H["1"], fourth, "2026-10-04"),
                (H["2"], second, "2026-10-02"),
                (H["1"], first, "2026-10-01"),
            ],
        )

    def test_uncommitted_and_untracked_are_said_not_hidden(self):
        first = self.commit(dvc("1"), "2026-10-01", "run.log.dvc")
        (self.dir / "run.log.dvc").write_text(dvc("3"))
        history = rp.dvc_history(self.dir / "run.log.dvc")
        self.assertEqual(history["state"], "uncommitted")
        self.assertEqual(history["now"], [("run.log", H["3"])])
        self.assertEqual(history["history"]["run.log"], [(H["1"], first, "2026-10-01")])
        (self.dir / "other.log.dvc").write_text(dvc("4"))
        lone = rp.dvc_history(self.dir / "other.log.dvc")
        self.assertEqual(lone["state"], "unknown")
        self.assertEqual(lone["history"], {"run.log": []})

    def test_the_history_stops_where_the_dvc_file_was_removed(self):
        self.commit(dvc("1"), "2026-10-01", "run.log.dvc")
        self.git("rm", "-q", "run.log.dvc")
        self.git("-c", "core.hooksPath=/dev/null", "commit", "-q", "-m", "chore: rm")
        again = self.commit(dvc("2"), "2026-10-03", "run.log.dvc")
        self.assertEqual(
            rp.dvc_history(self.dir / "run.log.dvc")["history"]["run.log"],
            [(H["2"], again, "2026-10-03")],
        )

    def test_main_takes_dvc_files_and_prints_the_history(self):
        import contextlib
        import io

        first = self.commit(dvc("1"), "2026-10-01", "run.log.dvc")
        second = self.commit(dvc("2"), "2026-10-02", "run.log.dvc")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            status = rp.main([str(self.dir / "run.log.dvc")])
        self.assertEqual(status, 0)
        lines = out.getvalue().splitlines()
        self.assertIn(f"  {H['2'][:8]}  run.log", lines)
        self.assertIn(f"    {H['2'][:8]}  {second} 2026-10-02", lines)
        self.assertIn(f"    {H['1'][:8]}  {first} 2026-10-01", lines)
        self.assertLess(
            lines.index(f"    {H['2'][:8]}  {second} 2026-10-02"),
            lines.index(f"    {H['1'][:8]}  {first} 2026-10-01"),
        )

    def test_main_says_uncommitted_on_the_out_line(self):
        import contextlib
        import io

        self.commit(dvc("1"), "2026-10-01", "run.log.dvc")
        (self.dir / "run.log.dvc").write_text(dvc("3"))
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rp.main([str(self.dir / "run.log.dvc")])
        self.assertIn(
            f"  {H['3'][:8]}  run.log  uncommitted", out.getvalue().splitlines()
        )


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
