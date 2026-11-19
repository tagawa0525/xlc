"""check-readable-outputs.py の検査（ADR-0258 の「形式を変える修正の検査」）。読み手は、`bad` で始まるファイルを読めない
とする偽物に替え、一時の git のリポジトリの dvc.lock、.dvc、キャッシュで実行する。本番の読み手での検査は
crates/outputs/tests/check.rs にある。

実行: python3 -m unittest discover -s scripts -p 'test_check_readable_outputs.py'
"""

import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "check-readable-outputs.py"

# 標準入力の行を {d}/read に書き足し、`bad` で始まるファイルを読めないとする
FAKE_READER = """\
#!{python}
import sys
bad = False
with open("{d}/read", "a") as log:
    for line in sys.stdin:
        kind, path = line.rstrip("\\n").split("\\t")
        log.write(f"{{kind}}\\t{{path}}\\n")
        with open(path, "rb") as f:
            if f.read().startswith(b"bad"):
                bad = True
                print(f"UNREADABLE\\t{{path}}\\tcannot read it")
sys.exit(int(open("{d}/code").read()) if __import__("os").path.exists("{d}/code") else 1 if bad else 0)
"""

ROOT_YAML = """\
stages:
  build:
    cmd: scripts/build-binaries.sh
    deps:
      - scripts/build-binaries.sh
      - crates
      - Cargo.lock
    outs:
      - bin/lattice
"""
ROOT_LOCK = """\
schema: '2.0'
stages:
  build:
    cmd: scripts/build-binaries.sh
    outs:
    - path: bin/lattice
      hash: md5
      md5: {bin}
"""
# lib と parts と axial は上流、lattice と assemble と core は xlc の実行ファイルを呼ぶステージ、compare は呼ばない。
# unrun は流していない（outs に md5 が無い）
LOCK = """\
schema: '2.0'
stages:
  lib:
    cmd: python3 make.py
    outs:
    - path: lib.msgpack
      hash: md5
      md5: {lib}
  gen:
    cmd: python3 gen.py
    outs:
    - path: gen.toml
      hash: md5
      md5: {gen}
    - path: gen.json
      hash: md5
      md5: {gen}
  range:
    cmd: ../../bin/xs range xs.toml xs-range.toml
    deps:
    - path: ../../binary-forms/xs.md5
      hash: md5
      md5: {form}
    outs:
    - path: xs-range.toml
      hash: md5
      md5: {range}
  nominal:
    cmd: ../../bin/lattice nominal.toml
    deps:
    - path: ../../binary-forms/lattice.md5
      hash: md5
      md5: {form}
    outs:
    - path: nominal-state.toml
      hash: md5
      md5: {gen}
  part@U235:
    cmd: ../../bin/xs nuclide
    outs:
    - path: parts/U235.msgpack
      hash: md5
      md5: {part}
  axial:
    cmd: ../../bin/lattice axial.toml
    outs:
    - path: axial
      hash: md5
      md5: {axial}.dir
  pins:
    cmd: ../../bin/lattice pins.toml
    outs:
    - path: pins.msgpack
      hash: md5
      md5: {pins}
  unrun:
    cmd: ../../bin/lattice unrun.toml
    outs:
    - path: unrun.msgpack
      hash: md5
  lattice:
    cmd: ../../bin/lattice lattice.toml
    deps:
    - path: ../../toml-values/examples/a/lattice.toml.json
      hash: md5
      md5: {form}
    - path: lib.msgpack
      hash: md5
      md5: {lib}
    - path: ../../bin/lattice
      hash: md5
      md5: {bin}
    - path: ../../binary-forms/lattice.md5
      hash: md5
      md5: {form}
    - path: art.msgpack
      hash: md5
      md5: {art}
    - path: unrun.msgpack
      hash: md5
      md5: {lib}
    - path: gen.toml
      hash: md5
      md5: {gen}
    - path: gen.json
      hash: md5
      md5: {gen}
    outs:
    - path: lattice.out
      hash: md5
      md5: {form}
  assemble:
    cmd: ../../bin/xs assemble xs.toml xs-range.toml parts
    deps:
    - path: xs-range.toml
      hash: md5
      md5: {range}
    - path: parts
      hash: md5
      md5: {form}.dir
    - path: ../../binary-forms/xs.md5
      hash: md5
      md5: {form}
    outs:
    - path: assembled.msgpack
      hash: md5
      md5: {form}
  core:
    cmd: ../../bin/core core.toml
    deps:
    - path: axial/t2.msgpack
      hash: md5
      md5: {t2}
    - path: nominal-state.toml
      hash: md5
      md5: {gen}
    - path: ../../binary-forms/core.md5
      hash: md5
      md5: {form}
    outs:
    - path: core.out
      hash: md5
      md5: {form}
  core-all:
    cmd: ../../bin/core core-all.toml
    deps:
    - path: axial
      hash: md5
      md5: {axial}.dir
    - path: ../../binary-forms/core.md5
      hash: md5
      md5: {form}
    outs:
    - path: core-all.out
      hash: md5
      md5: {form}
  compare:
    cmd: python3 compare.py pins.msgpack
    deps:
    - path: pins.msgpack
      hash: md5
      md5: {pins}
    - path: ../../python-ast/examples/a/compare.py.ast
      hash: md5
      md5: {form}
    outs:
    - path: compare.json
      hash: md5
      md5: {form}
"""
DVC_FILE = """\
outs:
- md5: {art}
  size: 3
  hash: md5
  path: art.msgpack
"""
MD5 = {
    name: f"{i:02x}" + "a" * 30
    for i, name in enumerate(
        [
            "lib",
            "range",
            "part",
            "axial",
            "pins",
            "form",
            "bin",
            "art",
            "t1",
            "t2",
            "gen",
        ],
        start=16,
    )
}
OUTPUTS = {
    "lib": "examples/a/lib.msgpack",
    "range": "examples/a/xs-range.toml",
    "part": "examples/a/parts/U235.msgpack",
    "t1": "examples/a/axial/t1.msgpack",
    "t2": "examples/a/axial/t2.msgpack",
    "art": "examples/a/art.msgpack",
}
ISSUE = """\
---
type: Issue
status: {status}
labels: [dvc]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# 読めない出力
"""


# マージの hook の中の run-checks.sh から呼ばれると、git が GIT_DIR や GIT_INDEX_FILE を渡す。一時のリポジトリの git と
# 検査がそれを継ぐと、本物のリポジトリを読み書きする
ENV = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def git(root, *args):
    # 利用者の git の hook（core.hooksPath）を一時のリポジトリに掛けない
    command = ["git", "-C", str(root), "-c", "core.hooksPath=/dev/null", *args]
    return subprocess.run(
        command, check=True, capture_output=True, text=True, env=ENV
    ).stdout


class CheckReadableOutputs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "repo"
        self.d = Path(self.tmp.name) / "fake"
        self.cache = Path(self.tmp.name) / "cache"
        self.d.mkdir()
        files = {
            "dvc.yaml": ROOT_YAML,
            "dvc.lock": ROOT_LOCK.format(**MD5),
            "examples/a/dvc.lock": LOCK.format(**MD5),
            "examples/a/art.msgpack.dvc": DVC_FILE.format(**MD5),
            "crates/x.rs": "fn main() {}\n",
            "Cargo.lock": "lock\n",
            "docs/issues/0001-open.md": ISSUE.format(status="open"),
            "docs/issues/0002-closed.md": ISSUE.format(status="closed"),
        }
        for path, text in files.items():
            (self.root / path).parent.mkdir(parents=True, exist_ok=True)
            (self.root / path).write_text(text)
        git(self.root, "init", "-q", "-b", "main")
        git(self.root, "add", ".")
        git(
            self.root,
            "-c",
            "user.name=t",
            "-c",
            "user.email=t@t",
            "commit",
            "-qm",
            "base",
        )
        for name in ["lib", "range", "part", "pins", "art", "t1", "t2", "bin", "gen"]:
            self.write_cache(MD5[name], b"ok")
        manifest = [
            {"md5": MD5["t1"], "relpath": "t1.msgpack"},
            {"md5": MD5["t2"], "relpath": "t2.msgpack"},
        ]
        self.write_cache(MD5["axial"] + ".dir", json.dumps(manifest).encode())
        reader = self.d / "reader"
        reader.write_text(FAKE_READER.format(python=sys.executable, d=self.d))
        reader.chmod(0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def cache_path(self, md5):
        return self.cache / "files" / "md5" / md5[:2] / md5[2:]

    def write_cache(self, md5, data):
        path = self.cache_path(md5)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)

    def run_check(self, *args):
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--cache",
                str(self.cache),
                "--reader",
                str(self.d / "reader"),
                *args,
            ],
            cwd=self.root,
            capture_output=True,
            text=True,
            check=False,
            env=ENV,
        )

    def read_lines(self):
        path = self.d / "read"
        return path.read_text().splitlines() if path.exists() else []

    def exempt(self, text):
        (self.root / "scripts").mkdir(exist_ok=True)
        (self.root / "scripts/unreadable-outputs.toml").write_text(
            textwrap.dedent(text)
        )

    def test_reads_every_upstream_output_an_xlc_stage_depends_on_once(self):
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        expected = {
            f"{'xs-range' if name == 'range' else 'msgpack'}\t{self.cache_path(MD5[name])}"
            for name in OUTPUTS
        }
        self.assertEqual(sorted(self.read_lines()), sorted(expected))
        self.assertIn("READ: 6 files", result.stdout)

    def test_an_unreadable_output_stops_with_its_path_and_the_stage_that_wrote_it(self):
        self.write_cache(MD5["lib"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            f"UNREADABLE: examples/a/lib.msgpack (md5 {MD5['lib']}, written by "
            "examples/a/dvc.yaml:lib): cannot read it",
            result.stdout,
        )

    def test_a_file_in_a_directory_output_names_the_directory_stage(self):
        self.write_cache(MD5["t1"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            f"UNREADABLE: examples/a/axial/t1.msgpack (md5 {MD5['t1']}, written by "
            "examples/a/dvc.yaml:axial)",
            result.stdout,
        )

    def test_a_file_under_a_directory_dependency_names_its_own_stage(self):
        self.write_cache(MD5["part"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("written by examples/a/dvc.yaml:part@U235", result.stdout)

    def test_an_artifact_names_its_dvc_file(self):
        self.write_cache(MD5["art"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("written by examples/a/art.msgpack.dvc", result.stdout)

    def test_outputs_only_other_tools_read_are_not_read(self):
        self.write_cache(MD5["pins"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn(
            str(self.cache_path(MD5["pins"])), "\n".join(self.read_lines())
        )

    def test_files_in_no_xlc_format_are_not_read(self):
        self.write_cache(MD5["gen"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn(str(self.cache_path(MD5["gen"])), "\n".join(self.read_lines()))
        # Python のステージが書く TOML と JSON、xlc のステージが書く範囲でない TOML
        self.assertIn("3 outputs in no xlc format", result.stdout)

    def test_an_output_tied_to_an_open_issue_passes_and_others_still_stop(self):
        self.write_cache(MD5["lib"], b"bad")
        self.exempt(f"""\
            [[output]]
            md5 = "{MD5["lib"]}"
            issue = "#0001"
            path = "examples/a/lib.msgpack"
            """)
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("ALLOWED: examples/a/lib.msgpack", result.stdout)
        self.assertIn("#0001", result.stdout)
        self.write_cache(MD5["art"], b"bad")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("UNREADABLE: examples/a/art.msgpack", result.stdout)

    def test_a_directory_md5_ties_every_file_in_it(self):
        self.write_cache(MD5["t1"], b"bad")
        self.write_cache(MD5["t2"], b"bad")
        self.exempt(f"""\
            [[output]]
            md5 = "{MD5["axial"]}.dir"
            issue = "#0001"
            path = "examples/a/axial"
            """)
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout.count("ALLOWED:"), 2, result.stdout)

    def test_an_entry_tied_to_a_closed_or_missing_issue_stops(self):
        for issue in ["#0002", "#0003", "0001"]:
            with self.subTest(issue=issue):
                self.exempt(f"""\
                    [[output]]
                    md5 = "{MD5["lib"]}"
                    issue = "{issue}"
                    path = "examples/a/lib.msgpack"
                    """)
                result = self.run_check()
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("STOP:", result.stdout)

    def test_a_malformed_entry_stops(self):
        for text in [
            'output = ["x"]\n',
            '[[output]]\nmd5 = 123\nissue = "#0001"\npath = "p"\n',
            '[[output]]\nmd5 = "{lib}"\nissue = 1\npath = "p"\n'.format(**MD5),
            '[[output]]\nmd5 = "{lib}"\nissue = "#0001"\npath = 1\n'.format(**MD5),
        ]:
            with self.subTest(text=text):
                self.exempt(text)
                result = self.run_check()
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertIn("STOP:", result.stdout)
                self.assertNotIn("Traceback", result.stderr)

    def test_two_entries_with_the_same_md5_stop(self):
        entry = '[[output]]\nmd5 = "{lib}"\nissue = "#0001"\npath = "p"\n'.format(**MD5)
        self.exempt(entry + entry)
        result = self.run_check()
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("STOP:", result.stdout)

    def test_an_output_missing_from_the_cache_stops(self):
        self.cache_path(MD5["t2"]).unlink()
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("MISSING: examples/a/axial/t2.msgpack", result.stdout)

    def test_a_reader_that_fails_to_run_stops(self):
        (self.d / "code").write_text("2")
        result = self.run_check("--all")
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("STOP:", result.stdout)

    def test_a_branch_that_does_not_change_the_build_deps_reads_nothing(self):
        (self.root / "docs/note.md").write_text("note\n")
        result = self.run_check("--base", "main")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("SKIP:", result.stdout)
        self.assertEqual(self.read_lines(), [])

    def test_moving_a_file_out_of_the_build_deps_reads_the_outputs(self):
        git(self.root, "mv", "crates/x.rs", "docs/x.rs")
        result = self.run_check("--base", "main")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("READ: 6 files", result.stdout)

    def test_an_untracked_file_in_the_build_deps_reads_the_outputs(self):
        (self.root / "crates/new.rs").write_text("fn new() {}\n")
        result = self.run_check("--base", "main")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("READ: 6 files", result.stdout)

    def test_a_branch_that_changes_the_build_deps_reads_the_outputs(self):
        for path in ["crates/x.rs", "Cargo.lock"]:
            with self.subTest(path=path):
                git(self.root, "checkout", "-q", "--", ".")
                (self.root / path).write_text("changed\n")
                result = self.run_check("--base", "main")
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("READ: 6 files", result.stdout)


if __name__ == "__main__":
    os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    unittest.main()
