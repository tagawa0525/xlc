"""normalized_forms.py の検査（pre-commit の hook が Python と TOML の正規形を書き直す。ADR-0258 の「正規形の置き場」）。
dvc は PATH の偽物で実行する。

実行: python3 -m unittest discover -s scripts -p 'test_normalized_forms.py'
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import normalized_forms as forms

REPRO_CALL = "dvc repro -s dvc.yaml:python-ast dvc.yaml:toml-values"

ROOT_YAML = """\
stages:
  build:
    cmd: scripts/build-binaries.sh
  # Python の正規形
  python-ast:
    foreach:
{python}
    do:
      cmd: python3 scripts/normalize_python.py ${{item}} python-ast/${{item}}.ast
      deps:
        - scripts/normalize_python.py
        - ${{item}}
      outs:
        - python-ast/${{item}}.ast:
            cache: false
  # TOML の正規形
  toml-values:
    foreach:
{toml}
    do:
      cmd: python3 scripts/normalize_toml.py ${{item}} toml-values/${{item}}.json
      deps:
        - python-ast/scripts/normalize_toml.py.ast
        - ${{item}}
      outs:
        - toml-values/${{item}}.json:
            cache: false
"""

EXAMPLE_YAML = """\
vars:
  - cases:
      a: {toml: a.toml}
stages:
  compare:
    cmd: python3 compare.py
    deps:
      - ../../python-ast/examples/v/compare.py.ast
  run:
    foreach: ${cases}
    do:
      cmd: lattice ${item.toml}
      deps:
        - ../../toml-values/examples/v/${item.toml}.json
"""

PYTHON = ["examples/v/compare.py", "scripts/normalize_toml.py"]
TOML = ["examples/v/a.toml"]


def lock_entry(name: str, source: str) -> str:
    """直下の dvc.lock の 1 項（dvc の書き方: 項は 2 文字、その中は 4 文字の字下げ、長い cmd は折り返す）"""
    return (
        f"  {name}:\n"
        f"    cmd: python3 scripts/normalize.py {source}\n"
        f"      out/{source}\n"
        "    deps:\n"
        f"    - path: {source}\n"
        "      hash: md5\n"
        "      md5: 0123\n"
        "      size: 1\n"
    )


def root_yaml(python: list[str], toml: list[str]) -> str:
    return ROOT_YAML.format(
        python="\n".join(f"      - {x}" for x in python),
        toml="\n".join(f"      - {x}" for x in toml),
    )


# repro は呼ばれたことを書き、{d}/lock があればそれを dvc.lock に写し、{d}/code があればその終了コードで終わる
FAKE_DVC = """\
#!/bin/sh
echo "dvc $*" >> "{d}/calls"
[ -f "{d}/lock" ] && cp "{d}/lock" dvc.lock
echo "Stage 'python-ast@x' didn't change, skipping"
exit "$(cat "{d}/code" 2>/dev/null || echo 0)"
"""


class Repository(unittest.TestCase):
    """正規形の hook が見るリポジトリ（直下の dvc.yaml と、例の dvc.yaml）を一時のディレクトリに作る"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        self.repo = self.d / "repo"
        (self.repo / "examples" / "v").mkdir(parents=True)
        (self.repo / "scripts").mkdir()
        subprocess.run(["git", "init", "-q", "-b", "main", str(self.repo)], check=True)
        self.write("dvc.yaml", root_yaml(PYTHON, TOML))
        self.write("examples/v/dvc.yaml", EXAMPLE_YAML)
        self.write(
            "dvc.lock",
            "schema: '2.0'\nstages:\n"
            + lock_entry("python-ast@examples/v/compare.py", "examples/v/compare.py")
            + lock_entry("toml-values@examples/v/a.toml", "examples/v/a.toml"),
        )
        self.write("python-ast/examples/v/compare.py.ast", "x = 1\n")
        self.write("python-ast/scripts/normalize_toml.py.ast", "y = 1\n")
        self.write("toml-values/examples/v/a.toml.json", "{}\n")
        self.script = self.repo / "scripts" / "normalized_forms.py"
        shutil.copy(SCRIPTS / "normalized_forms.py", self.script)
        self.git("add", ".")
        self.git("commit", "-q", "-m", "base")
        bindir = self.d / "bin"
        bindir.mkdir()
        fake = bindir / "dvc"
        fake.write_text(FAKE_DVC.format(d=self.d))
        fake.chmod(0o755)
        self.env = {**os.environ, "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}"}

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, path: str, text: str) -> None:
        target = self.repo / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    def read(self, path: str) -> str:
        return (self.repo / path).read_text(encoding="utf-8")

    def git(self, *args: str) -> str:
        # 利用者の git の hook（core.hooksPath）を一時のリポジトリに掛けない
        return subprocess.run(
            [
                "git",
                "-C",
                str(self.repo),
                "-c",
                "core.hooksPath=/dev/null",
                "-c",
                "user.name=t",
                "-c",
                "user.email=t@t",
                *args,
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout

    def run_hook(self, cwd: Path | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(self.script)],
            cwd=cwd or self.repo,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )

    def calls(self) -> list[str]:
        path = self.d / "calls"
        return path.read_text().splitlines() if path.exists() else []

    def rewritten(self, r: subprocess.CompletedProcess) -> list[str]:
        return [x for x in r.stdout.splitlines() if x.startswith("REWRITTEN:")]


class Items(Repository):
    def test_items_are_the_sources_of_the_forms_the_stages_depend_on(self):
        # 例の deps（vars の foreach を解いたものを含む）と、直下のほかのステージの deps（toml-values が読む
        # normalize_toml.py の正規形）から。正規形を書くステージ自身の deps は数えない
        self.assertEqual(
            forms.form_items(self.repo),
            {"python-ast": PYTHON, "toml-values": TOML},
        )


class Hook(Repository):
    def test_a_tree_in_step_passes_without_writing(self):
        r = self.run_hook()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.rewritten(r), [])
        self.assertEqual(self.git("status", "--porcelain"), "")

    def test_the_python_and_toml_stages_run_in_one_dvc_repro(self):
        # toml-values は python-ast の normalize_toml.py.ast を読むので後に挙げる（dvc repro -s は挙げた順に実行する）
        self.run_hook()
        self.assertEqual(self.calls(), [REPRO_CALL])

    def test_the_skipped_stages_are_not_printed_when_dvc_passes(self):
        # 正規形は数百のステージで、dvc は飛ばしたものも 1 行ずつ出す
        r = self.run_hook()
        self.assertNotIn("didn't change", r.stdout + r.stderr)

    def test_files_the_stages_rewrite_stop_the_commit_and_are_listed(self):
        (self.d / "lock").write_text(self.read("dvc.lock").replace("0123", "4567"))
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.rewritten(r), ["REWRITTEN: dvc.lock"])
        self.assertIn("git add", r.stderr)

    def test_the_commit_passes_again_once_the_rewritten_files_are_added(self):
        (self.d / "lock").write_text(self.read("dvc.lock").replace("0123", "4567"))
        self.assertEqual(self.run_hook().returncode, 1)
        self.git("add", "dvc.lock")
        r = self.run_hook()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.rewritten(r), [])

    def test_a_form_added_to_an_example_is_added_to_the_foreach(self):
        self.write(
            "examples/v/dvc.yaml",
            EXAMPLE_YAML.replace(
                "      - ../../python-ast/examples/v/compare.py.ast\n",
                "      - ../../python-ast/examples/v/compare.py.ast\n"
                "      - ../../python-ast/examples/common.py.ast\n",
            ),
        )
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.rewritten(r), ["REWRITTEN: dvc.yaml"])
        self.assertEqual(
            self.read("dvc.yaml"),
            root_yaml(["examples/common.py", *PYTHON], TOML),
        )
        self.assertEqual(self.calls(), [REPRO_CALL])

    def test_a_form_no_stage_depends_on_leaves_the_foreach_the_lock_and_the_tree(self):
        self.write("examples/v/dvc.yaml", EXAMPLE_YAML.replace("a.toml", "b.toml"))
        self.write("toml-values/examples/v/b.toml.json", "{}\n")
        self.git("add", ".")
        self.git("commit", "-q", "-m", "rename")
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertEqual(
            self.rewritten(r),
            [
                "REWRITTEN: dvc.lock",
                "REWRITTEN: dvc.yaml",
                "REWRITTEN: toml-values/examples/v/a.toml.json",
            ],
        )
        self.assertEqual(
            self.read("dvc.yaml"), root_yaml(PYTHON, ["examples/v/b.toml"])
        )
        self.assertEqual(
            self.read("dvc.lock"),
            "schema: '2.0'\nstages:\n"
            + lock_entry("python-ast@examples/v/compare.py", "examples/v/compare.py"),
        )
        self.assertFalse((self.repo / "toml-values/examples/v/a.toml.json").exists())

    def test_a_failing_dvc_repro_stops_and_shows_its_output(self):
        (self.d / "code").write_text("3\n")
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertIn("STOP:", r.stderr)
        self.assertIn("didn't change", r.stdout + r.stderr)

    def test_a_stop_still_lists_what_was_rewritten_before_it(self):
        self.write("examples/v/dvc.yaml", EXAMPLE_YAML.replace("compare.py", "c.py"))
        (self.d / "code").write_text("3\n")
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertIn("REWRITTEN: dvc.yaml", r.stdout)

    def test_a_stage_whose_forms_no_stage_depends_on_is_an_error(self):
        # 空の foreach を書かない（並びの無いステージは、手で消すか、deps を足して残すかを人が決める）
        self.write("examples/v/dvc.yaml", "stages:\n  compare:\n    cmd: echo\n")
        self.write("dvc.yaml", root_yaml(["examples/v/compare.py"], TOML))
        r = self.run_hook()
        self.assertEqual(r.returncode, 1)
        self.assertIn("toml-values", r.stderr)
        self.assertEqual(
            self.read("dvc.yaml"), root_yaml(["examples/v/compare.py"], TOML)
        )
        self.assertEqual(self.calls(), [])

    def test_running_outside_the_repository_root_is_refused(self):
        r = self.run_hook(cwd=self.repo / "examples" / "v")
        self.assertEqual(r.returncode, 1)
        self.assertIn("root of the repository", r.stderr)
        self.assertEqual(self.calls(), [])


if __name__ == "__main__":
    unittest.main()
