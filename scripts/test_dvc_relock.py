"""dvc-relock.py の検査（ADR-0258 の「rebase の後の付け直し」）。部品（dvc-checkout.py、check-references-shell.sh、
check-dvc-lock.sh、dvc-repro.py）は呼ばれたことと引数を書く偽物に替え、dvc は PATH の偽物で実行する。

実行: python3 -m unittest discover -s scripts -p 'test_dvc_relock.py'
"""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "dvc-relock.py"
SPEC = importlib.util.spec_from_file_location("dvc_relock", SCRIPT)
rel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rel)

V = "examples/v/dvc.yaml"
LATTICE = f"{V}:lattice-1B"
COMPARE = f"{V}:compare-1B"
OPENMC = f"{V}:openmc-1B"

SHELL = {f"{V}:references-shell": ["always changed"]}
DEPS = {
    **SHELL,
    f"{V}:xs": [{"changed deps": {"binary-forms/xs.md5": "modified"}}],
    LATTICE: [{"changed deps": {"binary-forms/lattice.md5": "modified"}}],
}
AFTER_LATTICE = {
    **SHELL,
    f"{V}:xs": [{"changed deps": {"binary-forms/xs.md5": "modified"}}],
    COMPARE: [{"changed deps": {"examples/v/p1/lattice-1B.out": "modified"}}],
}
FORMS_CALL = (
    "dvc repro -s build dvc.yaml:python-ast dvc.yaml:binary-forms dvc.yaml:toml-values"
)

# repro -s <直下のステージ> は呼ばれたことを書き、python-ast は追跡している正規形を書き直す。status は status.json を出す
FAKE_DVC = """\
#!/bin/sh
echo "dvc $*" >> "{d}/calls"
case "$1" in
  repro)
    echo "Stage '$3' didn't change, skipping"
    [ "$3" = "$(cat "{d}/fail-repro" 2>/dev/null)" ] && echo "dvc failed on $3" >&2 && exit 5
    case "$*" in *dvc.yaml:python-ast*) echo rewritten > python-ast/a.py.ast;; esac
    exit 0;;
  status) [ -f "{d}/status-fails" ] && exit 3; cat "{d}/status.json";;
  *) exit 99;;
esac
"""

# 部品の偽物: 引数を書き、{d}/<名前>.code があればその終了コードで終わる。dvc-repro.py は after.json を status にする
FAKE_PART = """\
#!/bin/sh
echo "{name} $*" >> "{d}/calls"
code=$(cat "{d}/{name}.code" 2>/dev/null || echo 0)
if [ "{name}" = dvc-repro.py ] && [ "$code" = 0 ] && [ -f "{d}/after.json" ]; then cp "{d}/after.json" "{d}/status.json"; fi
exit "$code"
"""
PARTS = (
    "dvc-checkout.py",
    "check-references-shell.sh",
    "check-dvc-lock.sh",
    "dvc-repro.py",
)


class Stops(unittest.TestCase):
    def test_a_stage_with_only_changed_outs_stops(self):
        status = {
            **DEPS,
            f"{V}:xs": [{"changed outs": {"examples/v/vera.msgpack": "modified"}}],
            "examples/v/table.msgpack.dvc": [
                {"changed outs": {"examples/v/table.msgpack": "modified"}}
            ],
        }
        self.assertEqual(
            rel.outs_only(status), [f"{V}:xs", "examples/v/table.msgpack.dvc"]
        )

    def test_changed_outs_with_changed_deps_is_left_as_changed(self):
        status = {
            **DEPS,
            f"{V}:xs": [
                {"changed deps": {"binary-forms/xs.md5": "modified"}},
                {"changed outs": {"examples/v/vera.msgpack": "modified"}},
            ],
        }
        self.assertEqual(rel.outs_only(status), [])

    def test_always_changed_alone_is_not_a_stage_left_changed(self):
        self.assertEqual(
            rel.left_changed(DEPS),
            [(f"{V}:lattice-1B", ["changed deps"]), (f"{V}:xs", ["changed deps"])],
        )

    def test_unknown_change_kind_is_an_error(self):
        with self.assertRaises(ValueError):
            rel.left_changed({f"{V}:xs": [{"changed something": {}}]})


class Process(unittest.TestCase):
    """偽の部品と dvc を置いたリポジトリで実行し、どこで止まり、何を実行したかを見る"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        repo = self.d / "repo"
        (repo / "scripts").mkdir(parents=True)
        (repo / "python-ast").mkdir()
        (repo / "examples" / "v").mkdir(parents=True)
        # 利用者の git の hook（core.hooksPath）を一時のリポジトリに掛けない
        git = [
            "git",
            "-C",
            str(repo),
            "-c",
            "core.hooksPath=/dev/null",
            "-c",
            "user.name=t",
            "-c",
            "user.email=t@t",
        ]
        subprocess.run(["git", "init", "-q", "-b", "main", str(repo)], check=True)
        (repo / "python-ast" / "a.py.ast").write_text("old\n")
        (repo / "dvc.lock").write_text("schema: '2.0'\n")
        (repo / "examples" / "v" / "dvc.lock").write_text("schema: '2.0'\n")
        (repo / "examples" / "v" / "verification.md").write_text("record\n")
        subprocess.run([*git, "add", "."], check=True)
        subprocess.run([*git, "commit", "-q", "-m", "base"], check=True)
        subprocess.run([*git, "checkout", "-q", "-b", "topic"], check=True)
        (repo / "examples" / "v" / "dvc.lock").write_text("schema: '2.0'\nstages: {}\n")
        subprocess.run([*git, "commit", "-q", "-am", "topic"], check=True)
        self.script = repo / "scripts" / SCRIPT.name
        shutil.copy(SCRIPT, self.script)
        shutil.copy(SCRIPT.parent / "count-changed-stages.py", repo / "scripts")
        for name in PARTS:
            part = repo / "scripts" / name
            part.write_text(FAKE_PART.format(name=name, d=self.d))
            part.chmod(0o755)
        bindir = self.d / "bin"
        bindir.mkdir()
        fake = bindir / "dvc"
        fake.write_text(FAKE_DVC.format(d=self.d))
        fake.chmod(0o755)
        self.env = {**os.environ, "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}"}
        self.repo = repo

    def tearDown(self):
        self.tmp.cleanup()

    def exit_with(self, name, code=1):
        (self.d / f"{name}.code").write_text(f"{code}\n")

    def run_script(self, status, *args, cwd=None):
        (self.d / "status.json").write_text(json.dumps(status))
        return subprocess.run(
            [sys.executable, str(self.script), *args],
            cwd=cwd or self.repo,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )

    def calls(self):
        path = self.d / "calls"
        return path.read_text().splitlines() if path.exists() else []

    def ran(self, prefix):
        return [c for c in self.calls() if c.startswith(prefix)]

    def test_without_named_stages_only_the_build_and_normalized_forms_run(self):
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.ran("dvc repro"), [FORMS_CALL])
        self.assertEqual(self.ran("dvc-repro.py"), [])
        self.assertIn("RAN: none", r.stdout)
        self.assertIn(f"CHANGED: {LATTICE}: changed deps", r.stdout)
        self.assertIn(f"CHANGED: {V}:xs: changed deps", r.stdout)
        self.assertNotIn(
            "references-shell",
            "\n".join(x for x in r.stdout.splitlines() if x.startswith("CHANGED:")),
        )

    def test_the_build_and_normalized_forms_run_in_one_dvc_repro(self):
        # dvc は目標を DAG の順に並べる（binary-forms は build の下流）
        self.run_script(DEPS)
        self.assertEqual(self.ran("dvc repro"), [FORMS_CALL])

    def test_parts_run_in_the_order_of_the_relock(self):
        self.run_script(DEPS)
        order = [
            c.split()[0] if not c.startswith("dvc ") else " ".join(c.split()[:2])
            for c in self.calls()
        ]
        self.assertEqual(
            order,
            [
                # 読めない lock は、作業ツリーを書き換える前に止める（直すコミットの fixup と rebase が汚れた作業ツリーで止まる）
                "check-dvc-lock.sh",
                "dvc-checkout.py",
                "dvc repro",
                "dvc status",
                "check-references-shell.sh",
                "check-dvc-lock.sh",
            ],
        )

    def test_named_stages_run_in_order_through_the_single_item_route(self):
        (self.d / "after.json").write_text(json.dumps(AFTER_LATTICE))
        r = self.run_script(DEPS, "--run", LATTICE, "--run", COMPARE)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(
            self.ran("dvc-repro.py"),
            [f"dvc-repro.py -s {LATTICE}", f"dvc-repro.py -s {COMPARE}"],
        )
        self.assertEqual(self.ran("dvc repro"), [FORMS_CALL])
        self.assertIn(f"RAN: {LATTICE}", r.stdout)
        self.assertIn(f"RAN: {COMPARE}", r.stdout)
        # 実行した後の dvc status で、残した changed を出す
        self.assertEqual(len(self.ran("dvc status")), 2)
        self.assertIn(f"CHANGED: {COMPARE}: changed deps", r.stdout)
        self.assertNotIn(f"CHANGED: {LATTICE}", r.stdout)

    def test_allow_reference_goes_only_with_its_named_stage(self):
        r = self.run_script(
            DEPS, "--run", LATTICE, "--run", OPENMC, "--allow-reference", OPENMC
        )
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(
            self.ran("dvc-repro.py"),
            [
                f"dvc-repro.py -s {LATTICE}",
                f"dvc-repro.py -s --allow-reference {OPENMC} {OPENMC}",
            ],
        )

    def test_allow_reference_for_a_stage_not_named_is_refused_before_anything_runs(
        self,
    ):
        r = self.run_script(DEPS, "--run", LATTICE, "--allow-reference", OPENMC)
        self.assertEqual(r.returncode, 2)
        self.assertEqual(self.calls(), [])

    def test_a_named_stage_that_stops_stops_the_relock_before_the_next(self):
        self.exit_with("dvc-repro.py")
        r = self.run_script(DEPS, "--run", LATTICE, "--run", COMPARE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.ran("dvc-repro.py"), [f"dvc-repro.py -s {LATTICE}"])
        self.assertIn("STOP:", r.stderr)

    def test_outputs_that_do_not_come_back_stop_before_anything_runs(self):
        self.exit_with("dvc-checkout.py")
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(
            self.calls(), ["check-dvc-lock.sh examples/v/dvc.lock", "dvc-checkout.py "]
        )
        self.assertIn("STOP:", r.stderr)

    def test_a_failing_build_or_normalized_form_stops(self):
        (self.d / "fail-repro").write_text("build")
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.ran("dvc repro"), [FORMS_CALL])
        self.assertEqual(self.ran("dvc-repro.py"), [])
        self.assertIn("Stage 'build' didn't change", r.stdout + r.stderr)
        self.assertIn("dvc failed on build", r.stdout + r.stderr)

    def test_the_normalized_forms_print_only_their_names_when_they_pass(self):
        # 正規形は数百のステージで、dvc は飛ばしたものも 1 行ずつ出す
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("didn't change", r.stdout + r.stderr)

    def test_failing_dvc_status_stops(self):
        (self.d / "status-fails").write_text("")
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.ran("dvc-repro.py"), [])
        self.assertIn("STOP:", r.stderr)

    def test_a_stage_with_only_changed_outs_stops_without_running_named_stages(self):
        status = {
            **DEPS,
            f"{V}:xs": [{"changed outs": {"examples/v/vera.msgpack": "modified"}}],
        }
        r = self.run_script(status, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"{V}:xs", r.stderr)
        self.assertEqual(self.ran("dvc-repro.py"), [])

    def test_stale_references_shell_stops_without_running_named_stages(self):
        self.exit_with("check-references-shell.sh")
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.ran("dvc-repro.py"), [])
        self.assertIn("STOP:", r.stderr)

    def test_unreadable_lock_stops_before_the_worktree_is_rewritten(self):
        self.exit_with("check-dvc-lock.sh")
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.calls(), ["check-dvc-lock.sh examples/v/dvc.lock"])
        self.assertIn("STOP:", r.stderr)

    def test_a_stop_at_the_first_lock_check_lists_nothing_to_commit(self):
        # 読めない lock は作業ツリーの変更でありうる。それをコミットするファイルとして出さない
        (self.repo / "examples" / "v" / "dvc.lock").write_text("stages: {}\n")
        self.exit_with("check-dvc-lock.sh")
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 1)
        self.assertNotIn("COMMIT:", r.stdout)

    def test_locks_changed_since_the_fork_point_are_checked(self):
        # examples/v/dvc.lock はブランチのコミットで変わった（rebase の pick は pre-commit の hook を通らない）
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 0, r.stderr)
        # 始めと、dvc が lock を書き換えた後の 2 度
        self.assertEqual(
            self.ran("check-dvc-lock.sh"), ["check-dvc-lock.sh examples/v/dvc.lock"] * 2
        )

    def test_files_to_commit_are_those_the_relock_wrote(self):
        # rebase の hook が付け替えた記録（コミットしていない）は、このコマンドが書いたものでない
        (self.repo / "examples" / "v" / "verification.md").write_text("remapped\n")
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 0, r.stderr)
        commit = [x for x in r.stdout.splitlines() if x.startswith("COMMIT:")]
        self.assertEqual(commit, ["COMMIT: python-ast/a.py.ast"])

    def test_a_rerun_after_a_stop_still_lists_the_forms_and_locks_to_commit(self):
        # 前の実行が書き直した正規形は、この実行では中身が変わらないが、まだコミットしていない
        (self.repo / "python-ast" / "a.py.ast").write_text("rewritten\n")
        (self.repo / "examples" / "v" / "verification.md").write_text("remapped\n")
        r = self.run_script(DEPS)
        self.assertEqual(r.returncode, 0, r.stderr)
        commit = [x for x in r.stdout.splitlines() if x.startswith("COMMIT:")]
        self.assertEqual(commit, ["COMMIT: python-ast/a.py.ast"])

    def test_a_stop_after_the_forms_still_lists_what_ran_and_what_to_commit(self):
        (self.d / "after.json").write_text(json.dumps(AFTER_LATTICE))
        (self.d / "dvc-repro.py.code").write_text("0\n")
        # 2 つ目のステージで止める: 偽の dvc-repro.py は 1 度目の後に終了コードを 1 にする
        part = self.repo / "scripts" / "dvc-repro.py"
        part.write_text(
            part.read_text().replace(
                'exit "$code"', f'echo 1 > "{self.d}/dvc-repro.py.code"; exit "$code"'
            )
        )
        r = self.run_script(DEPS, "--run", LATTICE, "--run", COMPARE)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"RAN: {LATTICE}", r.stdout)
        self.assertNotIn(f"RAN: {COMPARE}", r.stdout)
        self.assertIn("COMMIT: python-ast/a.py.ast", r.stdout)

    def test_a_stage_left_with_only_changed_outs_after_the_runs_stops(self):
        after = {
            **AFTER_LATTICE,
            f"{V}:xs": [{"changed outs": {"examples/v/vera.msgpack": "modified"}}],
        }
        (self.d / "after.json").write_text(json.dumps(after))
        r = self.run_script(DEPS, "--run", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"{V}:xs", r.stderr)

    def test_status_is_scoped_to_the_given_dvc_yaml(self):
        r = self.run_script(DEPS, V)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.ran("dvc status"), [f"dvc status --json --with-deps {V}"])

    def test_running_outside_the_repository_root_is_refused(self):
        r = self.run_script(DEPS, cwd=self.repo / "examples" / "v")
        self.assertEqual(r.returncode, 1)
        self.assertIn("root of the repository", r.stderr)
        self.assertEqual(self.calls(), [])


if __name__ == "__main__":
    unittest.main()
