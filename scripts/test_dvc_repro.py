"""dvc-repro.py の検査（ADR-0257、`-s` の経路は ADR-0258）。数えは dvc の出力の文字列を受け、プロセスの検査は PATH の偽の dvc で実行する。

実行: python3 -m unittest discover -s scripts -p 'test_dvc_repro.py'
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

SCRIPT = Path(__file__).resolve().parent / "dvc-repro.py"
SPEC = importlib.util.spec_from_file_location("dvc_repro", SCRIPT)
rep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rep)

V = "examples/v/dvc.yaml"

# dvc 3.67 の dvc repro --dry -f の形（直下の dvc.yaml のステージは yaml を付けない）
DRY = f"""\
Running stage 'build':
> scripts/build-binaries.sh

Running stage 'binary-forms@xs':
> python3 scripts/normalize_binary.py bin/xs binary-forms/xs.md5

Running stage '{V}:references-shell':
> nix eval --raw ../..#devShells.x86_64-linux.references.outPath > references-shell.txt

Running stage '{V}:xs-nuclide@U-235':
> ../../bin/xs nuclide xs.toml U-235 xs-nuclide/U-235.msgpack

Running stage '{V}:xs-nuclide@U-238':
> ../../bin/xs nuclide xs.toml U-238 xs-nuclide/U-238.msgpack

Running stage '{V}:xs':
> ../../bin/xs merge xs.toml xs-nuclide vera.msgpack

Running stage '{V}:openmc-1B':
> nix develop ../..#references -c python3 openmc_model.py run p1/openmc-1B.toml p1/openmc-1B-model p1/openmc-1B

Running stage '{V}:compare-1B':
> python3 compare_k.py p1/openmc-1B lattice-1B.out
Use `dvc push` to send your updates to remote storage.
"""

DOT = f"""\
strict digraph {{
"build";
"binary-forms@xs";
"build" -> "binary-forms@xs";
"binary-forms@xs" -> "{V}:xs-nuclide@U-235";
"binary-forms@xs" -> "{V}:xs-nuclide@U-238";
"{V}:xs-nuclide@U-235" -> "{V}:xs";
"{V}:xs-nuclide@U-238" -> "{V}:xs";
"{V}:references-shell" -> "{V}:openmc-1B";
"{V}:openmc-1B" -> "{V}:compare-1B";
"{V}:xs" -> "{V}:compare-1B";
}}
"""

SHELL_ONLY = {f"{V}:references-shell": ["always changed"]}
XS_CHANGED = {
    **SHELL_ONLY,
    f"{V}:xs-nuclide@U-235": [{"changed deps": {"binary-forms/xs.md5": "modified"}}],
}
OPENMC_CHANGED = {
    **SHELL_ONLY,
    f"{V}:openmc-1B": [{"changed outs": {"examples/v/p1/openmc-1B": "deleted"}}],
}


LATTICE = f"{V}:lattice-1B"
# `dvc repro --dry -f -s <目標>` は目標のステージだけを並べる（foreach の組なら組の全部）
DRY_LATTICE = f"""\
Running stage '{LATTICE}':
> ../../bin/lattice p1/pincell-1B.toml > p1/lattice-1B.out
Use `dvc push` to send your updates to remote storage.
"""
OPENMC = f"{V}:openmc-1B"
DRY_OPENMC = f"""\
Running stage '{OPENMC}':
> nix develop ../..#references -c python3 openmc_model.py run p1/openmc-1B.toml p1/openmc-1B-model p1/openmc-1B
Use `dvc push` to send your updates to remote storage.
"""
# 上流の changed が deps だけ（outs は dvc.lock のとおり）
UPSTREAM_DEPS = {
    "build": [{"changed deps": {"crates": "modified"}}],
    f"{V}:xs": [{"changed deps": {"binary-forms/xs.md5": "modified"}}],
    LATTICE: [{"changed deps": {"binary-forms/lattice.md5": "modified"}}],
}
UPSTREAM_OUTS = {
    **UPSTREAM_DEPS,
    f"{V}:xs": [{"changed outs": {"examples/v/vera.msgpack": "deleted"}}],
}
NORMALIZED_STALE = {
    **UPSTREAM_DEPS,
    "python-ast@examples/v/compare_k.py": [
        {"changed deps": {"examples/v/compare_k.py": "modified"}}
    ],
}


def counted(status, shell_matches=True, allow=()):
    return rep.count(
        status, rep.parse_dot(DOT), rep.parse_dry_run(DRY), shell_matches, list(allow)
    )


class Parse(unittest.TestCase):
    def test_dry_run_lists_stages_in_order_with_their_commands(self):
        cmds = rep.parse_dry_run(DRY)
        self.assertEqual(
            list(cmds)[:3], ["build", "binary-forms@xs", f"{V}:references-shell"]
        )
        self.assertEqual(
            cmds[f"{V}:xs"], ["../../bin/xs merge xs.toml xs-nuclide vera.msgpack"]
        )

    def test_dot_gives_edges_from_upstream_to_downstream(self):
        edges = rep.parse_dot(DOT)
        self.assertIn(("build", "binary-forms@xs"), edges)
        self.assertEqual(len(edges), 8)

    def test_only_commands_in_the_references_shell_are_reference_stages(self):
        cmds = rep.parse_dry_run(DRY)
        self.assertEqual(
            [k for k in cmds if rep.is_reference(cmds[k])], [f"{V}:openmc-1B"]
        )


class CountStages(unittest.TestCase):
    def test_downstream_of_a_changed_stage_is_counted(self):
        c = counted(XS_CHANGED)
        self.assertEqual(
            c.runs,
            {
                f"{V}:references-shell": "changed: always changed",
                f"{V}:xs-nuclide@U-235": "changed: changed deps",
                f"{V}:xs": f"downstream of {V}:xs-nuclide@U-235",
                f"{V}:compare-1B": f"downstream of {V}:xs-nuclide@U-235",
            },
        )
        self.assertEqual((c.references, c.refused), ([], []))

    def test_matching_references_shell_does_not_rerun_the_references(self):
        self.assertEqual(list(counted(SHELL_ONLY).runs), [f"{V}:references-shell"])

    def test_stale_references_shell_reruns_the_references(self):
        c = counted(SHELL_ONLY, shell_matches=False)
        self.assertEqual(
            list(c.runs), [f"{V}:references-shell", f"{V}:openmc-1B", f"{V}:compare-1B"]
        )
        self.assertEqual(c.refused, [f"{V}:openmc-1B"])

    def test_reference_stage_is_refused_unless_allowed(self):
        self.assertEqual(counted(OPENMC_CHANGED).refused, [f"{V}:openmc-1B"])
        allowed = counted(OPENMC_CHANGED, allow=[f"{V}:openmc-1B"])
        self.assertEqual(
            (allowed.references, allowed.refused), ([f"{V}:openmc-1B"], [])
        )

    def test_a_foreach_group_allows_all_its_stages(self):
        self.assertTrue(rep.allows(f"{V}:openmc-neutron", f"{V}:openmc-neutron@U-235"))
        self.assertFalse(rep.allows(f"{V}:openmc", f"{V}:openmc-neutron@U-235"))

    def test_allowing_a_stage_outside_the_run_is_an_error(self):
        with self.assertRaises(rep.ReproError):
            counted(OPENMC_CHANGED, allow=["openmc-1B"])

    def test_changed_keys_outside_the_run_are_listed_apart(self):
        c = counted({**SHELL_ONLY, "examples/v/x.msgpack.dvc": [{"changed outs": {}}]})
        self.assertEqual(c.outside, ["examples/v/x.msgpack.dvc"])
        self.assertNotIn("examples/v/x.msgpack.dvc", c.runs)


def planned(status, dry=DRY_LATTICE, allow=()):
    return rep.plan_single(status, rep.parse_dry_run(dry), list(allow))


class PlanSingleItem(unittest.TestCase):
    def test_upstream_with_only_changed_deps_runs_the_target_alone(self):
        p = planned(UPSTREAM_DEPS)
        self.assertEqual(
            (p.runs, p.blocked, p.run_first, p.refused), ([LATTICE], [], [], [])
        )

    def test_upstream_with_changed_outs_command_or_checksum_blocks(self):
        for change in (
            {"changed outs": {"examples/v/vera.msgpack": "modified"}},
            {"changed command": {}},
            {"changed checksum": {}},
        ):
            with self.subTest(change=change):
                p = planned({**UPSTREAM_DEPS, f"{V}:xs": [change]})
                self.assertEqual(p.blocked, [(f"{V}:xs", sorted(change))])

    def test_target_is_matched_by_the_stage_key_dvc_lists_not_by_its_spelling(self):
        own = {LATTICE: [{"changed outs": {"examples/v/p1/lattice-1B.out": "deleted"}}]}
        # 目標は書き方（./examples/...）でなく、dvc repro --dry -f -s が並べた key で決まる
        p = planned({**UPSTREAM_DEPS, **own})
        self.assertEqual(p.blocked, [])
        self.assertEqual(rep.stage_target(rep.parse_dry_run(DRY_LATTICE)), LATTICE)

    def test_stale_normalized_forms_are_run_first_rather_than_blocking(self):
        p = planned(NORMALIZED_STALE)
        self.assertEqual(p.blocked, [])
        self.assertEqual(p.run_first, ["python-ast@examples/v/compare_k.py"])

    def test_reference_target_is_refused_unless_allowed(self):
        self.assertEqual(planned(UPSTREAM_DEPS, DRY_OPENMC).refused, [OPENMC])
        allowed = planned(UPSTREAM_DEPS, DRY_OPENMC, allow=[OPENMC])
        self.assertEqual((allowed.references, allowed.refused), ([OPENMC], []))

    def test_allowing_a_stage_that_is_not_the_target_is_an_error(self):
        with self.assertRaises(rep.ReproError):
            planned(UPSTREAM_DEPS, allow=[OPENMC])

    def test_references_shell_upstream_is_reported_for_the_shell_check(self):
        p = planned({**UPSTREAM_DEPS, **SHELL_ONLY})
        self.assertEqual(p.blocked, [])
        self.assertTrue(p.needs_shell_check)
        self.assertFalse(planned(UPSTREAM_DEPS).needs_shell_check)


FAKE_DVC = """\
#!/bin/sh
echo "$*" >> "{d}/calls"
case "$1" in
  repro) case "$2" in
    --dry) cat "{d}/dry.txt";;
    *) if [ -f "{d}/after.json" ]; then mv "{d}/after.json" "{d}/status.json"; exit 0; fi; echo ran; exit 7;;
  esac;;
  dag) cat "{d}/dag.dot";;
  status) cat "{d}/status.json";;
  *) exit 99;;
esac
"""


class Process(unittest.TestCase):
    """偽の dvc を PATH に置き、止まるときは dvc repro を実行しないこと、通れば実行することを確かめる"""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        repo = self.d / "repo"
        (repo / "scripts").mkdir(parents=True)
        subprocess.run(["git", "init", "-q", str(repo)], check=True)
        self.script = repo / "scripts" / SCRIPT.name
        shutil.copy(SCRIPT, self.script)
        shutil.copy(SCRIPT.parent / "check-skip-upstream.py", repo / "scripts")
        self.check = repo / "scripts" / "check-references-shell.sh"
        self.set_shell_check(0)
        bindir = self.d / "bin"
        bindir.mkdir()
        fake = bindir / "dvc"
        fake.write_text(FAKE_DVC.format(d=self.d))
        fake.chmod(0o755)
        (self.d / "dry.txt").write_text(DRY)
        (self.d / "dag.dot").write_text(DOT)
        self.env = {**os.environ, "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}"}
        self.repo = repo

    def tearDown(self):
        self.tmp.cleanup()

    def set_shell_check(self, code):
        self.check.write_text(
            f'#!/bin/sh\necho shell >> "{self.d}/checks"\nexit {code}\n'
        )
        self.check.chmod(0o755)

    def run_script(self, status, *args):
        (self.d / "status.json").write_text(json.dumps(status))
        return subprocess.run(
            [sys.executable, str(self.script), *args],
            cwd=self.repo,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )

    def ran_repro(self):
        calls = (self.d / "calls").read_text().splitlines()
        return [c for c in calls if c.startswith("repro") and "--dry" not in c]

    def test_target_with_a_reference_stage_stops_without_running(self):
        r = self.run_script(OPENMC_CHANGED, "--max-stages", "10", V)
        self.assertEqual(r.returncode, 1)
        self.assertIn(f"RUN: {V}:openmc-1B [reference]", r.stdout)
        self.assertIn("not allowed by --allow-reference", r.stderr)
        self.assertEqual(self.ran_repro(), [])

    def test_target_over_the_limit_stops_without_running(self):
        r = self.run_script(XS_CHANGED, "--max-stages", "3", V)
        self.assertEqual(r.returncode, 1)
        self.assertIn("COUNT: stages=4 max=3", r.stdout)
        self.assertIn("more than --max-stages 3", r.stderr)
        self.assertEqual(self.ran_repro(), [])

    def test_stop_follows_the_list_in_a_log_that_merges_both_streams(self):
        (self.d / "status.json").write_text(json.dumps(XS_CHANGED))
        r = subprocess.run(
            [sys.executable, str(self.script), "--max-stages", "3", V],
            cwd=self.repo,
            env=self.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        self.assertLess(r.stdout.index("COUNT:"), r.stdout.index("STOP:"))

    def test_target_within_the_limit_runs_dvc_repro(self):
        r = self.run_script(XS_CHANGED, "--max-stages", "4", "--no-run-cache", V)
        self.assertEqual(r.returncode, 7)  # 偽の dvc repro の終了コード
        self.assertIn("ran", r.stdout)
        self.assertEqual(self.ran_repro(), [f"repro --no-run-cache {V}"])

    def test_allowed_reference_stage_runs(self):
        r = self.run_script(
            OPENMC_CHANGED,
            "--max-stages",
            "3",
            "--allow-reference",
            f"{V}:openmc-1B",
            V,
        )
        self.assertEqual(r.returncode, 7)
        self.assertEqual(self.ran_repro(), [f"repro {V}"])

    def test_count_only_does_not_run(self):
        r = self.run_script(XS_CHANGED, "--max-stages", "4", "--count-only", V)
        self.assertEqual(r.returncode, 0)
        self.assertEqual(self.ran_repro(), [])

    def test_stale_references_shell_stops_on_the_references(self):
        self.set_shell_check(1)
        r = self.run_script(SHELL_ONLY, "--max-stages", "10", V)
        self.assertEqual(r.returncode, 1)
        self.assertIn("matches the references devShell: no", r.stdout)
        self.assertEqual(self.ran_repro(), [])

    def test_single_item_runs_the_target_alone_when_upstream_outs_match_the_lock(self):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        r = self.run_script(UPSTREAM_DEPS, "-s", "--no-run-cache", LATTICE)
        self.assertEqual(r.returncode, 7)  # 偽の dvc repro の終了コード
        self.assertEqual(self.ran_repro(), [f"repro -s --no-run-cache {LATTICE}"])
        self.assertIn(f"RUN: {LATTICE}", r.stdout)

    def test_single_item_does_not_take_the_targets_own_change_for_an_upstream_one(self):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        own = {LATTICE: [{"changed outs": {"examples/v/p1/lattice-1B.out": "deleted"}}]}
        r = self.run_script({**UPSTREAM_DEPS, **own}, "-s", "./" + LATTICE)
        self.assertEqual(r.returncode, 7)
        self.assertEqual(self.ran_repro(), [f"repro -s ./{LATTICE}"])

    def test_single_item_stops_without_running_on_upstream_changes_other_than_deps(
        self,
    ):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        for change in (
            {"changed outs": {"examples/v/vera.msgpack": "deleted"}},
            {"changed command": {}},
            {"changed checksum": {}},
        ):
            with self.subTest(change=change):
                (self.d / "calls").unlink(missing_ok=True)
                r = self.run_script(
                    {**UPSTREAM_DEPS, f"{V}:xs": [change]}, "-s", LATTICE
                )
                self.assertEqual(r.returncode, 1)
                self.assertIn(f"BLOCKED: {V}:xs: {next(iter(change))}", r.stdout)
                self.assertEqual(self.ran_repro(), [])

    def test_single_item_runs_stale_normalized_forms_before_the_target(self):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        (self.d / "after.json").write_text(json.dumps(UPSTREAM_DEPS))
        r = self.run_script(NORMALIZED_STALE, "-s", LATTICE)
        self.assertEqual(r.returncode, 7)
        self.assertEqual(
            self.ran_repro(),
            ["repro -s python-ast@examples/v/compare_k.py", f"repro -s {LATTICE}"],
        )

    def test_single_item_keeps_no_run_cache_and_checks_the_shell_once_around_normalized_forms(
        self,
    ):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        (self.d / "after.json").write_text(json.dumps({**UPSTREAM_DEPS, **SHELL_ONLY}))
        r = self.run_script(
            {**NORMALIZED_STALE, **SHELL_ONLY}, "-s", "--no-run-cache", LATTICE
        )
        self.assertEqual(r.returncode, 7)
        self.assertEqual(
            self.ran_repro(),
            [
                "repro -s --no-run-cache python-ast@examples/v/compare_k.py",
                f"repro -s --no-run-cache {LATTICE}",
            ],
        )
        self.assertEqual((self.d / "checks").read_text().splitlines(), ["shell"])

    def test_single_item_runs_nothing_when_other_upstream_blocks_besides_normalized_forms(
        self,
    ):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        r = self.run_script({**NORMALIZED_STALE, **UPSTREAM_OUTS}, "-s", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(self.ran_repro(), [])

    def test_single_item_reference_target_needs_allow_reference(self):
        (self.d / "dry.txt").write_text(DRY_OPENMC)
        r = self.run_script(UPSTREAM_DEPS, "-s", OPENMC)
        self.assertEqual(r.returncode, 1)
        self.assertIn("not allowed by --allow-reference", r.stderr)
        self.assertNotIn("BLOCKED", r.stderr)
        self.assertEqual(self.ran_repro(), [])
        r = self.run_script(UPSTREAM_DEPS, "-s", "--allow-reference", OPENMC, OPENMC)
        self.assertEqual(r.returncode, 7)
        self.assertEqual(self.ran_repro(), [f"repro -s {OPENMC}"])

    def test_single_item_stops_when_the_references_shell_is_stale(self):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        self.set_shell_check(1)
        r = self.run_script({**UPSTREAM_DEPS, **SHELL_ONLY}, "-s", LATTICE)
        self.assertEqual(r.returncode, 1)
        self.assertIn("matches the references devShell: no", r.stdout)
        self.assertEqual(self.ran_repro(), [])

    def test_single_item_count_only_checks_without_running(self):
        (self.d / "dry.txt").write_text(DRY_LATTICE)
        r = self.run_script(NORMALIZED_STALE, "-s", "--count-only", LATTICE)
        self.assertEqual(r.returncode, 0)
        self.assertIn("RUN-FIRST: python-ast@examples/v/compare_k.py", r.stdout)
        self.assertEqual(self.ran_repro(), [])

    def test_single_item_takes_one_stage_target_and_no_stage_count(self):
        for args in (
            ["-s", "--max-stages", "1", LATTICE],
            ["-s", LATTICE, OPENMC],
            ["-s", V],
            ["-s", "-P"],
        ):
            with self.subTest(args=args):
                r = self.run_script(UPSTREAM_DEPS, *args)
                self.assertEqual(r.returncode, 2)
                self.assertFalse((self.d / "calls").exists())

    def test_running_outside_the_repository_root_is_refused(self):
        # dvc の目標は cwd からの相対で、ステージの key はリポジトリの直下からのパス
        sub = self.repo / "examples" / "v"
        sub.mkdir(parents=True)
        (self.d / "status.json").write_text(json.dumps(XS_CHANGED))
        r = subprocess.run(
            [sys.executable, str(self.script), "--max-stages", "4", "dvc.yaml"],
            cwd=sub,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(r.returncode, 1)
        self.assertIn("root of the repository", r.stderr)
        self.assertFalse((self.d / "calls").exists())

    def test_max_stages_is_required(self):
        r = self.run_script(XS_CHANGED, V)
        self.assertEqual(r.returncode, 2)
        self.assertFalse((self.d / "calls").exists())


if __name__ == "__main__":
    unittest.main()
