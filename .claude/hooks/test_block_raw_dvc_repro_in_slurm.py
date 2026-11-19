"""block-raw-dvc-repro-in-slurm.py の検査（ADR-0257）。

実行: cd .claude/hooks && python3 -m unittest test_block_raw_dvc_repro_in_slurm
"""

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HOOK = Path(__file__).resolve().parent / "block-raw-dvc-repro-in-slurm.py"
SPEC = importlib.util.spec_from_file_location("block_raw_dvc_repro_in_slurm", HOOK)
hook = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hook)

SLURM = (
    "~/.claude/scripts/slurm-run.sh --class trial --ref '#0001' --purpose 'p' -t 30 job"
)

BLOCKED = [
    f"{SLURM} 'dvc repro examples/vera/dvc.yaml'",
    f"{SLURM} 'nix develop -c dvc repro examples/c5g7/dvc.yaml'",
    f"cd /home/tagawa/github/xlc-x && {SLURM} 'dvc repro -P'",
    f"{SLURM} 'dvc repro -s examples/vera/dvc.yaml:xs'",
    f"{SLURM} 'dvc repro --no-run-cache examples/vera/dvc.yaml:compare-1C'",
    f"{SLURM} 'cd examples/vera && dvc repro dvc.yaml'",
    f"{SLURM} \"bash -c 'dvc repro examples/vera/dvc.yaml'\"",
    f"{SLURM} 'python3 examples/stage_runner.py examples/vera/dvc.yaml xs-nuclide && dvc repro examples/vera/dvc.yaml:xs'",
    f"{SLURM} 'python3 -m dvc repro examples/vera/dvc.yaml'",
    f"timeout 3600 {SLURM} 'dvc repro examples/vera/dvc.yaml'",
    f"bash -c \"{SLURM} 'dvc repro examples/vera/dvc.yaml'\"",
    "slurm-run.sh --class record --ref '- x' --purpose 'p' --timed -t 90 n 'dvc repro x'",
    f"{SLURM} 'dvc repro examples/vera/dvc.yaml' 2>&1 | tail -20",
    f"{SLURM} 'dvc repro examples/vera/dvc.yaml' 2>/dev/null",
]

ALLOWED = [
    "dvc repro examples/vera/dvc.yaml",
    "nix develop -c dvc repro -s dvc.yaml:python-ast",
    "cd /home/tagawa/github/xlc-x && nix develop -c dvc repro -s examples/vera/dvc.yaml:openmc-1B-model",
    f"{SLURM} 'scripts/dvc-repro.py --max-stages 10 examples/vera/dvc.yaml'",
    f"{SLURM} 'nix develop -c scripts/dvc-repro.py --max-stages 3 --allow-reference examples/vera/dvc.yaml:openmc-1B examples/vera/dvc.yaml'",
    f"{SLURM} 'python3 examples/stage_runner.py examples/vera/dvc.yaml xs-nuclide'",
    f"{SLURM} 'scripts/dvc-repro.py -s examples/vera/dvc.yaml:lattice-1B'",
    f"{SLURM} 'python3 examples/stage_runner.py examples/vera/dvc.yaml xs-nuclide && scripts/dvc-repro.py --max-stages 1 examples/vera/dvc.yaml:xs'",
    f"{SLURM} 'target/release/lattice x.toml'",
    f"{SLURM} 'echo dvc repro'",
    f"{SLURM} 'dvc status'",
    "grep -rn 'slurm-run.sh' .claude/skills",
    "echo \"slurm-run.sh -t 1 n 'dvc repro x'\"",
    "~/.claude/scripts/slurm-queue.sh",
    "~/.claude/scripts/slurm-run.sh",
    "",
]


class RawDvcReproInSlurmTest(unittest.TestCase):
    def test_slurm_jobs_that_run_dvc_repro_directly_are_blocked(self):
        for command in BLOCKED:
            with self.subTest(command=command):
                self.assertTrue(hook.submits_raw_dvc_repro(command))

    def test_wrapper_runner_and_local_runs_are_allowed(self):
        for command in ALLOWED:
            with self.subTest(command=command):
                self.assertFalse(hook.submits_raw_dvc_repro(command))

    def test_unbalanced_quotes_in_the_job_command_are_an_error(self):
        with self.assertRaises(ValueError):
            hook.submits_raw_dvc_repro(f'{SLURM} "echo \'dvc repro"')


def run_hook(payload):
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )


def bash(command):
    return {"tool_name": "Bash", "tool_input": {"command": command}}


class HookProcessTest(unittest.TestCase):
    def test_blocked_command_exits_2_with_the_wrapper_to_use(self):
        r = run_hook(bash(BLOCKED[0]))
        self.assertEqual(r.returncode, 2)
        self.assertIn("scripts/dvc-repro.py --max-stages", r.stderr)
        self.assertIn("--allow-reference", r.stderr)
        # 上流を飛ばす実行も、Slurm ではラッパーの -s の経路へ導く（ADR-0258）
        self.assertIn("scripts/dvc-repro.py -s", r.stderr)

    def test_allowed_command_exits_0_silently(self):
        r = run_hook(bash(ALLOWED[3]))
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_other_tools_pass(self):
        r = run_hook({"tool_name": "Read", "tool_input": {"file_path": "/x"}})
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_unparsable_command_is_reported_without_blocking(self):
        r = run_hook(bash(f'{SLURM} "echo \'dvc repro"'))
        self.assertEqual(r.returncode, 1)
        self.assertIn("読めない", r.stderr)

    def test_malformed_input_fails(self):
        r = subprocess.run(
            [sys.executable, str(HOOK)],
            input="not json",
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(r.returncode, 1)
        self.assertNotEqual(r.stderr, "")


SETTINGS = HOOK.parent.parent / "settings.json"
PROJECT = str(HOOK.parent.parent.parent)


def run_settings_command(payload, path=None):
    """.claude/settings.json の、この hook を呼ぶ command を、Claude Code と同じくシェルで実行する"""
    (entry,) = json.loads(SETTINGS.read_text())["hooks"]["PreToolUse"]
    (handler,) = [h for h in entry["hooks"] if HOOK.name in h["command"]]
    env = {**os.environ, "CLAUDE_PROJECT_DIR": PROJECT}
    if path is not None:
        env["PATH"] = path
    return subprocess.run(
        ["sh", "-c", handler["command"]],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


class SettingsCommandTest(unittest.TestCase):
    def test_settings_command_blocks_raw_dvc_repro_in_slurm(self):
        self.assertEqual(run_settings_command(bash(BLOCKED[0])).returncode, 2)

    def test_commands_without_slurm_run_skip_python(self):
        with tempfile.TemporaryDirectory() as d:
            fake = Path(d) / "python3"
            fake.write_text("#!/bin/sh\nexit 99\n")
            fake.chmod(0o755)
            path = f"{d}{os.pathsep}{os.environ['PATH']}"
            r = run_settings_command(bash("dvc repro examples/vera/dvc.yaml"), path)
            self.assertEqual((r.returncode, r.stderr), (0, ""))
            self.assertEqual(
                run_settings_command(bash(BLOCKED[0]), path).returncode, 99
            )


if __name__ == "__main__":
    unittest.main()
