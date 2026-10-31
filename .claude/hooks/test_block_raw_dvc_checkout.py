"""block-raw-dvc-checkout.py の検査（ADR-0223）。

実行: cd .claude/hooks && python3 -m unittest test_block_raw_dvc_checkout
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

HOOK = Path(__file__).resolve().parent / "block-raw-dvc-checkout.py"
SPEC = importlib.util.spec_from_file_location("block_raw_dvc_checkout", HOOK)
hook = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hook)

BLOCKED = [
    "dvc checkout",
    "dvc checkout examples/vera/dvc.yaml",
    "nix develop -c dvc checkout dvc.yaml",
    "nix develop .#references --command dvc checkout dvc.yaml",
    "cd /home/tagawa/github/xlc-x && dvc checkout",
    "cd /home/tagawa/github/xlc-x; dvc checkout dvc.yaml",
    "dvc status\ndvc checkout",
    "git status || dvc checkout",
    "dvc checkout 2>&1 | tail -5",
    "dvc -q checkout",
    "dvc --cd examples/vera checkout dvc.yaml",
    "/nix/store/abc-dvc/bin/dvc checkout",
    "FOO=1 dvc checkout",
    "env FOO=1 dvc checkout",
    "time dvc checkout",
    "taskset -c 0 dvc checkout",
    "timeout 60 dvc checkout",
    "bash -c 'dvc checkout dvc.yaml'",
    "nix develop -c bash -lc 'cd x && dvc checkout'",
    "python3 -m dvc checkout",
    "(cd x && dvc checkout)",
    "echo $(dvc checkout)",
    "echo ok # comment\ndvc checkout",
    "cat <<'EOF' > f\nhello\nEOF\ndvc checkout",
    "dvc \\\n  checkout",
    "for y in examples/*/dvc.yaml; do dvc checkout $y; done",
    "if dvc checkout; then echo ok; fi",
    "if true; then echo ok; else dvc checkout; fi",
    "while false; do :; done; { dvc checkout; } 2>&1 | tail",
    "! dvc checkout",
    "bash -o pipefail -c 'dvc checkout'",
    "bash -euo pipefail -c 'dvc checkout dvc.yaml'",
    "bash --norc -c 'dvc checkout'",
    "while read l; do echo $l; done < <(dvc checkout 2>&1)",
    "grep -c failed <(dvc checkout 2>&1)",
    "nix develop --command bash -c 'dvc checkout'",
    "nix develop .#references --command bash -c 'cd x && dvc checkout'",
    "nix develop -c -- dvc checkout",
    "env -iC /x dvc checkout",
    "echo $((1<<2))\ndvc checkout",
    "((x = 1 << 2))\ndvc checkout",
    "echo $(( (1 + 2) << 1 )); dvc checkout",
    "echo $'don\\'t'; dvc checkout",
]

ALLOWED = [
    "scripts/dvc-checkout.py",
    "scripts/dvc-checkout.py examples/vera/dvc.yaml",
    "nix develop -c scripts/dvc-checkout.py dvc.yaml",
    "cd /home/tagawa/github/xlc-x && nix develop -c scripts/dvc-checkout.py dvc.yaml examples/vera/dvc.yaml",
    "python3 scripts/dvc-checkout.py",
    "scripts/worktree-restore-outputs.sh",
    "grep 'dvc checkout' -r .claude/skills",
    'grep -rn "dvc checkout" docs',
    "rg dvc\\ checkout",
    "echo dvc checkout",
    "git log --grep='dvc checkout'",
    "dvc status",
    "dvc repro examples/vera/dvc.yaml",
    "dvc commit examples/vera/dvc.yaml",
    "git checkout main",
    "# dvc checkout",
    "echo ok # dvc checkout",
    "git commit -m \"$(cat <<'EOF'\ndocs: replace dvc checkout\nEOF\n)\"",
    "cat > f <<'EOF'\ndvc checkout\nEOF",
    "cat > f <<-EOF\n\tdvc checkout\n\tEOF",
    "bash -o pipefail -c 'echo dvc checkout'",
    "echo $'dvc checkout\\n'",
    "for f in a b; do echo dvc checkout $f; done",
    "echo $((1 + 2)) dvc",
    "nix develop",
    "dvc",
    "",
]


class RawDvcCheckoutTest(unittest.TestCase):
    def test_commands_that_run_dvc_checkout_directly_are_blocked(self):
        for command in BLOCKED:
            with self.subTest(command=command):
                self.assertTrue(hook.runs_raw_dvc_checkout(command))

    def test_script_calls_and_quoted_text_are_allowed(self):
        for command in ALLOWED:
            with self.subTest(command=command):
                self.assertFalse(hook.runs_raw_dvc_checkout(command))

    def test_unbalanced_quotes_are_an_error(self):
        with self.assertRaises(ValueError):
            hook.runs_raw_dvc_checkout("echo 'dvc checkout")


def run_hook(payload):
    return subprocess.run(
        [sys.executable, str(HOOK)],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        check=False,
    )


class HookProcessTest(unittest.TestCase):
    def test_blocked_command_exits_2_with_the_script_to_use(self):
        r = run_hook(
            {"tool_name": "Bash", "tool_input": {"command": "dvc checkout dvc.yaml"}}
        )
        self.assertEqual(r.returncode, 2)
        self.assertIn("scripts/dvc-checkout.py", r.stderr)
        self.assertIn("255", r.stderr)

    def test_allowed_command_exits_0_silently(self):
        r = run_hook(
            {
                "tool_name": "Bash",
                "tool_input": {
                    "command": "nix develop -c scripts/dvc-checkout.py dvc.yaml"
                },
            }
        )
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_other_tools_pass(self):
        r = run_hook({"tool_name": "Read", "tool_input": {"file_path": "/x"}})
        self.assertEqual((r.returncode, r.stderr), (0, ""))

    def test_unparsable_command_is_reported_without_blocking(self):
        r = run_hook(
            {"tool_name": "Bash", "tool_input": {"command": "echo 'dvc checkout"}}
        )
        self.assertEqual(r.returncode, 1)
        self.assertIn("読めない", r.stderr)

    def test_non_string_command_fails_without_blocking(self):
        r = run_hook({"tool_name": "Bash", "tool_input": {"command": None}})
        self.assertEqual(r.returncode, 1)
        self.assertIn("hook の入力を読めない", r.stderr)
        self.assertNotIn("Traceback", r.stderr)

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
# PATH を絞るテストでも同じシェルで実行するため、絶対パスで持つ
SH = shutil.which("sh")


def run_settings_command(payload, project_dir, path=None):
    """.claude/settings.json の PreToolUse の Bash の command を、Claude Code と同じくシェルで実行する。"""
    (entry,) = json.loads(SETTINGS.read_text())["hooks"]["PreToolUse"]
    (handler,) = entry["hooks"]
    env = {k: v for k, v in os.environ.items() if k != "CLAUDE_PROJECT_DIR"}
    if project_dir is not None:
        env["CLAUDE_PROJECT_DIR"] = project_dir
    if path is not None:
        env["PATH"] = path
    return subprocess.run(
        [SH, "-c", handler["command"]],
        input=json.dumps(payload),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )


BLOCKED_PAYLOAD = {"tool_name": "Bash", "tool_input": {"command": "dvc checkout"}}


class SettingsCommandTest(unittest.TestCase):
    def test_settings_command_runs_the_hook(self):
        r = run_settings_command(BLOCKED_PAYLOAD, str(HOOK.parent.parent.parent))
        self.assertEqual(r.returncode, 2)

    def test_settings_command_blocks_usual_forms(self):
        for command in (
            "nix develop -c dvc checkout dvc.yaml",
            "cd /home/tagawa/github/xlc-x && dvc checkout",
            "bash -c 'dvc checkout'",
            "for y in examples/*/dvc.yaml; do dvc checkout $y; done",
        ):
            with self.subTest(command=command):
                payload = {"tool_name": "Bash", "tool_input": {"command": command}}
                self.assertEqual(run_settings_command(payload, PROJECT).returncode, 2)

    def test_commands_mentioning_check_reach_python(self):
        with tempfile.TemporaryDirectory() as d:
            fake = Path(d) / "python3"
            fake.write_text("#!/bin/sh\nexit 99\n")
            fake.chmod(0o755)
            path = f"{d}{os.pathsep}{os.environ['PATH']}"
            payload = {"tool_name": "Bash", "tool_input": {"command": "dvc 'check'out"}}
            self.assertEqual(
                run_settings_command(payload, PROJECT, path).returncode, 99
            )

    def test_missing_python_reports_without_blocking(self):
        with tempfile.TemporaryDirectory() as d:
            for tool in ("sh", "cat"):
                os.symlink(shutil.which(tool), Path(d) / tool)
            r = run_settings_command(BLOCKED_PAYLOAD, PROJECT, d)
            self.assertEqual(r.returncode, 127)
            self.assertNotEqual(r.stderr, "")

    def test_commands_without_checkout_skip_python(self):
        with tempfile.TemporaryDirectory() as d:
            fake = Path(d) / "python3"
            fake.write_text("#!/bin/sh\nexit 99\n")
            fake.chmod(0o755)
            path = f"{d}{os.pathsep}{os.environ['PATH']}"
            project = str(HOOK.parent.parent.parent)
            payload = {"tool_name": "Bash", "tool_input": {"command": "dvc status"}}
            r = run_settings_command(payload, project, path)
            self.assertEqual((r.returncode, r.stderr), (0, ""))
            r = run_settings_command(BLOCKED_PAYLOAD, project, path)
            self.assertEqual(r.returncode, 99)

    def test_missing_project_dir_reports_without_blocking(self):
        for project_dir in (None, "", "/nonexistent"):
            with self.subTest(project_dir=project_dir):
                r = run_settings_command(BLOCKED_PAYLOAD, project_dir)
                self.assertEqual(r.returncode, 1)
                self.assertNotEqual(r.stderr, "")


if __name__ == "__main__":
    unittest.main()
