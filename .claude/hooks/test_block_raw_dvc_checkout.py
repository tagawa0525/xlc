"""block-raw-dvc-checkout.py の検査（ADR-0223）。

実行: cd .claude/hooks && python3 -m unittest test_block_raw_dvc_checkout
"""

import importlib.util
import json
import subprocess
import sys
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


if __name__ == "__main__":
    unittest.main()
