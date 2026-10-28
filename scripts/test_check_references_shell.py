"""check-references-shell.sh の検査（一時の git のリポジトリの最小の flake で、nix eval を本物で実行する）。

実行: python3 -m unittest discover -s scripts -p 'test_check_references_shell.py'
"""

import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parent / "check-references-shell.sh"

# nixpkgs を使わない derivation。nix eval は store のパスを計算するだけで、ビルドしない
FLAKE = """\
{{
  outputs = {{ self }}: {{
    devShells.x86_64-linux.references = derivation {{
      name = "nix-shell";
      system = "x86_64-linux";
      builder = "/bin/sh";
      args = [ "-c" "echo {tag} > $out" ];
    }};
  }};
}}
"""


def run(*args, cwd):
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True, check=False)


class CheckReferencesShell(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name)
        run("git", "init", "-q", cwd=self.repo).check_returncode()
        self.write_flake("a")

    def tearDown(self):
        self.tmp.cleanup()

    def write_flake(self, tag):
        (self.repo / "flake.nix").write_text(FLAKE.format(tag=tag))
        run("git", "add", "flake.nix", cwd=self.repo).check_returncode()

    def shell_path(self):
        result = run(
            "nix",
            "eval",
            "--raw",
            ".#devShells.x86_64-linux.references.outPath",
            cwd=self.repo,
        )
        result.check_returncode()
        return result.stdout

    def record(self, rel, path):
        file = self.repo / rel
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text(path)
        run("git", "add", rel, cwd=self.repo).check_returncode()

    def check(self):
        return run("bash", str(SCRIPT), cwd=self.repo)

    def test_passes_when_every_recorded_path_matches_the_shell(self):
        path = self.shell_path()
        self.record("examples/a/references-shell.txt", path)
        self.record("examples/a/b/references-shell.txt", path)
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_fails_after_the_flake_changes_the_shell(self):
        old = self.shell_path()
        self.record("examples/a/references-shell.txt", old)
        self.write_flake("b")
        new = self.shell_path()
        self.assertNotEqual(old, new)
        result = self.check()
        self.assertEqual(result.returncode, 1)
        self.assertIn("examples/a/references-shell.txt", result.stderr)
        self.assertIn(old, result.stderr)
        self.assertIn(new, result.stderr)
        self.assertIn("ADR-0212", result.stderr)

    def test_names_only_the_mismatched_files(self):
        path = self.shell_path()
        self.record("examples/a/references-shell.txt", path)
        self.record("examples/b/references-shell.txt", "/nix/store/old-nix-shell")
        result = self.check()
        self.assertEqual(result.returncode, 1)
        self.assertIn("examples/b/references-shell.txt", result.stderr)
        self.assertNotIn("examples/a/references-shell.txt", result.stderr)

    def test_ignores_untracked_records(self):
        path = self.shell_path()
        self.record("examples/a/references-shell.txt", path)
        stray = self.repo / "examples/c/references-shell.txt"
        stray.parent.mkdir(parents=True)
        stray.write_text("/nix/store/old-nix-shell")
        result = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_names_a_tracked_record_missing_from_the_work_tree(self):
        path = self.shell_path()
        self.record("examples/a/references-shell.txt", path)
        self.record("examples/b/references-shell.txt", path)
        (self.repo / "examples/b/references-shell.txt").unlink()
        result = self.check()
        self.assertEqual(result.returncode, 1)
        self.assertIn("examples/b/references-shell.txt", result.stderr)
        self.assertIn("ADR-0212", result.stderr)

    def test_fails_without_tracked_records(self):
        result = self.check()
        self.assertEqual(result.returncode, 2)
        self.assertIn("references-shell.txt", result.stderr)


if __name__ == "__main__":
    unittest.main()
