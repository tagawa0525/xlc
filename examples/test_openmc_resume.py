"""openmc_resume.py の `run_model` が、実行の間に statepoint を最新の `KEEP` 個に保つこと（#0106）。

OpenMC の Python を import するので、python-tests の devShell で実行する（scripts/run-python-tests.sh）。
"""

import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

import openmc_resume


class ChdirModel:
    """`openmc.Model.run(cwd=...)` と同じく作業ディレクトリを `cwd` に移し、statepoint を世代の順に書く。"""

    def __init__(self, generations: int, pause: float) -> None:
        self.generations = generations
        self.pause = pause
        self.left_at_end: list[str] = []

    def run(self, cwd, threads, output, restart_file):
        back = Path.cwd()
        os.chdir(cwd)
        try:
            for n in range(1, self.generations + 1):
                Path(f"statepoint.{n:03d}.h5").write_bytes(b"")
                time.sleep(self.pause)
            self.left_at_end = sorted(p.name for p in Path().glob("statepoint.*.h5"))
        finally:
            os.chdir(back)


class RunModel(unittest.TestCase):
    def test_keeps_only_the_latest_statepoints_while_running_with_a_relative_work(self):
        """`work` が相対パス（run_remote.sh の `results/<名前>/run`）でも、実行の間に古い statepoint を消す。"""
        model = ChdirModel(generations=6, pause=0.2)
        with tempfile.TemporaryDirectory() as tmp:
            back = Path.cwd()
            os.chdir(tmp)
            try:
                work = Path("results/case/run")
                work.mkdir(parents=True)
                with mock.patch.object(openmc_resume, "POLL_SECONDS", 0.01):
                    openmc_resume.run_model(
                        model, {"threads": 1}, work, openmc_resume.Options()
                    )
            finally:
                os.chdir(back)
        self.assertEqual(model.left_at_end, ["statepoint.005.h5", "statepoint.006.h5"])


if __name__ == "__main__":
    unittest.main()
