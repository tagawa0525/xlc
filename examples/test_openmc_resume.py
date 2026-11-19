"""openmc_resume.py の `run_model` が、実行の間に statepoint を最新の `KEEP` 個に保つこと（#0106）と、`check_model` が
作業ディレクトリの model.xml と入力のモデルの違いを見つけること（#0104）。

OpenMC の Python を import するので、python-tests の devShell で実行する（scripts/run-python-tests.sh）。
"""

import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

import openmc
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


def sphere(density: float, library: str, statepoints: list[int] | None) -> openmc.Model:
    """水の球の小さなモデル。ID は毎回 1 から振る（run と extract は別のプロセスで、どちらも 1 から振るため）。"""
    openmc.reset_auto_ids()
    water = openmc.Material()
    water.add_nuclide("H1", 2.0)
    water.add_nuclide("O16", 1.0)
    water.set_density("g/cm3", density)
    materials = openmc.Materials([water])
    materials.cross_sections = library
    cell = openmc.Cell(
        fill=water, region=-openmc.Sphere(r=10.0, boundary_type="vacuum")
    )
    settings = openmc.Settings()
    settings.particles = 100
    settings.batches = 10
    settings.inactive = 2
    if statepoints is not None:
        settings.statepoint = {"batches": statepoints}
    return openmc.Model(openmc.Geometry([cell]), materials, settings)


class CheckModel(unittest.TestCase):
    def ran(self, work: Path) -> None:
        """run が `work` に書いた model.xml（別の worktree のライブラリ、`--statepoint-interval` で変えた間隔）。"""
        work.mkdir()
        sphere(1.0, "/elsewhere/cross_sections.xml", [5, 10]).export_to_model_xml(
            str(work / "model.xml")
        )

    def test_accepts_the_same_model_run_with_another_library_path_and_interval(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp) / "work"
            self.ran(work)
            openmc_resume.check_model(
                sphere(1.0, "/here/cross_sections.xml", None), work
            )

    def test_stops_when_the_model_differs(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp) / "work"
            self.ran(work)
            with self.assertRaises(SystemExit) as stop:
                openmc_resume.check_model(
                    sphere(0.7, "/elsewhere/cross_sections.xml", [5, 10]), work
                )
        self.assertIn("model.xml", str(stop.exception))

    def test_stops_when_the_work_has_no_model_xml(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(SystemExit):
            openmc_resume.check_model(
                sphere(1.0, "/here/cross_sections.xml", None), Path(tmp)
            )


if __name__ == "__main__":
    unittest.main()
