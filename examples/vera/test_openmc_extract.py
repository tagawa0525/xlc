"""手で流す参照解（openmc_core_2d.py、openmc_model_3d.py）の extract が、作業ディレクトリの model.xml と入力のモデルの
一致を確かめること（#0104）。

OpenMC の Python を import するので、python-tests の devShell で実行する（scripts/run-python-tests.sh）。
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import openmc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import openmc_core_2d
import openmc_model_3d


def write_model_xml(module, source: Path, work: Path) -> None:
    """`module` の run が `work` に書く model.xml を、OpenMC を実行せずに書く。"""

    def export(model, run, work, options):
        model.export_to_model_xml(str(work / "model.xml"))

    openmc.reset_auto_ids()
    with mock.patch.object(module, "run_model", export):
        module.solve(source, work, module.Options())


# (スクリプト, model.xml を書いた入力, extract に渡す入力)
CASES = (
    (openmc_core_2d, "p5/openmc-5B-2D.toml", "p5/openmc-5A-2D.toml"),
    (openmc_model_3d, "p3/openmc-3B.toml", "p3/openmc-3A.toml"),
)


class ExtractChecksTheModel(unittest.TestCase):
    def test_stops_when_the_work_model_was_built_from_another_input(self):
        """statepoint を出した model.xml が入力のモデルと違えば、statepoint を読む前に止まる。"""
        for module, ran, given in CASES:
            with self.subTest(given), tempfile.TemporaryDirectory() as tmp:
                work = Path(tmp) / "work"
                write_model_xml(module, HERE / ran, work)
                openmc.reset_auto_ids()
                with self.assertRaises(SystemExit) as stop:
                    module.extract(HERE / given, work, Path(tmp) / "results")
                self.assertIn("model.xml", str(stop.exception))


if __name__ == "__main__":
    unittest.main()
