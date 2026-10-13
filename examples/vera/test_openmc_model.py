"""openmc_model.py の `model` が書く model.xml の性質（ADR-0136）。

OpenMC のある環境で実行する: nix develop .#references -c python3 -m unittest discover -s examples/vera -p 'test_openmc_model.py'
"""

import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "openmc_model.py"
INPUT = HERE / "p1" / "openmc-1B.toml"


def write_model(source: Path, directory: Path) -> bytes:
    subprocess.run(
        [sys.executable, str(SCRIPT), "model", str(source), str(directory)],
        check=True,
    )
    return (directory / "model.xml").read_bytes()


class ModelXml(unittest.TestCase):
    def test_does_not_depend_on_the_library_path_or_the_thread_count(self):
        """ライブラリのパスとスレッド数は、model.xml に入らない（入ると、実行の環境でモデルが変わって見える）。"""
        text = INPUT.read_text()
        other = re.sub(r"(?m)^threads = .*$", "threads = 3", text)
        other = re.sub(
            r"(?m)^cross-sections = .*$", 'cross-sections = "elsewhere.xml"', other
        )
        self.assertNotEqual(text, other)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "other.toml").write_text(other)
            a = write_model(INPUT, tmp / "a")
            b = write_model(tmp / "other.toml", tmp / "b")
        self.assertEqual(a, b)
        self.assertNotIn(b"cross_sections", a)

    def test_is_the_same_bytes_every_time(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self.assertEqual(
                write_model(INPUT, tmp / "a"), write_model(INPUT, tmp / "b")
            )


if __name__ == "__main__":
    unittest.main()
