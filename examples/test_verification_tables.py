"""verification_tables.py の検査（表の描き方、目印の置き換え、--check の判定）。

実行: python3 -m unittest discover -s examples -p 'test_verification_tables.py'
"""

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verification_tables as vt


def summary(problem: str, **k) -> dict:
    return {
        "format": "xlc-compare-p3",
        "problem": problem,
        "k": {
            "xlc": 1.0020840,
            "keno": 1.000154,
            "keno_std": 6e-06,
            "openmc": 1.0013894,
            "openmc_std": 1.9e-05,
            "diff_openmc_pcm": 69.46,
            "diff_openmc_std_pcm": 1.9,
            "openmc_keno_diff_pcm": 123.54,
            **k,
        },
        "axial_openmc_vs_keno": {"ao_subject": -0.0021075, "ao_reference": -0.00062},
        "axial_vs_openmc": {"rms": 0.0030664, "relative_rms": 0.0095799},
        "pins_vs_openmc": {"rms": 0.00071399, "max": 0.00191603},
        "layer_pins": {"max": -0.0977339, "max_level": 1, "levels": [{}] * 49},
    }


class Draw(unittest.TestCase):
    def test_reference_results_row(self):
        table = vt.p3_reference_results([summary("3B")]).split("\n")
        self.assertEqual(len(table), 3)
        row = [c.strip() for c in table[2].strip("|").split("|")]
        self.assertEqual(
            row,
            [
                "3B",
                "1.001389 ± 0.000019",
                "1.000154 ± 0.000006",
                "+123.5 pcm",
                "\u22120.211%",
                "\u22120.062%",
            ],
        )

    def test_openmc_comparison_row_names_the_layer(self):
        row = [
            c.strip()
            for c in vt.p3_openmc_comparison([summary("3B")])
            .split("\n")[2]
            .strip("|")
            .split("|")
        ]
        self.assertEqual(
            row,
            [
                "3B",
                "1.002084",
                "1.001389",
                "+69 ± 2 pcm",
                "0.0031（0.96%）",
                "0.071%、+0.192%",
                "\u22129.8%（下端の層）",
            ],
        )

    def test_columns_are_padded_by_display_width(self):
        """全角は幅 2。どの行も同じ幅に並ぶ（prettier と同じ）"""
        table = vt.markdown_table(["問題", "k"], [["3A", "1.0"], ["3B 長い", "1"]])
        widths = {vt.width(line) for line in table.split("\n")}
        self.assertEqual(len(widths), 1)

    def test_a_summary_without_openmc_is_an_error(self):
        s = summary("3A")
        del s["k"]["openmc"]
        with self.assertRaises(vt.TableError):
            vt.p3_reference_results([s])


class Replace(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)
        (self.root / "t").mkdir()
        (self.root / "t/a.json").write_text(json.dumps(summary("3A")))
        self.saved = vt.TABLES.copy()
        vt.TABLES.clear()
        vt.TABLES["demo"] = (["t/a.json"], lambda s: vt.p3_openmc_comparison(s))
        self.addCleanup(lambda: (vt.TABLES.clear(), vt.TABLES.update(self.saved)))

    def lock(self, md5: str | None = None):
        md5 = md5 or hashlib.md5((self.root / "t/a.json").read_bytes()).hexdigest()
        (self.root / "t/dvc.lock").write_text(
            f"schema: '2.0'\nstages:\n  s:\n    outs:\n    - path: a.json\n      hash: md5\n      md5: {md5}\n"
        )

    def test_replaces_only_between_the_marks_and_is_idempotent(self):
        self.lock()
        text = "before\n<!-- table: demo -->\n\nold\n\nprose in the block\n<!-- /table -->\nafter\n"
        new = vt.render(text, self.root)
        self.assertTrue(new.startswith("before\n<!-- table: demo -->\n\n| 問題"))
        self.assertTrue(new.endswith("\n\n<!-- /table -->\nafter\n"))
        self.assertNotIn("old", new)
        self.assertEqual(vt.render(new, self.root), new)

    def test_unknown_unclosed_and_nested_marks_are_errors(self):
        self.lock()
        for text in [
            "<!-- table: nope -->\n<!-- /table -->\n",
            "<!-- table: demo -->\n",
            "<!-- table: demo -->\n<!-- table: demo -->\n<!-- /table -->\n",
            "<!-- table: demo -->\n<!-- /table -->\n<!-- table: demo -->\n<!-- /table -->\n",
        ]:
            with self.assertRaises(vt.TableError, msg=text):
                vt.render(text, self.root)

    def test_json_differing_from_dvc_lock_is_an_error(self):
        self.lock("0" * 32)
        with self.assertRaises(vt.TableError):
            vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)

    def test_missing_json_is_an_error_naming_dvc_repro(self):
        (self.root / "t/a.json").unlink()
        with self.assertRaisesRegex(vt.TableError, "dvc repro"):
            vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)

    def test_check_fails_on_a_hand_edited_table(self):
        self.lock()
        record = self.root / "record.md"
        record.write_text("<!-- table: demo -->\n<!-- /table -->\n")
        vt.ROOT, saved = self.root, vt.ROOT
        self.addCleanup(lambda: setattr(vt, "ROOT", saved))
        self.assertEqual(vt.main([str(record)]), 0)
        self.assertEqual(vt.main(["--check", str(record)]), 0)
        record.write_text(record.read_text().replace("1.002084", "1.002085"))
        self.assertEqual(vt.main(["--check", str(record)]), 1)


if __name__ == "__main__":
    unittest.main()
