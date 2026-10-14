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
            "diff_keno_pcm": 193.0,
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


def p5_summary(problem: str, kind: str, k: float, ref: float, **extra) -> dict:
    """compare_p5.py の要約（表に使う項目だけ）"""
    s = {
        "format": "xlc-compare-p5",
        "kind": kind,
        "problem": problem,
        "k": {
            "xlc": k,
            "reference": ref,
            "reference_std": 8e-06,
            "diff_pcm": 1e5 * (k / ref - 1),
            **extra,
        },
        "assemblies": {
            "relative": {
                "n": 193,
                "rms": 0.0166,
                "max": 0.0463,
                "at": [1, 11],
                "name": "D02",
            },
            "peripheral": {
                "n": 56,
                "rms": 0.0191,
                "max": 0.0463,
                "at": [1, 11],
                "name": "D02",
            },
        },
    }
    if kind == "openmc":
        s["pins"] = {
            "all": {"n": 50952, "rms": 0.0317, "max": 0.5473, "at": [3, 68]},
            "peripheral": {"n": 14784, "rms": 0.0535, "max": 0.5473, "at": [3, 68]},
            "max_pin": {"xlc": 1.48451, "reference": 1.43999},
        }
    else:
        s["max_pin"] = {"xlc": 1.48451, "reference": 1.44402}
    return s


class DrawP5(unittest.TestCase):
    def test_keno_comparison_row_has_rod_worth_from_5a(self):
        """制御棒の価値は (1/k_5B − 1/k_5A) × 1e5 pcm。5A の行は空欄"""
        sources = [
            p5_summary("5A", "keno", 1.00540075, 1.004085),
            p5_summary("5B", "keno", 0.992924, 0.991496),
        ]
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in vt.p5_keno_comparison(sources).split("\n")[2:]
        ]
        self.assertEqual(
            rows[0],
            [
                "5A-2D",
                "1.005401",
                "1.004085",
                "+131.0",
                "",
                "1.66%、+4.6%（D02）",
                "1.485、1.444",
            ],
        )
        self.assertEqual(rows[1][4], "1250、1265 pcm")

    def test_openmc_comparison_names_the_assembly(self):
        sources = [
            p5_summary("5A", "openmc", 1.005401, 1.005459, sigma=2.9),
        ]
        row = [
            c.strip()
            for c in vt.p5_openmc_comparison(sources)
            .split("\n")[2]
            .strip("|")
            .split("|")
        ]
        self.assertEqual(
            row,
            [
                "5A-2D",
                "\u22125.8",
                "1.66%、+4.6%（D02）",
                "3.17%、+55%",
                "5.35%",
                "1.485、1.440",
            ],
        )

    def test_reference_results_leave_missing_openmc_blank(self):
        """OpenMC が 5A だけのとき、5B と 5C の OpenMC の列と制御棒の価値（OpenMC）は空欄。k の差は JSON の絶対差 ± σ"""
        sources = [
            p5_summary("5A", "keno", 1.005401, 1.004085),
            p5_summary("5B", "keno", 0.992924, 0.991496),
            p5_summary(
                "5A",
                "openmc",
                1.005401,
                1.005459,
                reference_std=2e-05,
                diff_keno_pcm=136.84,
                diff_keno_std_pcm=2.15,
            ),
        ]
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in vt.p5_reference_results(sources).split("\n")[2:]
        ]
        self.assertEqual(
            rows[0][:4],
            ["5A-2D", "1.005459 ± 0.000020", "1.004085 ± 0.000008", "+137 ± 2 pcm"],
        )
        self.assertEqual(rows[1][1], "")
        self.assertEqual(rows[1][5], "1265 ± 1 pcm")

def lattice_summary(problem: str, diff: float, pins: bool = True) -> dict:
    return {
        "format": "xlc-compare-keno-lattice",
        "problem": problem,
        "k": {
            "xlc": 1.18,
            "keno": 1.18336,
            "keno_std": 2.4e-05,
            "diff_pcm": diff,
        },
        "pins": {"positions": 39, "rms": 0.0005, "max": -0.00123, "at": [1, 0]}
        if pins
        else None,
    }


class Endfb70(unittest.TestCase):
    def test_lattice_row_puts_the_two_libraries_side_by_side(self):
        table = vt.endfb70_lattices(
            [lattice_summary("2B", -336.0), lattice_summary("2B", 45.4)]
        ).split("\n")
        row = [c.strip() for c in table[2].strip("|").split("|")]
        self.assertEqual(
            row,
            [
                "2B",
                "1.180000",
                "1.183360 ± 0.000024",
                "\u2212336 pcm",
                "+45 pcm",
                "0.050%、\u22120.123%",
                "0.050%、\u22120.123%",
            ],
        )

    def test_a_pin_cell_has_no_pin_columns(self):
        table = vt.endfb70_lattices(
            [lattice_summary("1B", 1.0, False), lattice_summary("1B", 2.0, False)]
        ).split("\n")
        row = [c.strip() for c in table[2].strip("|").split("|")]
        self.assertEqual(row[5:], ["", ""])

    def test_the_two_libraries_must_list_the_same_problems(self):
        with self.assertRaises(vt.TableError):
            vt.endfb70_lattices(
                [lattice_summary("2B", 1.0), lattice_summary("2C", 1.0)]
            )

    def test_p3_k_row(self):
        table = vt.endfb70_p3_k([summary("3B"), summary("3B")]).split("\n")
        row = [c.strip() for c in table[2].strip("|").split("|")]
        self.assertEqual(row[0], "3B")
        self.assertEqual(len(row), 7)


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

    def test_malformed_and_stray_marks_are_errors(self):
        """目印の書き損じを文章として読み飛ばすと、手で直した表が検査を通る"""
        self.lock()
        for text in [
            "<!-- table: demo  -->\n<!-- /table -->\n",
            "<!--table: demo-->\n<!-- /table -->\n",
            " <!-- table: demo -->\n<!-- /table -->\n",
            "<!-- /table -->\n",
            "<!-- table: demo -->\n<!-- /table  -->\n",
        ]:
            with self.assertRaises(vt.TableError, msg=text):
                vt.render(text, self.root)

    def test_an_invalid_or_old_json_is_an_error_not_a_traceback(self):
        (self.root / "t/a.json").write_text('{"format": "xlc-compare-p3"')
        self.lock()
        with self.assertRaisesRegex(vt.TableError, "not valid JSON"):
            vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)
        (self.root / "t/a.json").write_text(
            '{"format": "xlc-compare-p3", "k": {"openmc": 1}}'
        )
        self.lock()
        with self.assertRaisesRegex(vt.TableError, "lacks"):
            vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)

    def test_a_value_rounding_to_zero_has_no_sign(self):
        self.assertEqual(vt.signed(-0.3, 0), "0")
        self.assertEqual(vt.signed(-0.0004, 3), "0.000")
        self.assertEqual(vt.signed(-1.5, 1), "\u22121.5")

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
