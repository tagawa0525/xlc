"""verification_tables.py の検査（表の描き方、目印の置き換え、--check の判定）。

実行: python3 -m unittest discover -s examples -p 'test_verification_tables.py'
"""

import hashlib
import json
import re
import subprocess
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


def openmc_summary(problem: str, pins: bool) -> dict:
    pin = {"positions": 39, "rms": 0.0006, "max": 0.0015, "at": [1, 0]}
    return {
        "format": "xlc-compare-openmc-keno",
        "problem": problem,
        "k": {
            "xlc": 1.184112,
            "openmc": 1.18300,
            "openmc_std": 4e-05,
            "keno": 1.182149,
            "keno_std": 6.8e-05,
            "diff_xlc_openmc_pcm": 111.2,
            "diff_xlc_openmc_std_pcm": 4.0,
            "diff_openmc_keno_pcm": 85.1,
            "diff_openmc_keno_std_pcm": 7.9,
            "diff_xlc_keno_pcm": 196.3,
        },
        "pins_openmc_vs_keno": pin if pins else None,
        "pins_xlc_vs_openmc": {**pin, "positions": 264} if pins else None,
    }


class OpenmcEndfb70(unittest.TestCase):
    def test_row_has_the_three_k_differences_with_errors(self):
        table = vt.endfb70_openmc([openmc_summary("1B", False)]).split("\n")
        row = [c.strip() for c in table[2].strip("|").split("|")]
        self.assertEqual(
            row,
            [
                "1B",
                "1.184112",
                "1.183000 ± 0.000040",
                "1.182149 ± 0.000068",
                "+111 ± 4 pcm",
                "+85 ± 8 pcm",
                "+196 pcm",
                "",
                "",
            ],
        )

    def test_assembly_row_has_pin_columns(self):
        table = vt.endfb70_openmc([openmc_summary("2B", True)]).split("\n")
        row = [c.strip() for c in table[2].strip("|").split("|")]
        # 比べる点の数が違う（KENO-VI の 1/8 の表と、燃料ピンすべて）ので、数を添える
        self.assertEqual(
            row[7:], ["0.060%、+0.150%（39 点）", "0.060%、+0.150%（264 点）"]
        )


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

    def git(self, *args: str):
        subprocess.run(
            ["git", "-C", str(self.root), *args], check=True, capture_output=True
        )

    def test_an_untracked_json_is_an_error_naming_cache_false(self):
        self.lock()
        self.git("init", "-q", "-b", "main")
        with self.assertRaisesRegex(vt.TableError, "cache: false"):
            vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)

    def test_a_tracked_json_passes_in_a_git_repository(self):
        self.lock()
        self.git("init", "-q", "-b", "main")
        self.git("add", "t/a.json")
        vt.render("<!-- table: demo -->\n<!-- /table -->\n", self.root)

    def test_the_git_tracked_json_must_be_tracked_too(self):
        source = "examples/depletion/aws-check/compare-aws-1C-rings1.json"
        (self.root / source).parent.mkdir(parents=True)
        (self.root / source).write_text(json.dumps(summary("3A")))
        self.git("init", "-q", "-b", "main")
        with self.assertRaisesRegex(vt.TableError, "cache: false"):
            vt.load(self.root, source, True)

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


def depletion(changes: list[float], rings: int = 1) -> dict:
    """compare_depletion.py の要約（2 点、U-235 だけ OpenMC との差がある）"""
    return {
        "format": "xlc-compare-depletion",
        "burnup": [0.0, 1.0],
        "transported": 64,
        "k": [
            {
                "xlc": 1.17,
                "openmc": 1.171,
                "openmc_std_pcm": 20.0,
                "diff_pcm": -100.0 + c,
                "change_pcm": c,
            }
            for c in changes
        ],
        "pin": {
            "U-235": [0.0, 0.0123],
            "U-238": [0.0, -0.0001],
            "Pu-239": [None, 0.02],
            "Xe-135": [None, -0.005],
            "Sm-149": [None, 0.0],
        },
        "seconds": {
            "xlc": {
                "shielding": 90.1,
                "moc": 3.52,
                "rates": 0.018,
                "matrix": 0.186,
                "cram": 0.283,
                "solves": 13,
            },
            "openmc": {
                "transport": 5041.6,
                "fastest": 383.8,
                "slowest": 394.7,
                "solves": 13,
                "cram": 5.18,
                "total": 5047.0,
            },
        },
        "rings": [
            {
                "U-235": 0.01,
                "U-238": 0.0,
                "Pu-239": 0.03,
                "Xe-135": None,
                "Sm-149": -0.001,
            }
        ]
        * rings,
    }


class DrawDepletion(unittest.TestCase):
    def test_k_rows_list_each_case_by_burnup(self):
        table = vt.depletion_k(
            [depletion([0.0, -30.0]), depletion([0.0, 12.0], rings=4)]
        )
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(
            rows[1],
            ["1", "\u2212130 ± 20", "\u221230", "\u221288 ± 20", "+12"],
        )

    def test_density_rows_in_percent_with_a_dash_for_no_reference(self):
        table = vt.depletion_densities([depletion([0.0, 0.0])])
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(rows[0], ["0", "0.00", "0.00", "–", "–", "–"])
        self.assertEqual(
            rows[1], ["1", "+1.23", "\u22120.01", "+2.00", "\u22120.50", "0.00"]
        )

    def test_ring_rows_number_the_rings_from_the_inside(self):
        table = vt.depletion_rings([depletion([0.0, 0.0], rings=2)])
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(
            rows[1], ["2（外）", "+1.00", "0.00", "+3.00", "–", "\u22120.10"]
        )

    def test_time_rows_put_xlc_and_openmc_of_each_case(self):
        table = vt.depletion_times(
            [depletion([0.0, 0.0]), depletion([0.0, 0.0], rings=4)]
        )
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(len(rows), 4)
        self.assertEqual(
            rows[0],
            ["xlc、1 リング", "13", "94.1", "90.1", "3.5", "0.02", "0.19", "0.28"],
        )
        self.assertEqual(
            rows[1],
            [
                "OpenMC、1 リング",
                "13",
                "5,047",
                "5,042（1 回 384〜395）",
                "",
                "",
                "",
                "5.18",
            ],
        )
        self.assertEqual(rows[2][0], "xlc、4 リング")

    def test_cases_with_other_burnup_points_are_an_error(self):
        other = depletion([0.0, 0.0], rings=4)
        other["burnup"] = [0.0, 2.0]
        with self.assertRaises(vt.TableError):
            vt.depletion_k([depletion([0.0, 0.0]), other])

    def test_a_single_ring_has_no_inner_or_outer_label(self):
        table = vt.depletion_rings([depletion([0.0, 0.0], rings=1)])
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(rows[0][0], "1")

    def test_a_summary_of_another_format_is_an_error(self):
        with self.assertRaises(vt.TableError):
            vt.depletion_k([{"format": "xlc-compare-p3"}])


class DrawAwsDepletion(unittest.TestCase):
    def summary(self) -> dict:
        k = {
            "local": 1.134038,
            "local_std_pcm": 17.7,
            "aws": 1.1335,
            "aws_std_pcm": 21.8,
            "diff_pcm": -53.8,
            "sigma_pcm": 28.0,
            "ratio": -1.92,
        }
        return {
            "format": "xlc-compare-aws-depletion",
            "burnup": [0.1],
            "k": [k],
            "density": {n: [None] for n in vt.DEPLETION_NUCLIDES}
            | {"Pu-239": [-0.0007]},
        }

    def test_k_row_has_the_difference_and_the_ratio_to_the_combined_sigma(self):
        row = [
            c.strip()
            for c in vt.aws_depletion_k([self.summary()])
            .split("\n")[2]
            .strip("|")
            .split("|")
        ]
        self.assertEqual(row[0], "0.1")
        self.assertEqual(row[3:], ["\u221253.8", "28.0", "\u22121.92"])

    def test_density_row_is_in_percent_and_dashes_a_zero_base(self):
        row = [
            c.strip()
            for c in vt.aws_depletion_densities([self.summary()])
            .split("\n")[2]
            .strip("|")
            .split("|")
        ]
        self.assertEqual(row, ["0.1", "–", "–", "\u22120.070", "–", "–"])

    def test_the_git_tracked_json_is_not_checked_against_a_dvc_lock(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = "examples/depletion/aws-check/compare-aws-1C-rings1.json"
            (root / source).parent.mkdir(parents=True)
            (root / source).write_text(json.dumps(self.summary()))
            self.assertEqual(vt.load(root, source, True)["burnup"], [0.1])


class DrawAssemblyPins(unittest.TestCase):
    def summary(self) -> dict:
        """compare_depletion_assembly.py の要約（2 点）"""
        return {
            "format": "xlc-compare-depletion-assembly",
            "problem": "2A",
            "burnup": [0.0, 0.1],
            "pin_std_rms": [0.00062, 0.00065],
            "pins": [
                {
                    "orbits": 39,
                    "rms": 0.0012,
                    "max": 0.0031,
                    "at": [8, 2],
                    "max_std": 0.0007,
                },
                {
                    "orbits": 39,
                    "rms": 0.0015,
                    "max": -0.0042,
                    "at": [0, 0],
                    "max_std": 0.0009,
                },
            ],
        }

    def rows(self) -> list[list[str]]:
        table = vt.depletion_assembly_pins([self.summary()])
        return [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]

    def test_rows_have_the_rms_the_largest_difference_and_the_openmc_error_in_percent(
        self,
    ):
        # 位置は上と左から 1 で数える
        self.assertEqual(
            self.rows(),
            [
                ["0", "0.120", "+0.310 ± 0.070", "(9, 3)", "0.062"],
                ["0.1", "0.150", "\u22120.420 ± 0.090", "(1, 1)", "0.065"],
            ],
        )

    def test_each_problem_has_its_own_table(self):
        sources, draw = vt.TABLES["depletion-2a-pins"]
        self.assertEqual(sources, ["examples/depletion/p2/compare-2A.json"])
        self.assertIs(draw, vt.depletion_assembly_pins)
        self.assertIn("depletion-2a-rods-pins", vt.TABLES)

    def test_a_summary_of_another_format_is_an_error(self):
        with self.assertRaises(vt.TableError):
            vt.depletion_assembly_pins([{"format": "xlc-compare-depletion"}])

    def test_points_without_pins_are_a_table_error(self):
        s = self.summary()
        s["pins"].pop()
        with self.assertRaises(vt.TableError):
            vt.depletion_assembly_pins([s])

    def test_the_problems_are_the_jobs_of_the_assembly_stages(self):
        dvc = (vt.ROOT / "examples/depletion/p2/dvc.yaml").read_text()
        (jobs,) = re.findall(r"^  - jobs: \[(.*)\]$", dvc, re.MULTILINE)
        self.assertEqual(
            vt.ASSEMBLY_DEPLETION_JOBS, [j.strip() for j in jobs.split(",")]
        )


class DrawAssemblyDepletion(unittest.TestCase):
    def summary(self) -> dict:
        """compare_depletion_assembly.py の要約（2 点、燃料と B-10 の吸収体 2 種）"""
        return {
            "format": "xlc-compare-depletion-assembly",
            "problem": "2N",
            "burnup": [0.0, 0.1],
            "k": [
                {
                    "xlc": 1.1,
                    "openmc": 1.0991,
                    "openmc_std_pcm": 11.6,
                    "diff_pcm": 39.3,
                    "change_pcm": 0.0,
                },
                {
                    "xlc": 1.05,
                    "openmc": 1.0495,
                    "openmc_std_pcm": 12.4,
                    "diff_pcm": -59.6,
                    "change_pcm": -98.9,
                },
            ],
            "kinds": {
                "fuel": {
                    "U-235": [0.0, 0.0001],
                    "Pu-239": [None, 0.0012],
                },
                "ifba": {"B-10": [0.0, -0.0213]},
                "waba": {"B-10": [0.0, -0.0139]},
            },
            "seconds": {
                "xlc": {
                    "shielding": 100.04,
                    "moc": 10.0,
                    "rates": 1.0,
                    "matrix": 0.5,
                    "cram": 0.25,
                    "solves": 3,
                },
                "openmc": {},
            },
        }

    def cells(self, table: str) -> list[list[str]]:
        return [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")
        ]

    def test_k_rows_have_the_difference_with_its_error_and_the_change(self):
        rows = self.cells(vt.depletion_assembly_k([self.summary()]))
        self.assertEqual(rows[0], ["燃焼度 (MWd/kg)", "差 (pcm)", "変化 (pcm)"])
        self.assertEqual(rows[2], ["0", "+39 \u00b1 12", "0"])
        self.assertEqual(rows[3], ["0.1", "\u221260 \u00b1 12", "\u221299"])

    def test_fuel_densities_are_in_percent_with_a_dash_for_zero_in_openmc(self):
        rows = self.cells(vt.depletion_assembly_fuel_densities([self.summary()]))
        self.assertEqual(rows[0], ["燃焼度 (MWd/kg)", "U-235", "Pu-239"])
        self.assertEqual(rows[2], ["0", "0.00", "\u2013"])
        self.assertEqual(rows[3], ["0.1", "+0.01", "+0.12"])

    def test_absorber_densities_have_a_column_per_kind_and_nuclide(self):
        rows = self.cells(vt.depletion_assembly_absorber_densities([self.summary()]))
        self.assertEqual(rows[0], ["燃焼度 (MWd/kg)", "ifba B-10", "waba B-10"])
        self.assertEqual(rows[3], ["0.1", "\u22122.13", "\u22121.39"])

    def test_times_have_a_row_per_problem_with_the_sum_as_the_total(self):
        rows = self.cells(vt.depletion_assembly_times([self.summary()]))
        self.assertEqual(
            rows[2], ["2N", "3", "111.8", "100.0", "10.0", "1.00", "0.50", "0.25"]
        )

    def test_each_problem_has_its_k_and_density_tables(self):
        sources, draw = vt.TABLES["depletion-2f-k"]
        self.assertEqual(sources, ["examples/depletion/p2/compare-2F.json"])
        self.assertIs(draw, vt.depletion_assembly_k)
        self.assertIs(
            vt.TABLES["depletion-2n-densities"][1],
            vt.depletion_assembly_fuel_densities,
        )
        self.assertIs(
            vt.TABLES["depletion-2n-absorbers"][1],
            vt.depletion_assembly_absorber_densities,
        )

    def test_a_summary_of_another_format_is_an_error(self):
        for draw in (
            vt.depletion_assembly_k,
            vt.depletion_assembly_fuel_densities,
            vt.depletion_assembly_absorber_densities,
            vt.depletion_assembly_times,
        ):
            with self.assertRaises(vt.TableError):
                draw([{"format": "xlc-compare-depletion"}])


class DrawCoreDepletion(unittest.TestCase):
    def summary(self) -> dict:
        """compare_core_depletion.py の要約（テーブルの点 1 つと点の間 1 つ）"""
        return {
            "format": "xlc-compare-core-depletion",
            "specific_power": 38.41,
            "table_burnups": [0.0, 10.0],
            "points": [
                {
                    "efpd": 0.0,
                    "burnup": 0.0,
                    "on_table": True,
                    "history_k": 1.18283604,
                    "k": 1.18283548,
                    "rho_pcm": -0.0400123,
                    "branch_rho_pcm": -0.00061,
                },
                {
                    "efpd": 130.17443374,
                    "burnup": 5.0,
                    "on_table": False,
                    "history_k": 1.10057334,
                    "k": 1.1021,
                    "rho_pcm": 125.71,
                    "branch_rho_pcm": None,
                },
            ],
        }

    def test_rows_mark_the_table_points_and_have_the_reactivity_difference(self):
        table = vt.core_depletion_k([self.summary()])
        rows = [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]
        self.assertEqual(
            rows,
            [
                ["0", "0.00", "○", "1.18283604", "1.18283548", "−0.040", "−0.0006"],
                ["5", "130.17", "", "1.10057334", "1.10210000", "+125.710", ""],
            ],
        )

    def test_a_summary_of_another_format_is_an_error(self):
        with self.assertRaises(vt.TableError):
            vt.core_depletion_k([{"format": "xlc-compare-depletion"}])

    def test_the_2a_comparison_is_a_git_tracked_json(self):
        (sources, draw) = vt.TABLES["depletion-2a-core-k"]
        self.assertIs(draw, vt.core_depletion_k)
        self.assertTrue(set(sources) <= vt.GIT_TRACKED)


class DrawBurnupInterpolation(unittest.TestCase):
    @staticmethod
    def summary():
        """compare_burnup_interpolation.py の要約（テーブル 1 つ、選び方 2 つ、方式 2 つ）"""

        def result(rho, overshoot):
            return {
                "rho_pcm": rho,
                "max_abs_pcm": max(abs(r) for r in rho),
                "rms_pcm": (sum(r * r for r in rho) / len(rho)) ** 0.5,
                "overshoot": overshoot,
            }

        return {
            "format": "xlc-compare-burnup-interpolation",
            "methods": {"linear": "線形", "hermite": "Hermite"},
            "tables": [
                {
                    "name": "2A",
                    "burnup": [0.0, 1.0, 2.0, 3.0, 4.0],
                    "table_k_max_abs_pcm": 1e-5,
                    "selections": [
                        {
                            "name": "1 つおき",
                            "burnup": [0.0, 2.0, 4.0],
                            "removed": [1.0, 3.0],
                            "results": {
                                "linear": result([-30.0, 40.0], 0.0),
                                "hermite": result([-3.0, 4.0], 3.2e-4),
                            },
                        },
                        {
                            "name": "全点の 3 点",
                            "burnup": [0.0, 1.0, 4.0],
                            "removed": [2.0, 3.0],
                            "results": {
                                "linear": result([60.0, 50.0], 0.0),
                                "hermite": result([-6.0, -5.0], 0.0),
                            },
                        },
                    ],
                }
            ],
        }

    @staticmethod
    def rows(table):
        return [
            [c.strip() for c in line.strip("|").split("|")]
            for line in table.split("\n")[2:]
        ]

    def test_method_rows_have_the_largest_and_rms_errors_of_each_selection(self):
        table = vt.burnup_interpolation_methods([self.summary()])
        header = [c.strip() for c in table.split("\n")[0].strip("|").split("|")]
        self.assertEqual(
            header,
            [
                "方式",
                "2A 1 つおき（3 点）",
                "2A 全点の 3 点（3 点）",
                "オーバーシュートの最大",
            ],
        )
        self.assertEqual(
            self.rows(table),
            [
                ["線形", "40.0 / 35.4", "60.0 / 55.2", "0"],
                ["Hermite", "4.0 / 3.5", "6.0 / 5.5", "3.2e-04"],
            ],
        )

    def test_point_rows_are_the_cubic_hermite_of_each_selection(self):
        table = vt.burnup_interpolation_points([self.summary()])
        self.assertEqual(
            self.rows(table),
            [
                ["2A", "1 つおき", "3", "2", "+4.0（3）", "3.5", "−3.0〜+4.0"],
                ["2A", "全点の 3 点", "3", "2", "−6.0（2）", "5.5", "−6.0〜−5.0"],
            ],
        )

    def test_a_selection_of_every_point_has_no_worst_point(self):
        s = self.summary()
        sel = s["tables"][0]["selections"][0]
        sel["removed"] = []
        for r in sel["results"].values():
            r.update(rho_pcm=[], max_abs_pcm=0.0, rms_pcm=0.0)
        rows = self.rows(vt.burnup_interpolation_points([s]))
        self.assertEqual(rows[0], ["2A", "1 つおき", "3", "0", "", "0.0", ""])

    def test_a_summary_without_selections_draws_the_methods_only(self):
        s = self.summary()
        s["tables"][0]["selections"] = []
        rows = self.rows(vt.burnup_interpolation_methods([s]))
        self.assertEqual(rows, [["線形", "0"], ["Hermite", "0"]])

    def test_a_summary_of_another_format_is_an_error(self):
        for draw in (vt.burnup_interpolation_methods, vt.burnup_interpolation_points):
            with self.assertRaises(vt.TableError):
                draw([{"format": "xlc-compare-core-depletion"}])

    def test_the_comparison_is_a_git_tracked_json(self):
        for name in (
            "depletion-burnup-interpolation-methods",
            "depletion-burnup-interpolation-points",
        ):
            (sources, _) = vt.TABLES[name]
            self.assertTrue(set(sources) <= vt.GIT_TRACKED)


if __name__ == "__main__":
    unittest.main()
