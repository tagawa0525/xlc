"""normalize_toml.py の検査（DVC の deps に置く TOML の値の正規形）。

実行: python3 -m unittest discover -s scripts -p 'test_normalize_toml.py'
"""

import hashlib
import json
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
import normalize_toml as norm


def form(text: str) -> str:
    return norm.normalize(textwrap.dedent(text))


BASE = """
title = "pin"
order = ["U-235", "U-238", "O-16"]

[geometry]
pitch = 1.26
rings = 3
reflective = true

[[material]]
name = "fuel"
density = 10.257

[[material]]
name = "water"
density = 0.743
"""


class InvariantTest(unittest.TestCase):
    def test_comment_changes_keep_the_form(self):
        changed = (
            "# 先頭のコメント\n"
            + BASE.replace("rings = 3", "rings = 3  # 行末のコメント")
            + "# 末尾のコメント\n"
        )
        self.assertEqual(form(changed), form(BASE))

    def test_whitespace_and_line_breaks_keep_the_form(self):
        changed = BASE.replace("pitch = 1.26", "pitch   =   1.26").replace(
            'order = ["U-235", "U-238", "O-16"]',
            'order = [\n  "U-235",\n  "U-238",\n  "O-16",\n]\n\n',
        )
        self.assertEqual(form(changed), form(BASE))

    def test_quote_styles_keep_the_form(self):
        self.assertEqual(
            form('a = "x\\\\y"\nb = "\\u00e9"\n'), form("a = 'x\\y'\nb = \"é\"\n")
        )

    def test_number_spellings_of_the_same_value_keep_the_form(self):
        self.assertEqual(
            form("a = 1e10\nb = 0.10\nc = 1_000\nd = 0x10\n"),
            form("a = 1.0e10\nb = 0.1\nc = 1000\nd = 16\n"),
        )

    def test_table_spellings_of_the_same_value_keep_the_form(self):
        self.assertEqual(
            form("[a]\nb = 1\nc.d = 2\n"), form("a = { b = 1, c = { d = 2 } }\n")
        )


class ValueChangeTest(unittest.TestCase):
    def test_value_change_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE.replace("1.26", "1.27")))

    def test_integer_and_float_differ(self):
        self.assertNotEqual(form("a = 1\n"), form("a = 1.0\n"))

    def test_bool_and_integer_differ(self):
        self.assertNotEqual(form("a = true\n"), form("a = 1\n"))

    def test_signed_zero_differs(self):
        self.assertNotEqual(form("a = 0.0\n"), form("a = -0.0\n"))

    def test_key_order_changes_the_form(self):
        # 核種の並びの順で出力を書くプログラムがあるので、表のキーの順も値に数える（ADR-0212 の (b) の条件 3）
        self.assertNotEqual(form("a = 1\nb = 2\n"), form("b = 2\na = 1\n"))

    def test_array_order_changes_the_form(self):
        self.assertNotEqual(
            form(BASE), form(BASE.replace('"U-235", "U-238"', '"U-238", "U-235"'))
        )

    def test_string_and_number_differ(self):
        self.assertNotEqual(form('a = "1"\n'), form("a = 1\n"))

    def test_added_key_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE + "extra = 0\n"))

    def test_form_reads_back_to_the_parsed_values_in_order(self):
        values = tomllib.loads(BASE)
        back = json.loads(norm.normalize(BASE))
        self.assertEqual(back, values)
        self.assertEqual(list(back), list(values))
        self.assertEqual(list(back["geometry"]), list(values["geometry"]))


class FailureTest(unittest.TestCase):
    def test_syntax_error_raises(self):
        with self.assertRaises(tomllib.TOMLDecodeError):
            form("a = \n")

    def test_datetime_raises(self):
        with self.assertRaises(TypeError):
            form("a = 2026-11-19\n")

    def test_nan_and_infinity_raise(self):
        for text in ["a = nan\n", "a = inf\n", "a = -inf\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                form(text)


class MainTest(unittest.TestCase):
    def test_writes_the_form_creating_parent_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.toml")
            src.write_text("# コメント\na = 1\n", encoding="utf-8")
            dst = Path(tmp, "out", "deep", "a.toml.json")
            norm.main([str(src), str(dst)])
            self.assertEqual(dst.read_text(encoding="utf-8"), norm.normalize("a = 1\n"))


class InputMd5Test(unittest.TestCase):
    def test_toml_input_is_the_md5_of_its_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.toml")
            src.write_text("a = 1  # コメント\n", encoding="utf-8")
            self.assertEqual(
                norm.input_md5(src),
                hashlib.md5(norm.normalize("a = 1\n").encode()).hexdigest(),
            )

    def test_other_input_is_the_md5_of_its_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.msgpack")
            src.write_bytes(b"\x81\xa1a\x01")
            self.assertEqual(
                norm.input_md5(src), hashlib.md5(b"\x81\xa1a\x01").hexdigest()
            )


if __name__ == "__main__":
    unittest.main()
