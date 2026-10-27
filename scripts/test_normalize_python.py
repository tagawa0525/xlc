"""normalize_python.py の検査（DVC の deps に置く Python の正規形）。

実行: python3 -m unittest discover -s scripts -p 'test_normalize_python.py'
"""

import ast
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import normalize_python as norm


def form(source: str) -> str:
    return norm.normalize(textwrap.dedent(source))


BASE = '''
"""モジュールの説明。"""

import math


class Model:
    """クラスの説明。"""

    def k(self, x):
        """メソッドの説明。"""

        def inner(y):
            """ネストした関数の説明。"""
            return y * 2

        return inner(x) + math.pi


async def fetch(path):
    """非同期の関数の説明。"""
    return path
'''


class InvariantTest(unittest.TestCase):
    def test_docstring_changes_keep_the_form(self):
        changed = (
            BASE.replace("モジュールの説明。", "モジュールの説明を書き直した。")
            .replace("クラスの説明。", "別の文。")
            .replace(
                "メソッドの説明。", "メソッドの説明。\n\n        段落を足す。\n        "
            )
            .replace("ネストした関数の説明。", "x")
            .replace("非同期の関数の説明。", "")
        )
        self.assertNotEqual(changed, BASE)
        self.assertEqual(form(changed), form(BASE))

    def test_comment_changes_keep_the_form(self):
        changed = BASE.replace(
            "import math\n", "# 先頭のコメント\nimport math  # 行末のコメント\n"
        ).replace("return y * 2", "return y * 2  # 2 倍する")
        self.assertEqual(form(changed), form(BASE))

    def test_blank_lines_and_formatting_keep_the_form(self):
        a = """
        x = f(1, 2)
        y = {'a': 1}
        """
        b = """
        x = f(
            1,
            2,
        )


        y = {"a": (1)}
        """
        self.assertEqual(form(a), form(b))


class CodeChangeTest(unittest.TestCase):
    def test_expression_change_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE.replace("y * 2", "y * 3 - 1")))

    def test_name_change_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE.replace("inner", "outer")))

    def test_constant_change_changes_the_form(self):
        self.assertNotEqual(form("x = 112.5\n"), form("x = 112.47\n"))


class DocstringTest(unittest.TestCase):
    def test_drops_module_class_function_and_nested_docstrings(self):
        text = form(BASE)
        for doc in ["モジュール", "クラス", "メソッド", "ネストした", "非同期"]:
            self.assertNotIn(doc, text)
        self.assertIn("return inner(x) + math.pi", text)

    def test_keeps_string_expressions_that_are_not_docstrings(self):
        text = form(
            """
            def f(x):
                f"{x} は docstring でない"
                label = "代入の右辺"
                "2 番目の文"
                return label
            """
        )
        self.assertIn("docstring でない", text)
        self.assertIn("代入の右辺", text)
        self.assertIn("2 番目の文", text)

    def test_body_with_only_a_docstring_stays_valid(self):
        text = form(
            '''
            class Empty:
                """説明だけ。"""

            def f():
                """説明だけ。"""
            '''
        )
        ast.parse(text)
        self.assertNotIn("説明だけ", text)


class FailureTest(unittest.TestCase):
    def test_syntax_error_raises(self):
        with self.assertRaises(SyntaxError):
            form("def f(:\n")


class MainTest(unittest.TestCase):
    def test_writes_the_form_creating_parent_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.py")
            src.write_text('"""説明。"""\nx = 1  # コメント\n', encoding="utf-8")
            dst = Path(tmp, "out", "deep", "a.py.ast")
            norm.main([str(src), str(dst)])
            self.assertEqual(dst.read_text(encoding="utf-8"), "x = 1\n")


if __name__ == "__main__":
    unittest.main()
