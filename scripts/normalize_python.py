"""DVC のステージの deps に置く、Python のファイルの正規形を書く（ADR-0210）。

使い方: python3 scripts/normalize_python.py <元の .py> <正規形の出力>

標準ライブラリの ast で読み、docstring（モジュール、クラス、関数の本体の先頭の文字列の定数の式）を落として
ast.unparse で書き出す。コメント、空行、改行と引用符の書き方は AST に残らないので、正規形に表れない。構文の誤りは
SyntaxError で止まる（元のファイルの md5 に戻さない）。
"""

import ast
import sys
from pathlib import Path

DOCUMENTED = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def is_docstring(node: ast.stmt) -> bool:
    # f-string（JoinedStr）とバイト列は Python の docstring にならない（__doc__ に入らない）ので残す
    return (
        isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    )


def normalize(source: str, filename: str = "<source>") -> str:
    tree = ast.parse(source, filename=filename)
    for node in ast.walk(tree):
        if isinstance(node, DOCUMENTED) and node.body and is_docstring(node.body[0]):
            node.body = node.body[1:]
            # 本体が docstring だけのクラスと関数は、空の本体を書き出せない
            if not node.body and not isinstance(node, ast.Module):
                node.body = [ast.Pass()]
    return ast.unparse(tree) + "\n"


def main(argv: list[str]) -> None:
    if len(argv) != 2:
        raise SystemExit(
            "usage: python3 scripts/normalize_python.py <source.py> <output>"
        )
    src, dst = Path(argv[0]), Path(argv[1])
    text = normalize(src.read_text(encoding="utf-8"), filename=str(src))
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1:])
