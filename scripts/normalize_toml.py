"""DVC のステージの deps に置く、入力の TOML の値の正規形を書く（ADR-0254）。

使い方: python3 scripts/normalize_toml.py <元の .toml> <正規形の出力>

標準ライブラリの tomllib で読んだ値を、キーの順を保って JSON で書き出す。コメント、空白、改行、引用符と数の書き方は
値に残らないので、正規形に表れない。整数と浮動小数、真偽値と数は JSON でも別の書き方になる。構文の誤りは
TOMLDecodeError、日時（TOML の date と time）は TypeError、nan と inf は ValueError で止まる（元のファイルの md5 に戻さない）。

比較のスクリプトが JSON の `inputs` に書く入力の md5 は、TOML なら正規形の md5（`input_md5`）にする。ステージの
deps が記録する md5 と同じ値になり、コメントだけの修正で JSON が古くならない。
"""

import hashlib
import json
import sys
from pathlib import Path

import tomllib


def normalize(text: str) -> str:
    # 日時は JSON の型にないので、default を渡さず TypeError で止める（文字列にすると、同じ文字列の値と区別できない）
    values = tomllib.loads(text)
    return json.dumps(values, ensure_ascii=False, indent=1, allow_nan=False) + "\n"


def input_md5(path: Path, data: bytes | None = None) -> str:
    """`path` の入力の md5。`data` があれば、ファイルを読み直さずにそのバイトを使う（スクリプトが解いたバイト）"""
    if data is None:
        data = path.read_bytes()
    if path.suffix == ".toml":
        data = normalize(data.decode("utf-8")).encode("utf-8")
    return hashlib.md5(data).hexdigest()


def main(argv: list[str]) -> None:
    if len(argv) != 2:
        raise SystemExit(
            "usage: python3 scripts/normalize_toml.py <source.toml> <output>"
        )
    src, dst = Path(argv[0]), Path(argv[1])
    text = normalize(src.read_text(encoding="utf-8"))
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1:])
