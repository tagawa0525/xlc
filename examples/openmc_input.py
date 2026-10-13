"""OpenMC の参照解のスクリプトが入力の TOML を読む部品。

核種ごとのライブラリの段階（openmc_library.py）が、モデルの部品（openmc_common.py）の変更で走り直さないよう、
入力を読む部分だけを分ける（ADR-0136）。
"""

from pathlib import Path
from typing import Any

import tomllib


def check_keys(
    table: dict[str, Any], where: str, required: set[str], optional: set[str]
) -> None:
    """`table` に `required` のキーがすべてあり、`required` と `optional` のほかのキーがないことを確かめる。

    計算の条件は入力に書き、既定値を持たない（ADR-0122、ADR-0130）ので、足りないキーも綴りを誤ったキーも止める。
    """
    missing = required - table.keys()
    unknown = table.keys() - required - optional
    if missing or unknown:
        raise SystemExit(
            f"{where}: missing keys {sorted(missing)}, unknown keys {sorted(unknown)}"
        )


def read_input(path: Path, required: set[str], optional: set[str]) -> dict[str, Any]:
    """入力の TOML を読み、最上位のキーを `check_keys` で確かめる。"""
    table = tomllib.loads(path.read_text())
    check_keys(table, str(path), required, optional)
    return table


def input_path(input_file: Path, value: str) -> Path:
    """入力に書いたパス（入力のファイルからの相対）を、絶対パスにする。"""
    return (input_file.parent / value).resolve()
