#!/usr/bin/env python3
"""Claude Code の PreToolUse の hook: Bash のコマンドが dvc checkout を直接呼べば止める（ADR-0223）。

stdin の JSON（`tool_name`、`tool_input.command`）を読む。コマンドが dvc の checkout を直接呼べば、理由と
scripts/dvc-checkout.py の使い方を stderr に出して終了コード 2（Claude Code がツールの呼び出しを止める）。そうでなければ 0。
コマンドを読めなければ（閉じていない引用）止めずに 1 で知らせる。

hook はシェルから devShell の外で実行されるので、標準ライブラリだけを使う（コマンドの読み方は同じディレクトリの shell_words）。
"""

import json
import sys

from shell_words import runs_dvc

MESSAGE = """\
生の dvc checkout は止めている（ADR-0223）。dvc.lock に hash のまだ無いステージ（流していない参照解）の出力が戻らず、
生の dvc checkout はいつも終了コード 255 で終わるので、hash があるのにキャッシュに無い出力（本当の失敗）が同じ
「Checkout failed」の一覧に埋もれる（#0184）。

代わりにリポジトリ（worktree）の直下で次を使う（default の devShell。direnv で入っていれば nix develop -c は要らない）:

  nix develop -c scripts/dvc-checkout.py [<dvc.yaml か .dvc の指定>...]

指定はリポジトリの直下からのパスで、引数は dvc checkout にそのまま渡る。未実行の出力は UNRUN: に数えるだけで、hash が
あるのに戻らない出力があれば MISSING: に並べて終了コード 1 で止まる。worktree にすべての出力を戻すなら
scripts/worktree-restore-outputs.sh（中で同じスクリプトを呼ぶ）。
"""


def runs_raw_dvc_checkout(command):
    """コマンドが dvc checkout を直接呼ぶか。コマンドを読めなければ ValueError。"""
    return runs_dvc(command, "checkout")


def main():
    try:
        payload = json.load(sys.stdin)
        command = (
            payload.get("tool_input", {}).get("command", "")
            if payload.get("tool_name") == "Bash"
            else ""
        )
        if not isinstance(command, str):
            raise TypeError(f"tool_input.command が文字列でない: {command!r}")
    except (json.JSONDecodeError, AttributeError, TypeError) as e:
        print(f"block-raw-dvc-checkout: hook の入力を読めない: {e}", file=sys.stderr)
        return 1
    try:
        blocked = runs_raw_dvc_checkout(command)
    except ValueError as e:
        print(
            f"block-raw-dvc-checkout: コマンドを読めない（{e}）ので dvc checkout を確かめていない",
            file=sys.stderr,
        )
        return 1
    if blocked:
        print(MESSAGE, file=sys.stderr, end="")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
