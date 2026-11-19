#!/usr/bin/env python3
"""Claude Code の PreToolUse の hook: slurm-run.sh に渡すジョブのコマンドが dvc repro を直接呼べば止める（ADR-0257）。

stdin の JSON（`tool_name`、`tool_input.command`）を読む。コマンドの中の `slurm-run.sh` の呼び出しの最後の引数（ジョブの
コマンド）が dvc の repro を直接呼べば、理由と scripts/dvc-repro.py の使い方を stderr に出して終了コード 2（Claude Code が
ツールの呼び出しを止める）。そうでなければ 0。`slurm-run.sh` の外の dvc repro（手元の `-s` の再実行など）は止めない。
コマンドを読めなければ（閉じていない引用）止めずに 1 で知らせる。

hook はシェルから devShell の外で実行されるので、標準ライブラリだけを使う（コマンドの読み方は同じディレクトリの shell_words）。
"""

import json
import os
import sys

from shell_words import runs_dvc, simple_commands, unwrap

MESSAGE = """\
Slurm に投げる生の dvc repro は止めている（ADR-0257）。走るステージを数えずに投げると、changed になった参照解（数時間の
OpenMC）や目標の上流のライブラリの作り直しまでが走る（2026-10-22 の Slurm 454、2026-11-19 の Slurm 960）。

代わりに、ジョブのコマンドを次の形にする（default の devShell。数えるのはジョブの中なので、投げる側を待たせない）:

  scripts/dvc-repro.py --max-stages <見込む本数> [--allow-reference <yaml>:<ステージか foreach の組>]... <dvc repro の目標>

走りうるステージ（changed とその下流）が --max-stages を超えるか、--allow-reference に挙げていない参照解のステージ
（cmd を references の devShell で包んだ NJOY と OpenMC）が入れば、何も実行せずに終了コード 1 で止まる。投げる前に数を
見るなら、手元で同じ引数に --count-only を足す。上流を飛ばす数秒の dvc repro -s は Slurm に投げず手元で実行する
（/dvc-relock、ADR-0184）。核種ごとの xs のステージは examples/stage_runner.py で流す。
"""


def submits_raw_dvc_repro(command):
    """コマンドが、dvc repro を直接呼ぶジョブを slurm-run.sh に渡すか。コマンドを読めなければ ValueError。"""
    for argv in simple_commands(command):
        unwrapped = unwrap(argv)
        if isinstance(unwrapped, str):  # bash -c '<文字列>'
            if submits_raw_dvc_repro(unwrapped):
                return True
            continue
        # slurm-run.sh の位置引数は <name> <command> の 2 つで、ジョブのコマンドは最後の引数（slurm-run.sh の Usage）
        if (
            len(unwrapped) >= 3
            and os.path.basename(unwrapped[0]) == "slurm-run.sh"
            and runs_dvc(unwrapped[-1], "repro")
        ):
            return True
    return False


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
        print(
            f"block-raw-dvc-repro-in-slurm: hook の入力を読めない: {e}", file=sys.stderr
        )
        return 1
    try:
        blocked = submits_raw_dvc_repro(command)
    except ValueError as e:
        print(
            f"block-raw-dvc-repro-in-slurm: コマンドを読めない（{e}）ので Slurm に渡す dvc repro を確かめていない",
            file=sys.stderr,
        )
        return 1
    if blocked:
        print(MESSAGE, file=sys.stderr, end="")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
