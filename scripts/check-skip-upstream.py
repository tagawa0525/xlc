#!/usr/bin/env python3
"""`dvc repro -s <yaml>:<ステージ>` の 2 つ目の例外（ADR-0184 の決定の `-s` の項）のうち、`dvc status` で分かる条件を判定する（#0220）。

使い方: scripts/check-skip-upstream.py <yaml>:<ステージ>

`dvc status --with-deps --json <目標>` の目標の上流のステージが、どれも `changed deps` だけを持つ（`always changed` の項目は数えない）
ときに 0、ほかの種類（`changed outs`、`changed command`、`changed checksum`）を持つものがあれば、その一覧を出して 1 で終わる。
目標のステージ自身の変更は再実行の理由なので数えない。知らない種類の変更は ValueError で止める。

`dvc status` に出ない条件は判定しない: `references-shell.txt` は `scripts/check-references-shell.sh` で確かめ、`python-ast` が
changed なら先に実行する（ADR-0184）。
"""

import json
import subprocess
import sys

KINDS = (
    "changed deps",
    "changed outs",
    "changed command",
    "changed checksum",
    "always changed",
)
# `always changed` は always_changed のステージ（references-shell）がいつも出す項目で、deps と outs の差を表さない
ALLOWED_UPSTREAM = {"changed deps", "always changed"}


def status_key(target):
    # 直下の dvc.yaml のステージは、dvc status --json のキーが yaml を持たない（`build`）
    return target.removeprefix("dvc.yaml:")


def blockers(status, target):
    """目標を除く上流のうち、許さない種類の変更を持つステージと、その種類（並べたもの）の一覧。"""
    found = []
    for key, changes in sorted(status.items()):
        # 変更は {種類: 中身} の辞書か、中身の無い種類（always changed）の文字列
        kinds = {
            kind
            for change in changes
            for kind in ([change] if isinstance(change, str) else change)
        }
        unknown = kinds - set(KINDS)
        if unknown:
            raise ValueError(f"{key}: unknown change kind {sorted(unknown)}")
        if key == status_key(target):
            continue
        blocking = sorted(kinds - ALLOWED_UPSTREAM)
        if blocking:
            found.append((key, blocking))
    return found


def main(argv):
    if len(argv) != 2 or ":" not in argv[1]:
        print("usage: scripts/check-skip-upstream.py <yaml>:<stage>", file=sys.stderr)
        return 2
    target = argv[1]
    result = subprocess.run(
        ["dvc", "status", "--with-deps", "--json", target],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        return result.returncode
    found = blockers(json.loads(result.stdout), target)
    for key, kinds in found:
        print(f"BLOCKED: {key}: {', '.join(kinds)}")
    if found:
        print(
            "ERROR: 上流に changed deps でない変更がある。-s を使わず dvc repro <目標> で上流から実行する（ADR-0184）",
            file=sys.stderr,
        )
        return 1
    print(
        "OK: 上流は changed deps だけ。references-shell.txt は scripts/check-references-shell.sh で確かめ、"
        "python-ast が changed なら先に実行する（ADR-0184）"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
