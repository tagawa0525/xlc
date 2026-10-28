#!/usr/bin/env python3
"""dvc checkout を実行し、戻らなかった出力を、未実行のステージの出力とそれ以外に分ける（#0184）。

使い方: scripts/dvc-checkout.py [<dvc.yaml の指定>...]（default の devShell で。指定はリポジトリの直下からのパス）

引数は dvc checkout にそのまま渡し、リポジトリの直下で実行する。dvc.lock に hash がまだ無いステージ（流していない参照解）の
出力は戻らず、dvc checkout は終了コード 255 で終わる。dvc がその出力について `No file hash info found` と警告したものは
未実行として数えるだけにし、hash があるのに戻らなかった出力（キャッシュに無い、ブランチが変えた lock の出力）が 1 つでも
あれば終了コード 1 で止める。`dvc checkout --allow-missing` を使わないのは、hash があってキャッシュに無い出力も、警告なしに
飛ばして 0 で終わるため（dvc 3.67.1。#0184）。dvc の出力の形を読めなければ、推測せずに止める。
"""

import subprocess
import sys
from pathlib import Path

NO_HASH_PREFIX = "WARNING: No file hash info found for '"
NO_HASH_SUFFIX = "'. It won't be created."
FAILED_HEADER = "ERROR: Checkout failed for following targets:"
FAILED_END = "Is your cache up to date?"


def classify(output, root):
    """dvc checkout の出力から、戻らなかった出力を (未実行, それ以外) の 2 つの並べたリストにする。

    パスはリポジトリの直下からの形。失敗の一覧が無ければ ValueError。
    """
    lines = output.splitlines()
    no_hash = set()
    for line in lines:
        if line.startswith(NO_HASH_PREFIX) and line.endswith(NO_HASH_SUFFIX):
            path = Path(line[len(NO_HASH_PREFIX) : -len(NO_HASH_SUFFIX)])
            try:
                no_hash.add(path.relative_to(root).as_posix())
            except ValueError as e:
                raise ValueError(
                    f"hash の無い出力がリポジトリの外にある: {path}"
                ) from e
    try:
        start = lines.index(FAILED_HEADER) + 1
        end = lines.index(FAILED_END, start)
    except ValueError as e:
        raise ValueError(
            "dvc checkout の出力に、戻らなかった出力の一覧が無い（上の出力を読む）"
        ) from e
    failed = set(lines[start:end])
    return sorted(failed & no_hash), sorted(failed - no_hash)


def main(argv):
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    )
    result = subprocess.run(
        ["dvc", "checkout", *argv[1:]],
        cwd=root,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    sys.stdout.write(result.stdout)
    if result.returncode == 0:
        return 0
    try:
        unrun, missing = classify(result.stdout, root)
    except ValueError as e:
        print(
            f"ERROR: dvc checkout が終了コード {result.returncode} で終わった: {e}",
            file=sys.stderr,
        )
        return 1
    print(f"UNRUN: {len(unrun)}（dvc.lock に hash の無い、未実行のステージの出力）")
    if missing:
        print(
            "ERROR: hash があるのに戻らなかった出力がある。キャッシュに無いか、ブランチが変えた"
            " dvc.lock の出力。手で dvc repro で埋めずに原因を調べる:",
            file=sys.stderr,
        )
        for path in missing:
            print(f"MISSING: {path}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
