#!/usr/bin/env python3
"""2 つの git の rev の間で変わった dvc.lock の、deps 以外が同じことを確かめる（#0182）。

使い方: scripts/compare-dvc-lock.py <base> <head>

`git diff --name-only <base>..<head> -- '*dvc.lock'` の各 lock について、ステージの集合と、各ステージの cmd、outs、
params を比べる。`dvc commit` で lock を付け直したブランチは、deps の hash だけが変わるはずである。違いはファイル、
ステージ名、欄を 1 行ずつ出し、1 つでもあれば終了コード 1 で終わる。ステージの集合の差も違いとして出す（増えて当然の
ときは、人が出力を読む。緩めるフラグは作らない）。base か head の一方に無い lock は、ステージが空の lock と比べる。
"""

import subprocess
import sys

import yaml

COMPARED = ("cmd", "outs", "params")


def parse(text):
    if text is None:
        return {}
    return (yaml.safe_load(text) or {}).get("stages") or {}


def compare_locks(name, base_text, head_text):
    base, head = parse(base_text), parse(head_text)
    diffs = [f"{name}: stage {s}: removed" for s in sorted(base.keys() - head.keys())]
    diffs += [f"{name}: stage {s}: added" for s in sorted(head.keys() - base.keys())]
    for stage in sorted(base.keys() & head.keys()):
        diffs += [
            f"{name}: stage {stage}: {field} differs"
            for field in COMPARED
            if base[stage].get(field) != head[stage].get(field)
        ]
    return diffs


def show(rev, path):
    exists = subprocess.run(
        ["git", "cat-file", "-e", f"{rev}:{path}"], capture_output=True, check=False
    )
    if exists.returncode != 0:
        return None
    return subprocess.run(
        ["git", "show", f"{rev}:{path}"], capture_output=True, text=True, check=True
    ).stdout


def main(argv):
    if len(argv) != 3:
        print("usage: compare-dvc-lock.py <base> <head>", file=sys.stderr)
        return 2
    base, head = argv[1], argv[2]
    listed = subprocess.run(
        ["git", "diff", "--name-only", f"{base}..{head}", "--", "*dvc.lock"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    diffs = []
    for path in listed:
        diffs += compare_locks(path, show(base, path), show(head, path))
    print("\n".join(diffs))
    print(f"{len(listed)} lock files compared, {len(diffs)} differences")
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
