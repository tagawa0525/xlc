#!/usr/bin/env python3
"""dvc status --json から、yaml ごとの changed deps と changed outs のステージの数を出す（#0182）。

使い方: scripts/count-changed-stages.py [<dvc.yaml>...]

引数は dvc status にそのまま渡す。`dvc status` の文字列を grep で数えると、タブで始まる `changed deps:` の行を拾って
誤るので、JSON の「yaml:ステージ」のキーを yaml ごとに数える。stages は
changed のステージの数で、種類ごとの数は 1 つのステージが複数の種類を持てば、それぞれに数える。知らない種類の変更
（not in cache は changed outs の中身で、種類ではない）は ValueError で止める。
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


def count_changed(status):
    counts = {}
    for key, changes in status.items():
        yaml_path = key.split(":", 1)[0] if ":" in key else "dvc.yaml"
        per_yaml = counts.setdefault(
            yaml_path, {"stages": 0, **dict.fromkeys(KINDS, 0)}
        )
        per_yaml["stages"] += 1
        # 変更は {種類: 中身} の辞書か、中身の無い種類（always changed）の文字列
        kinds = {
            kind
            for change in changes
            for kind in ([change] if isinstance(change, str) else change)
        }
        unknown = kinds - set(KINDS)
        if unknown:
            raise ValueError(f"{key}: unknown change kind {sorted(unknown)}")
        for kind in kinds:
            per_yaml[kind] += 1
    return counts


def main(argv):
    result = subprocess.run(
        ["dvc", "status", "--json", *argv[1:]],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        return result.returncode
    counts = count_changed(json.loads(result.stdout))
    for yaml_path in sorted(counts):
        print(
            f"{yaml_path}: "
            + ", ".join(f"{k} {v}" for k, v in counts[yaml_path].items())
        )
    columns = ("stages", *KINDS)
    print(
        "total: "
        + ", ".join(f"{k} {sum(c[k] for c in counts.values())}" for k in columns)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
