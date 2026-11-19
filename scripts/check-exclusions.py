#!/usr/bin/env python3
"""除外と仮の前提の印が、開いている issue を指すことを確かめる（#0294）。

使い方: scripts/check-exclusions.py [--list]

追跡しているファイルから次の印を集める:

- textlintignore: `.textlintignore` の項目。直前の注釈（空行まで）が指す issue を見る
- rust-ignore: crates の `#[ignore]`。その行が指す issue か ADR を見る
- python-skip: 無条件の `@unittest.skip(` と `@pytest.mark.skip(`。条件つき（`skipUnless`、`skipIf`）は集めない
- provisional: 例の入力（examples の TOML と YAML）の注釈の「確かめていない」「暫定」「仮の」「仮に」。同じ注釈の
  かたまり（注釈の行の連なり）が指す issue を見る
- excluded: 例の入力の注釈で、同じ行に「外す」「外した」「外して」と issue の番号があるもの

印は、開いている issue（docs/issues/NNNN-*.md の `status: open`）を 1 つ以上指す。rust-ignore と python-skip は
ADR（決めた除外）でもよい。excluded は issue を引く行だけを集めるので、止めるのは、引いた issue がどれも閉じている
行（理由が消えたのに外したまま）。issue を引かない「外す」は、一般の語（不連続因子を外した定数）が多いので集めない。

試行の入力（パスに `issue-NNNN` か `trial` で始まるディレクトリを含むもの）は見ない: 閉じた issue の試行が、その
issue を引くのは経緯の記録で、本番の形の入力ではない。

`--list` は印の一覧と種類ごとの数を出して 0 で終わる（本番の形の計算を投げる前に、残っている除外を数えるため）。
"""

import re
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

Marker = namedtuple("Marker", "path line kind text issues adrs")

ISSUE = re.compile(r"#(\d{4})(?!\d)")
ADR = re.compile(r"ADR-\d{4}")
TRIAL_DIRECTORY = re.compile(r"(^|/)(issue-\d{4}|trial[^/]*)/")
RUST_IGNORE = re.compile(r"^\s*#\[ignore\b")
PYTHON_SKIP = re.compile(r"^\s*@(unittest\.skip|pytest\.mark\.skip)\(")
PROVISIONAL = re.compile(r"確かめていない|暫定|仮の|仮に")
EXCLUDED = re.compile(r"外す|外した|外して")
# 開いている issue のほかに、ADR を理由にできる印（決めた除外。参照の出力が要るテストなど）
ADR_ALLOWED = {"rust-ignore", "python-skip"}


def references(text):
    return tuple(ISSUE.findall(text)), tuple(ADR.findall(text))


def textlintignore_markers(path, text):
    comment = ""
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            comment = ""
        elif stripped.startswith("#"):
            comment += stripped
        else:
            yield Marker(path, number, "textlintignore", stripped, *references(comment))


def line_markers(path, text, pattern, kind):
    for number, line in enumerate(text.splitlines(), start=1):
        if pattern.match(line):
            yield Marker(path, number, kind, line.strip(), *references(line))


def input_markers(path, text):
    lines = text.splitlines()
    block_start = None
    for index, line in enumerate(lines + [""]):
        if line.strip().startswith("#"):
            if block_start is None:
                block_start = index
            continue
        if block_start is not None:
            block = lines[block_start:index]
            block_references = references("".join(block))
            for offset, comment in enumerate(block):
                number = block_start + offset + 1
                if PROVISIONAL.search(comment):
                    yield Marker(
                        path, number, "provisional", comment.strip(), *block_references
                    )
                issues, adrs = references(comment)
                if EXCLUDED.search(comment) and issues:
                    yield Marker(
                        path, number, "excluded", comment.strip(), issues, adrs
                    )
            block_start = None


def find_markers(files):
    """パスから中身への dict を受け、印をパスと行の順に返す。"""
    markers = []
    for path in sorted(files):
        text = files[path]
        if TRIAL_DIRECTORY.search(path):
            continue
        if path == ".textlintignore":
            markers.extend(textlintignore_markers(path, text))
        elif path.startswith("crates/") and path.endswith(".rs"):
            markers.extend(line_markers(path, text, RUST_IGNORE, "rust-ignore"))
        elif path.endswith(".py"):
            markers.extend(line_markers(path, text, PYTHON_SKIP, "python-skip"))
        elif path.startswith("examples/") and path.endswith((".toml", ".yaml", ".yml")):
            markers.extend(input_markers(path, text))
    return markers


def violations(markers, statuses):
    """印ごとに、開いている issue（種類によっては ADR）を指さないものを、修正方法つきの文で返す。"""
    found = []
    for marker in markers:
        unknown = [issue for issue in marker.issues if issue not in statuses]
        if unknown:
            raise ValueError(
                f"{marker.path}:{marker.line}: issue #{unknown[0]} が docs/issues/ にありません"
            )
        if any(statuses[issue] == "open" for issue in marker.issues):
            continue
        if marker.kind in ADR_ALLOWED and marker.adrs:
            continue
        where = f"{marker.path}:{marker.line}: [{marker.kind}] {marker.text}"
        if marker.issues:
            closed = "、".join(f"#{issue}" for issue in marker.issues)
            found.append(
                f"{where}\n  指している issue（{closed}）は閉じています。理由が消えたなら外したものを戻し、"
                "別の理由で残すなら、その理由の開いている issue を書いてください"
            )
        else:
            found.append(
                f"{where}\n  開いている issue の番号（#NNNN）がありません。原因を直すのが既定です。"
                "残すなら、解消の条件を書いた issue を起こして番号を書いてください"
            )
    return found


def count_by_kind(markers):
    counts = {}
    for marker in markers:
        counts[marker.kind] = counts.get(marker.kind, 0) + 1
    return counts


def issue_statuses(root):
    statuses = {}
    for path in sorted((root / "docs" / "issues").glob("[0-9][0-9][0-9][0-9]-*.md")):
        match = re.search(
            r"^status:\s*(\w+)", path.read_text(encoding="utf-8"), re.MULTILINE
        )
        if not match:
            raise ValueError(f"{path}: フロントマターに status がありません")
        statuses[path.name[:4]] = match.group(1)
    return statuses


def tracked_files(root):
    names = subprocess.run(
        ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True, text=True
    ).stdout.split("\0")
    wanted = (".rs", ".py", ".toml", ".yaml", ".yml")
    files = {}
    for name in names:
        if name == ".textlintignore" or name.endswith(wanted):
            path = root / name
            # index にあって作業ツリーに無いファイル（消す途中）は印を持たない
            if path.is_file():
                files[name] = path.read_text(encoding="utf-8")
    return files


def main(arguments):
    if arguments not in ([], ["--list"]):
        print(__doc__, file=sys.stderr)
        return 2
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )
    markers = find_markers(tracked_files(root))
    statuses = issue_statuses(root)
    if arguments == ["--list"]:
        for marker in markers:
            cited = " ".join(
                [f"#{issue}({statuses.get(issue, '?')})" for issue in marker.issues]
                + list(marker.adrs)
            )
            print(f"{marker.path}:{marker.line}\t{marker.kind}\t{cited}")
        for kind, count in sorted(count_by_kind(markers).items()):
            print(f"COUNT: {kind} {count}")
        return 0
    found = violations(markers, statuses)
    for problem in found:
        print(problem, file=sys.stderr)
    if found:
        print(f"check-exclusions.py: {len(found)} 件（一覧は --list）", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
