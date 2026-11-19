#!/usr/bin/env python3
"""除外と仮の前提の印が、開いている issue を指すことを確かめる（#0294）。

使い方: scripts/check-exclusions.py [--list]

追跡しているファイルから次の印を集める:

- textlintignore: `.textlintignore` の項目。直前の注釈（空行か、前の項目まで）が指す issue を見る
- rust-ignore: crates の `#[ignore]`。その行が指す issue か ADR を見る
- python-skip: 無条件に飛ばすか失敗を許すデコレータ（`@unittest.skip`、`@pytest.mark.skip`、`@pytest.mark.xfail`、
  `@unittest.expectedFailure`）。条件つき（`skipUnless`、`skipIf`、`skipif`）と、テストの本文の中の skip は集めない
- provisional: 例の入力（examples の TOML と YAML）の注釈の「確かめていない」「暫定」「仮の」「仮に」。その行と、
  折り返した次の注釈の行が指す issue を見る（同じかたまりの離れた行の issue は、別の話のことがある）
- excluded: 例の入力の注釈で、除外の語（外す、除いた など）のすぐ前後に issue の番号があるもの

印は、開いている issue（docs/issues/NNNN-*.md の `status: open`）を 1 つ以上指す。rust-ignore と python-skip は
ADR（決めた除外）でもよい。excluded は issue を引く除外だけを集めるので、止めるのは、引いた issue がどれも閉じている
行（理由が消えたのに外したまま）。issue を引かない除外の語は、一般の語（不連続因子を外した定数）なので集めない。
`.textlintignore` が issue のファイルそのものを、その issue の番号で外す項目は、閉じた後も通す（語の調査の表のように、
ファイルが残るかぎり外す）。

試行の入力（パスに `issue-NNNN` か `trial` で始まるディレクトリを含むもの）は見ない: 閉じた issue の試行が、その
issue を引くのは経緯の記録で、本番の形の入力ではない。

`--list` は印の一覧と種類ごとの数を出して 0 で終わる（本番の形の計算を投げる前に、残っている除外を数えるため）。
開いている issue の数と、そのうち名前に review か leftover を含むもの（レビューで直さずに残した指摘）の数も出す。
"""

import re
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

Marker = namedtuple("Marker", "path line kind text issues adrs")

ISSUE = re.compile(r"#(\d{4})(?!\d)")
ADR = re.compile(r"ADR-\d{4}")
ISSUE_FILE = re.compile(r"^docs/issues/(\d{4})-[^/]*\.md$")
TRIAL_DIRECTORY = re.compile(r"(^|/)(issue-\d{4}|trial[^/]*)/")
RUST_IGNORE = re.compile(r"^\s*#\[ignore\b")
# `skip\b` は skipIf、skipUnless、skipif に当たらない（語の途中）
PYTHON_SKIP = re.compile(
    r"^\s*@(unittest\.(skip|expectedFailure)|pytest\.mark\.(skip|xfail))\b"
)
PROVISIONAL = re.compile(r"確かめていない|暫定|仮の|仮に")
# 除外の語と issue の番号が、文をまたがずに 16 文字以内で並ぶもの。離れた issue は、同じ行の別の話のことがある
# （「H-1 を除く全核種に…履歴（#NNNN）の温度」）
DROPPING = r"(外す|外した|外して|除く|除いた|除いて|除外)"
EXCLUDED = re.compile(
    rf"{DROPPING}[^。#]{{0,16}}#\d{{4}}(?!\d)|#\d{{4}}(?!\d)[^。#]{{0,16}}{DROPPING}"
)
REVIEW_LEFTOVER = re.compile(r"review|leftover")
# 開いている issue のほかに、ADR を理由にできる印（決めた除外。参照の出力が要るテストなど）
ADR_ALLOWED = {"rust-ignore", "python-skip"}


def references(text):
    return tuple(ISSUE.findall(text)), tuple(ADR.findall(text))


def textlintignore_markers(path, text):
    comment = ""
    after_entry = False
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            comment = ""
        elif stripped.startswith("#"):
            # 項目の後の注釈は、次のかたまりの理由（前の注釈の issue を引き継がない）
            comment = stripped if after_entry else comment + stripped
            after_entry = False
        else:
            after_entry = True
            yield Marker(path, number, "textlintignore", stripped, *references(comment))


def line_markers(path, text, pattern, kind):
    for number, line in enumerate(text.splitlines(), start=1):
        if pattern.match(line):
            yield Marker(path, number, kind, line.strip(), *references(line))


def input_markers(path, text):
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.strip().startswith("#"):
            continue
        if PROVISIONAL.search(line):
            following = lines[index + 1] if index + 1 < len(lines) else ""
            wrapped = following if following.strip().startswith("#") else ""
            yield Marker(
                path,
                index + 1,
                "provisional",
                line.strip(),
                *references(line + wrapped),
            )
        if EXCLUDED.search(line):
            yield Marker(path, index + 1, "excluded", line.strip(), *references(line))


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


def excludes_its_own_issue(marker):
    match = ISSUE_FILE.match(marker.text) if marker.kind == "textlintignore" else None
    return match is not None and match.group(1) in marker.issues


def violations(markers, statuses):
    """印ごとに、開いている issue（種類によっては ADR）を指さないものを、修正方法つきの文で返す。"""
    found = []
    for marker in markers:
        where = f"{marker.path}:{marker.line}: [{marker.kind}] {marker.text}"
        unknown = [issue for issue in marker.issues if issue not in statuses]
        if unknown:
            missing = "、".join(f"#{issue}" for issue in unknown)
            found.append(
                f"{where}\n  指している issue（{missing}）が、追跡している docs/issues/ にありません。"
                "番号を確かめるか、issue のファイルを git add してください"
            )
            continue
        if any(statuses[issue] == "open" for issue in marker.issues):
            continue
        if marker.kind in ADR_ALLOWED and marker.adrs:
            continue
        if excludes_its_own_issue(marker):
            continue
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


def count_open_issues(statuses, names):
    """issue の番号から状態とファイル名への dict を受け、開いている issue と、そのうちレビューの残りの数を返す。"""
    opened = [names[number] for number, status in statuses.items() if status == "open"]
    return {
        "open-issues": len(opened),
        "open-review-leftover-issues": sum(
            1 for name in opened if REVIEW_LEFTOVER.search(name)
        ),
    }


def issue_statuses(files):
    """パスから中身への dict のうち issue のファイルから、番号から状態への dict と、番号からファイル名への dict を返す。"""
    statuses = {}
    names = {}
    for path in sorted(files):
        match = ISSUE_FILE.match(path)
        if not match:
            continue
        status = re.search(r"^status:\s*(\w+)", files[path], re.MULTILINE)
        if not status:
            raise ValueError(f"{path}: フロントマターに status がありません")
        statuses[match.group(1)] = status.group(1)
        names[match.group(1)] = path.rsplit("/", 1)[-1]
    return statuses, names


def tracked_files(root):
    """追跡しているファイルのうち、印を持ちうるものと issue を読む（追跡していない issue は、コミットされた木に無い）。"""
    names = subprocess.run(
        ["git", "ls-files", "-z"], cwd=root, check=True, capture_output=True, text=True
    ).stdout.split("\0")
    wanted = (".rs", ".py", ".toml", ".yaml", ".yml")
    files = {}
    for name in names:
        if name == ".textlintignore" or name.endswith(wanted) or ISSUE_FILE.match(name):
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
    files = tracked_files(root)
    markers = find_markers(files)
    statuses, names = issue_statuses(files)
    if arguments == ["--list"]:
        for marker in markers:
            cited = " ".join(
                [f"#{issue}({statuses.get(issue, '?')})" for issue in marker.issues]
                + list(marker.adrs)
            )
            print(f"{marker.path}:{marker.line}\t{marker.kind}\t{cited}")
        counts = {**count_by_kind(markers), **count_open_issues(statuses, names)}
        for kind, count in sorted(counts.items()):
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
