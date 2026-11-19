#!/usr/bin/env python3
"""除外と仮の前提の印が、開いている issue を指すことを確かめる（#0294）。

使い方: scripts/check-exclusions.py [--list]

追跡しているファイルから次の印を集める:

- textlintignore: `.textlintignore` の項目。直前の注釈（空行か、前の項目まで）が指す issue を見る
- rust-ignore: crates の `#[ignore]`（`cfg_attr` の中のものを含む）。その行が指す issue か ADR を見る
- python-skip: 無条件に飛ばすか失敗を許すデコレータ（`skip`、`xfail`、`expectedFailure`。`unittest.`、`pytest.mark.`、
  `mark.` の前置きの有無によらない）。括弧が閉じる行までの issue か ADR を見る。条件つき（`skipUnless`、`skipIf`、
  `skipif`）と、テストの本文の中の skip は集めない
- provisional: 例の入力（examples の TOML と YAML）の注釈（値の後ろの注釈を含む）の「確かめていない」「暫定」「仮の」
  「仮に」「仮置き」「未確認」。同じ行が指す issue を見る（前後の行の issue は、別の話のことがある）
- excluded: 例の入力の注釈で、「外す」「外した」「外して」を含むもの（issue の番号の有無によらない。否定は除く）と、
  ほかの除外の語（除いた、省く など）のすぐ前後に issue の番号があるもの。語のすぐ前後（同じ文の NEAR 文字の中）の番号を
  見る（同じ行の離れた番号と、次の行の番号は数えない）

印は、開いている issue（追跡している docs/issues/NNNN-*.md の `status: open`）を 1 つ以上指す。rust-ignore と
python-skip は ADR（決めた除外）でもよい。前提を並べる一覧の issue（LIST_ISSUE）は、項目ごとの解消の条件を持たないので
理由に数えない。「除く」「省く」は形の説明（Hf は除き、1 点ずつ除く）に多く使うので、issue を引く行だけを集める。
`.textlintignore` が issue のファイルそのものを、その issue の番号で外す項目は、閉じた後も通す
（語の調査の表のように、ファイルが残るかぎり外す）。docs/issues/ に無い番号（ほかのリポジトリの issue、色の符号）は
理由に数えず、ほかに理由が無いときに知らせる。

試行の入力（パスに `issue-NNNN` か `trial` で始まるディレクトリを含む TOML と YAML）の注釈は見ない: 閉じた issue の
試行が、その issue を引くのは経緯の記録で、本番の形の入力ではない。

`--list` は印の一覧と種類ごとの数を出して 0 で終わる（本番の形の計算を投げる前に、残っている除外を数えるため）。
開いている issue の数と、そのうち名前に review か leftover を含むもの（レビューで直さずに残した指摘）の数も出す。
"""

import re
import subprocess
import sys
from collections import namedtuple
from pathlib import Path

Marker = namedtuple("Marker", "path line kind text issues adrs")

ISSUE = re.compile(r"#(\d{4})(?!\w)")
ADR = re.compile(r"ADR-\d{4}")
ISSUE_FILE = re.compile(r"^docs/issues/(\d{4})-[^/]*\.md$")
TRIAL_DIRECTORY = re.compile(r"(^|/)(issue-\d{4}|trial)[^/]*/")
RUST_IGNORE = re.compile(r"^\s*#\[(cfg_attr\(.*,\s*)?ignore\b")
# `skip\b` は skipIf、skipUnless、skipif に当たらない（語の途中）
PYTHON_SKIP = re.compile(
    r"^\s*@((unittest|pytest\.mark|mark)\.)?(skip|xfail|expectedFailure)\b"
)
PROVISIONAL = re.compile(r"確かめていない|暫定|仮の|仮に|仮置き|未確認")
# 番号の有無によらず集める除外の語（外す、外した、外して、外し、 など）。否定（外さない、外していない、外してない、
# 外しません）は除外でないので含めない
REMOVING = re.compile(r"外(?:す|し(?!(?:てい?)?(?:な|ませ)))")
DROPPING = re.compile(rf"{REMOVING.pattern}|除[くいき]|省[くいき]|除外")
# 除外の語の前後で、文をまたがずにこの文字数の中にある issue を、その除外の理由とする（16 文字の間と番号の 5 文字）。
# 離れた issue は、同じ行の別の話のことがある（「H-1 を除く全核種に…履歴（#NNNN）の温度」）
NEAR = 21
REVIEW_LEFTOVER = re.compile(r"review|leftover")
# 前提を並べる一覧の issue。項目を指すだけでは、その項目を解消する条件にならない
LIST_ISSUE = "0293"
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
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not pattern.match(line):
            continue
        # 理由の文字列が次の行に折り返したデコレータは、括弧が閉じる行まで（5 行まで）を読む
        end = index + 1
        while end < min(index + 6, len(lines)) and unbalanced(
            "".join(lines[index:end])
        ):
            end += 1
        yield Marker(
            path, index + 1, kind, line.strip(), *references("".join(lines[index:end]))
        )


def unbalanced(text):
    return text.count("(") > text.count(")")


def comment_of(line):
    """行の注釈の部分（引用符の外の最初の `#` から）。無ければ空の文字列。"""
    quote = None
    escaped = False
    for position, character in enumerate(line):
        if quote:
            # 基本文字列（"…"）の中の \" は閉じない。リテラル文字列（'…'）にエスケープは無い
            if escaped:
                escaped = False
            elif character == "\\" and quote == '"':
                escaped = True
            elif character == quote:
                quote = None
        elif character in "\"'":
            quote = character
        elif character == "#":
            return line[position:]
    return ""


def input_markers(path, text):
    for number, line in enumerate(text.splitlines(), start=1):
        comment = comment_of(line)
        if PROVISIONAL.search(comment):
            yield Marker(
                path, number, "provisional", line.strip(), *references(comment)
            )
        adjacent = issues_near_dropping(comment)
        if adjacent or REMOVING.search(comment):
            yield Marker(path, number, "excluded", line.strip(), adjacent, ())


def issues_near_dropping(comment):
    """除外の語ごとに、前後 NEAR 文字の同じ文の中の issue を集める（「#0293 の 1、#0292」のような並びを含む）。"""
    issues = []
    for word in DROPPING.finditer(comment):
        before = comment[max(0, word.start() - NEAR) : word.start()].split("。")[-1]
        after = comment[word.end() : word.end() + NEAR].split("。")[0]
        issues += ISSUE.findall(before) + ISSUE.findall(after)
    return tuple(dict.fromkeys(issues))


def find_markers(files):
    """パスから中身への dict を受け、印をパスと行の順に返す。"""
    markers = []
    for path in sorted(files):
        text = files[path]
        if path == ".textlintignore":
            markers.extend(textlintignore_markers(path, text))
        elif path.startswith("crates/") and path.endswith(".rs"):
            markers.extend(line_markers(path, text, RUST_IGNORE, "rust-ignore"))
        elif path.endswith(".py"):
            markers.extend(line_markers(path, text, PYTHON_SKIP, "python-skip"))
        elif (
            path.startswith("examples/")
            and path.endswith((".toml", ".yaml", ".yml"))
            and not TRIAL_DIRECTORY.search(path)
        ):
            markers.extend(input_markers(path, text))
    return markers


def excludes_its_own_issue(marker):
    match = ISSUE_FILE.match(marker.text) if marker.kind == "textlintignore" else None
    return match is not None and match.group(1) in marker.issues


def violations(markers, statuses):
    """印ごとに、開いている issue（種類によっては ADR）を指さないものを、修正方法つきの文で返す。"""
    found = []
    for marker in markers:
        known = [
            issue
            for issue in marker.issues
            if issue in statuses and issue != LIST_ISSUE
        ]
        if any(statuses[issue] == "open" for issue in known):
            continue
        if marker.kind in ADR_ALLOWED and marker.adrs:
            continue
        if excludes_its_own_issue(marker):
            continue
        where = f"{marker.path}:{marker.line}: [{marker.kind}] {marker.text}"
        place = (
            "除外の語のすぐ隣（同じ文の 16 文字以内）"
            if marker.kind == "excluded"
            else "同じ行"
        )
        unknown = "、".join(
            f"#{issue}" for issue in marker.issues if issue not in statuses
        )
        if known:
            closed = "、".join(f"#{issue}" for issue in known)
            found.append(
                f"{where}\n  指している issue（{closed}）は閉じています。理由が消えたなら外したものを戻し、"
                f"別の理由で残すなら、その理由の開いている issue を{place}に書いてください"
            )
        elif unknown:
            found.append(
                f"{where}\n  指している番号（{unknown}）の issue が、追跡している docs/issues/ にありません。"
                "番号を確かめるか、issue のファイルを git add してください"
            )
        elif LIST_ISSUE in marker.issues:
            found.append(
                f"{where}\n  #{LIST_ISSUE} は前提を並べる一覧で、項目を解消する条件を持たないので理由に数えません。"
                f"その項目を解消する開いている issue を{place}に書いてください"
            )
        else:
            found.append(
                f"{where}\n  開いている issue の番号（#NNNN）が{place}にありません。原因を直すのが既定です。"
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
        number = match.group(1)
        if number in names:
            raise ValueError(
                f"{path}: 番号 {number} の issue が {names[number]} と重なっています"
            )
        front_matter = re.match(r"---\n(.*?)\n---\n", files[path], re.DOTALL)
        status = front_matter and re.search(
            r"^status:[ \t]*(\w+)", front_matter.group(1), re.MULTILINE
        )
        if not status:
            raise ValueError(f"{path}: フロントマターに status の値がありません")
        statuses[number] = status.group(1)
        names[number] = path.rsplit("/", 1)[-1]
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
    try:
        statuses, names = issue_statuses(files)
    except ValueError as error:
        print(f"check-exclusions.py: {error}", file=sys.stderr)
        return 1
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
