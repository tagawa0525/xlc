"""検証の記録に書く「どの結果か」を、dvc.lock から読んで出す（ADR-0176 の帰結。記録を書く担当が lock を読み直さずに済む）。

使い方: python3 examples/record_provenance.py <dvc.yaml> <段階>...

段階ごとに、同じ場所の dvc.lock から次を出す。
- deps と outs の相対パスと md5（ディレクトリの出力は `.dir` の md5）
- outs ごとに、その md5 を dvc.lock に入れたコミット（8 桁）と日付。コミット済みの dvc.lock の履歴を新しい方から読み、
  段階のその出力の md5 が今と同じ間をさかのぼって、いちばん古いコミットを取る。relock で md5 が変わらなければコミットは
  変わらない。`git log -S<md5>` は同じ md5 を別の段階が入れたときにも当たるので使わない。取れないときは `unknown`
  （履歴がない）か `uncommitted`（作業ツリーの lock が HEAD と違う）と出し、黙らない。取れるのは、dvc.lock のコミットで、
  実行したコミット（その前）ではない
- 比較の段階（outs の JSON の `format` が REPORTS にあるもの）は、そのツールの `--report` の文章を続けて出す。JSON が `inputs`
  を持てば、deps の md5 と一致することを確かめ、違えば `MISMATCH` と出して 1 で終わる
- 標準ライブラリだけ
"""

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# JSON の `format` → `--report` を持つツール（リポジトリの直下から）
REPORTS = {
    "xlc-compare-p3": "examples/vera/compare_p3.py",
    "xlc-compare-p5": "examples/vera/compare_p5.py",
}


class ProvenanceError(Exception):
    pass


def parse_stage(text: str, stage: str) -> dict:
    """dvc.lock の 1 段階の deps と outs（(パス, md5) の並び）"""
    lines = text.split("\n")
    try:
        start = lines.index(f"  {stage}:") + 1
    except ValueError:
        raise ProvenanceError(f"stage {stage} is not in dvc.lock") from None
    result: dict[str, list[tuple[str, str]]] = {"deps": [], "outs": []}
    section = None
    path = None
    for line in lines[start:]:
        if re.match(r"  \S", line):
            break
        if m := re.match(r"    (deps|outs):$", line):
            section = m.group(1)
        elif re.match(r"    [^\s-]", line):
            section = None
        elif section and (m := re.match(r"    - path: (.+)$", line)):
            path = m.group(1)
        elif section and re.match(r"    - path:\s*$", line):
            path = ""  # dvc は長いパスを次の行に折る
        elif section and path == "" and line.strip():
            # 続きの行は `- path:` より深く字下げした最初の非空白の行。キーの行なら、パスが無いまま md5 に進むので止める
            if len(line) - len(line.lstrip()) <= 4 or re.match(
                r"\s*(hash|md5|size|nfiles):", line
            ):
                raise ProvenanceError(
                    f"dvc.lock stage {stage}: `- path:` has no path on the next line: {line!r}"
                )
            path = line.strip()
        elif (
            section
            and path
            and (m := re.match(r"      md5: ([0-9a-f]{32}(\.dir)?)$", line))
        ):
            result[section].append((path, m.group(1)))
            path = None
    return result


def git(lock: Path, *args: str) -> str | None:
    run = subprocess.run(
        ["git", "-C", str(lock.parent), *args],
        capture_output=True,
        text=True,
        check=False,
    )
    return run.stdout if run.returncode == 0 else None


def out_md5(text: str, stage: str, out: str) -> str | None:
    try:
        outs = dict(parse_stage(text, stage)["outs"])
    except ProvenanceError:
        return None
    return outs.get(out)


def locked_since(lock: Path, stage: str, out: str) -> tuple[str, str]:
    """段階 `stage` の出力 `out` の今の md5 を、dvc.lock に入れたコミットと日付。取れなければ ("unknown" か "uncommitted", "")"""
    now = out_md5(lock.read_text(), stage, out)
    if now is None:
        raise ProvenanceError(f"{lock}: stage {stage} has no output {out}")
    history = git(
        lock, "log", "--abbrev=8", "--format=%h %ad", "--date=short", "--", lock.name
    )
    if not history or not history.strip():
        return "unknown", ""
    commits = [line.split() for line in history.strip().splitlines()]
    head = git(lock, "show", f"{commits[0][0]}:./{lock.name}")
    # 出力の md5 だけでなく、段階の deps と outs のすべてが HEAD と同じこと（コミット前の deps を出さない）
    if head is None or parse_stage(head, stage) != parse_stage(lock.read_text(), stage):
        return "uncommitted", ""
    since = commits[0]
    for commit, date in commits:
        old = git(lock, "show", f"{commit}:./{lock.name}")
        if old is None or out_md5(old, stage, out) != now:
            break
        since = [commit, date]
    return since[0], since[1]


def input_problems(summary: dict, stage: dict) -> list[str]:
    """JSON の `inputs` が、段階の deps（パスと md5）と一致しない点"""
    deps = dict(stage["deps"])
    problems = []
    for item in summary.get("inputs", []):
        if item["path"] not in deps:
            problems.append(f"input {item['path']} is not a dep of the stage")
        elif deps[item["path"]] != item["md5"]:
            problems.append(
                f"input {item['path']}: JSON {item['md5']} but dvc.lock {deps[item['path']]}"
            )
    return problems


def show(dvc_yaml: Path, stage_name: str) -> int:
    lock = dvc_yaml.with_name("dvc.lock")
    if not lock.exists():
        raise ProvenanceError(f"{lock} does not exist")
    stage = parse_stage(lock.read_text(), stage_name)
    print(f"## {stage_name}  ({dvc_yaml.parent})")
    print("deps:")
    for path, md5 in stage["deps"]:
        print(f"  {md5[:8]}  {path}")
    print("outs (md5, the dvc.lock commit and date that first locked it):")
    for path, md5 in stage["outs"]:
        commit, date = locked_since(lock, stage_name, path)
        print(f"  {md5[:8]}  {path}  {commit} {date}".rstrip())
    status = 0
    for path, _ in stage["outs"]:
        file = dvc_yaml.parent / path
        if file.suffix != ".json":
            continue
        if not file.exists():
            print(f"MISSING: {file} (run the stage, or git checkout it)")
            status = 1
            continue
        if hashlib.md5(file.read_bytes()).hexdigest() != dict(stage["outs"])[path]:
            print(f"MISMATCH: {file} differs from the md5 in dvc.lock")
            status = 1
            continue
        summary = json.loads(file.read_text())
        tool = REPORTS.get(summary.get("format"))
        if tool is None:
            continue
        if "inputs" not in summary:
            print(f"MISMATCH: {file} has no inputs; rerun the stage")
            status = 1
            continue
        for problem in input_problems(summary, stage):
            print(f"MISMATCH: {problem}")
            status = 1
        run = subprocess.run(
            [sys.executable, str(ROOT / tool), "--report", str(file)],
            capture_output=True,
            text=True,
            check=False,
        )
        if run.returncode:
            print(f"report failed: {run.stderr.strip()}")
            status = 1
        else:
            print(f"report ({file.name}):\n{run.stdout.rstrip()}")
    return status


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    status = 0
    for name in argv[1:]:
        try:
            status |= show(Path(argv[0]), name)
        except ProvenanceError as e:
            print(f"FAILED: {e}", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
