"""検証の記録に書く「どの結果か」を、dvc.lock と .dvc から読んで出す（ADR-0176 の帰結。記録を書く担当が lock を読み直さずに済む）。

使い方: python3 examples/record_provenance.py <dvc.yaml> <ステージ>...
        python3 examples/record_provenance.py <x.dvc>...（最初の引数が .dvc で終わるとき）

ステージごとに、同じ場所の dvc.lock から次を出す。
- deps と outs の相対パスと md5（ディレクトリの出力は `.dir` の md5）
- outs ごとに、その md5 を dvc.lock に入れたコミット（8 桁）と日付。コミット済みの dvc.lock の履歴を新しい方から読み、
  ステージのその出力の md5 が今と同じ間をさかのぼって、いちばん古いコミットを取る。relock で md5 が変わらなければコミットは
  変わらない。`git log -S<md5>` は同じ md5 を別のステージが入れたときにも当たるので使わない。取れないときは `unknown`
  （履歴がない）か `uncommitted`（作業ツリーの lock が HEAD と違う）と出し、黙らない。取れるのは、dvc.lock のコミットで、
  実行したコミット（その前）ではない
- 比較のステージ（outs の JSON の `format` が REPORTS にあるもの）は、そのツールの `--report` の文章を続けて出す。JSON が `inputs`
  を持てば、deps の md5 と一致することを確かめ、違えば `MISMATCH` と出して 1 で終わる
- `.dvc`（`dvc add` で追う生成物。ADR-0222）を渡すと、outs ごとに今の md5 と、コミット済みの `.dvc` の履歴を新しい方から
  読んで、md5 が続く間ごとにそれを入れたいちばん古いコミット（8 桁）と日付を出す（同じ md5 に戻れば別の行）。`.dvc` を
  消したコミットか、その出力を持たない版でさかのぼるのをやめる。作業ツリーの `.dvc` の outs が HEAD と違えば今の md5 の行に
  `uncommitted`、履歴がなければ `unknown` と出し、黙らない。取れるのは `.dvc` のコミットで、生成物を作ったコミット
  ではない（md5 が前の版と同じまま作り直すと、前の版を入れたコミットが出る）。生成物の中身が `.dvc` の md5 と合うかは見ない
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
    "xlc-compare-p6": "examples/vera/compare_p6.py",
}


class ProvenanceError(Exception):
    pass


def parse_stage(text: str, stage: str) -> dict:
    """dvc.lock の 1 ステージの deps と outs（(パス, md5) の並び）"""
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
    """ステージ `stage` の出力 `out` の今の md5 を、dvc.lock に入れたコミットと日付。取れなければ ("unknown" か "uncommitted", "")"""
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
    # 出力の md5 だけでなく、ステージの deps と outs のすべてが HEAD と同じこと（コミット前の deps を出さない）
    if head is None or parse_stage(head, stage) != parse_stage(lock.read_text(), stage):
        return "uncommitted", ""
    since = commits[0]
    for commit, date in commits:
        old = git(lock, "show", f"{commit}:./{lock.name}")
        if old is None or out_md5(old, stage, out) != now:
            break
        since = [commit, date]
    return since[0], since[1]


def parse_dvc(text: str) -> list[tuple[str, str]]:
    """.dvc の outs（(パス, md5) の並び）。`dvc add` はキーの順を決めないので、項目ごとに集めてから読む"""
    lines = text.split("\n")
    try:
        start = lines.index("outs:") + 1
    except ValueError:
        raise ProvenanceError(".dvc has no outs") from None
    items: list[dict[str, str]] = []
    key = None
    for line in lines[start:]:
        if line and not line.startswith(("-", " ")):
            break
        if m := re.match(r"(- |  )(\w+):\s*(.*)$", line):
            if m.group(1) == "- ":
                items.append({})
            if not items:
                raise ProvenanceError(
                    f".dvc outs: a key before the first item: {line!r}"
                )
            key = m.group(2)
            items[-1][key] = m.group(3).strip()
            if key == "path" and items[-1][key].startswith(("'", '"')):
                # 引用符と YAML のエスケープを解かずにパスとして出さない（標準ライブラリだけで読む。ADR-0176）
                raise ProvenanceError(f".dvc out with a quoted path: {line!r}")
        elif key == "path" and line.startswith("   "):
            # dvc は長いパスを `path:` の次の、より深く字下げした行に折る。空白を含む名前は空白で何行にも折るが、
            # YAML の行の畳み方（空白の数、引用符）まで読まず、xlc の生成物の名前に無い形としてエラーにする
            if items[-1]["path"]:
                raise ProvenanceError(
                    f".dvc out with a path folded over several lines: {line!r}"
                )
            items[-1]["path"] = line.strip()
    outs = []
    for item in items:
        path, md5 = item.get("path"), item.get("md5", "")
        if not path or not re.fullmatch(r"[0-9a-f]{32}(\.dir)?", md5):
            raise ProvenanceError(f".dvc out without a path or an md5: {item}")
        outs.append((path, md5))
    if not outs:
        raise ProvenanceError(".dvc has no outs")
    return outs


def dvc_history(dvc: Path) -> dict:
    """.dvc の今の outs、状態（None、"uncommitted"、"unknown"）、outs ごとの (md5, それを入れたコミット, 日付) の新しい順"""
    now = parse_dvc(dvc.read_text())
    history: dict[str, list[tuple[str, str, str]]] = {path: [] for path, _ in now}
    log = git(
        dvc, "log", "--abbrev=8", "--format=%h %ad", "--date=short", "--", dvc.name
    )
    if not log or not log.strip():
        return {"now": now, "state": "unknown", "history": history}
    commits = [line.split() for line in log.strip().splitlines()]
    versions = []
    for commit, date in commits:
        old = git(dvc, "show", f"{commit}:./{dvc.name}")
        versions.append((dict(parse_dvc(old)) if old is not None else {}, commit, date))
    state = None if versions[0][0] == dict(now) else "uncommitted"
    for path, runs in history.items():
        # その出力を持たない版（.dvc を消したコミット、別のパス）より前は、別の生成物の履歴なのでさかのぼらない
        for outs, commit, date in versions:
            md5 = outs.get(path)
            if md5 is None:
                break
            if runs and runs[-1][0] == md5:
                runs[-1] = (md5, commit, date)
            else:
                runs.append((md5, commit, date))
    return {"now": now, "state": state, "history": history}


def show_dvc(dvc: Path) -> int:
    if not dvc.exists():
        raise ProvenanceError(f"{dvc} does not exist")
    result = dvc_history(dvc)
    print(f"## {dvc.name}  ({dvc.parent})")
    print(
        "outs (md5 now; then each md5 of the committed .dvc, newest first, with the commit and date that first put it in):"
    )
    for path, md5 in result["now"]:
        print(f"  {md5[:8]}  {path}  {result['state'] or ''}".rstrip())
        for old, commit, date in result["history"][path]:
            print(f"    {old[:8]}  {commit} {date}")
    return 0


def input_problems(summary: dict, stage: dict) -> list[str]:
    """JSON の `inputs` が、ステージの deps（パスと md5）と一致しない点"""
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
    if not argv or (len(argv) < 2 and not argv[0].endswith(".dvc")):
        print(__doc__, file=sys.stderr)
        return 2
    status = 0
    if argv[0].endswith(".dvc"):
        for name in argv:
            try:
                if not name.endswith(".dvc"):
                    raise ProvenanceError(f"{name} is not a .dvc file")
                status |= show_dvc(Path(name))
            except ProvenanceError as e:
                print(f"FAILED: {e}", file=sys.stderr)
                status = 1
        return status
    for name in argv[1:]:
        try:
            status |= show(Path(argv[0]), name)
        except ProvenanceError as e:
            print(f"FAILED: {e}", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
