#!/usr/bin/env python3
"""rebase の後に、出力を dvc.lock のとおりに戻し、実行ファイルと deps の正規形を書き直し、dvc.lock を確かめる（ADR-0258 の
「rebase の後の付け直し」）。

使い方（default の devShell で、リポジトリの直下から）:
    scripts/dvc-relock.py [--base <ブランチ>] [--run <yaml>:<ステージか foreach の組>]...
                          [--allow-reference <yaml>:<ステージか組>]... [<dvc.yaml>...]

次を順に行い、止まるところでは後を実行せずに終了コード 1 で止まる:

0. 分岐点（`git merge-base HEAD <base>`）から作業ツリーまでに変わった dvc.lock を scripts/check-dvc-lock.sh で確かめる。rebase の
   pick では pre-commit の hook が走らないので、ブランチのコミットの分も見る。作業ツリーを書き換える前に行う（読めない lock を
   入れたコミットを直す fixup と rebase は、汚れた作業ツリーでは動かない）。読めない lock があれば止まる
1. scripts/dvc-checkout.py。hash があるのに戻らない出力があれば止まる
2. 直下の dvc.yaml の build、python-ast、binary-forms、toml-values を 1 つの `dvc repro -s` で実行する（dvc が DAG の順に並べる）。
   実行ファイルと正規形は止まらずに書き直す（下流の changed が今のコードを示すように、dvc status の前に行う）。dvc repro が
   失敗すれば止まる
3. `dvc status --json`（<dvc.yaml> を渡せば `--with-deps` でその範囲）。失敗すれば止まる。`changed outs` だけのステージか .dvc が
   あれば止まる（出力が dvc.lock と違うのに deps は同じで、再実行の理由にならない。原因を調べる）
4. scripts/check-references-shell.sh。references-shell.txt は devShell が変わっても dvc status に出ないので、ここで比べ、食い違えば
   止まる
5. `--run` に挙げたステージだけを、挙げた順に `scripts/dvc-repro.py -s` で実行する（上流の確認は scripts/check-skip-upstream.py の
   判定）。挙げなければ、計算のステージは実行しない。止まれば、後のものは実行しない。実行したら dvc status を読み直し、3 と同じく
   `changed outs` だけのものがあれば止まる
6. dvc が書き換えた lock も含めて、0 をもう一度行う

終わりに（止まったときも）、実行したステージ（`RAN:`）と、コミットするファイル（`COMMIT:`）を出す。`COMMIT:` は、このコマンドが
書き換えたファイルと、コミットしていない dvc.lock、python-ast/、binary-forms/、toml-values/ のファイル（止まった前の実行が書き換えた
もの）である。始める前から変わっていたほかのファイル（rebase の hook が付け替えた記録など）は、このコマンドが書き換えて
いなければ出さない。通ったときは、残した changed のステージと変更の種類（`CHANGED:`。always changed だけの references-shell は 4 で
確かめたので出さない）と、yaml ごとの数（`COUNT:`。scripts/count-changed-stages.py の形）も出す。コミットは呼んだ側が行う。
"""

import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
FORMS = (
    "build",
    "dvc.yaml:python-ast",
    "dvc.yaml:binary-forms",
    "dvc.yaml:toml-values",
)
# このコマンドの部品が書き換える追跡するファイル（dvc.lock と正規形）
WRITTEN = ("python-ast/", "binary-forms/", "toml-values/")
ALWAYS = "always changed"

# 変更の種類の検査と yaml ごとの数は count-changed-stages.py と同じものを使う（名前に - があるので import 文では読めない）
_SPEC = importlib.util.spec_from_file_location(
    "count_changed_stages", SCRIPTS / "count-changed-stages.py"
)
counting = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(counting)


class Stop(Exception):
    pass


def kinds_of(status: dict) -> dict[str, list[str]]:
    """dvc status --json の key → 変更の種類（並べたもの）。知らない種類は ValueError"""
    counting.count_changed(status)
    return {
        key: sorted(counting.change_kinds(changes)) for key, changes in status.items()
    }


def outs_only(status: dict) -> list[str]:
    """`changed outs` のほかに変更を持たないステージと .dvc（always changed は数えない）"""
    return sorted(
        key
        for key, kinds in kinds_of(status).items()
        if set(kinds) - {ALWAYS} == {"changed outs"}
    )


def left_changed(status: dict) -> list[tuple[str, list[str]]]:
    """残した changed のステージと種類（always changed だけのものを除く）"""
    return sorted(
        (key, kinds)
        for key, kinds in kinds_of(status).items()
        if set(kinds) != {ALWAYS}
    )


def matches(entry: str, target: str) -> bool:
    """`--allow-reference` の 1 つが、`--run` の目標か、その foreach の組か、組の 1 つを指すか（dvc-repro.py の allows と同じ向き）"""
    return (
        entry == target
        or entry.split("@", 1)[0] == target
        or target.split("@", 1)[0] == entry
    )


def git(*args: str, stdin: str | None = None) -> str:
    return subprocess.run(
        ["git", *args], input=stdin, capture_output=True, text=True, check=True
    ).stdout


def worktree_files() -> dict[str, str]:
    """git status に出る追跡するファイルと追跡しないファイル → 中身の hash（消えたものは "deleted"）"""
    fields = git("status", "--porcelain=v1", "-z", "--untracked-files=all").split("\0")
    paths = []
    i = 0
    while i < len(fields):
        entry = fields[i]
        i += 1
        if not entry:
            continue
        paths.append(entry[3:])
        if entry[0] in "RC":  # 名前を変えたものは、次の欄が元の名前
            i += 1
    present = [p for p in paths if Path(p).is_file()]
    hashes = (
        git("hash-object", "--stdin-paths", stdin="\n".join(present)).split()
        if present
        else []
    )
    state = dict.fromkeys(paths, "deleted")
    state.update(zip(present, hashes))
    return state


def run(command: list[str], what: str, quiet: bool = False) -> None:
    """quiet なら、dvc の出力を失敗したときだけ出す（正規形の dvc repro は、飛ばした数百のステージを 1 行ずつ出す）"""
    print(f"RUNNING: {' '.join(command)}", flush=True)
    result = subprocess.run(command, capture_output=quiet, text=True, check=False)
    if result.returncode != 0:
        if quiet:
            sys.stdout.write(result.stdout)
            sys.stderr.write(result.stderr)
        raise Stop(what)


def dvc_status(targets: list[str]) -> dict:
    """dvc status --json を読み、`changed outs` だけのステージがあれば止める"""
    command = [
        "dvc",
        "status",
        "--json",
        *(["--with-deps", *targets] if targets else []),
    ]
    print(f"RUNNING: {' '.join(command)}", flush=True)
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise Stop(
            f"{' '.join(command)} failed with exit code {result.returncode}: a dvc.lock or dvc.yaml that dvc cannot read, "
            "or a broken reference to the cache. Find the cause in the dvc message above"
        )
    try:
        status = json.loads(result.stdout)
        stuck = outs_only(status)
    except (json.JSONDecodeError, ValueError) as e:
        raise Stop(f"cannot read the output of {' '.join(command)}: {e}") from e
    if stuck:
        raise Stop(
            "stages with only changed outs (their outputs differ from dvc.lock while their deps match it, so "
            "running them is not the fix; find out what changed the outputs): "
            + ", ".join(stuck)
        )
    return status


def check_locks(base: str) -> None:
    fork = git("merge-base", "HEAD", base).strip()
    locks = [
        p
        for p in git(
            "diff", "--name-only", fork, "--", "dvc.lock", "*/dvc.lock"
        ).split()
        if Path(p).is_file()
    ]
    if locks:
        run(
            [str(SCRIPTS / "check-dvc-lock.sh"), *locks],
            "a dvc.lock that dvc cannot read (see FAILED: above)",
        )


def written(path: str) -> bool:
    return path == "dvc.lock" or path.endswith("/dvc.lock") or path.startswith(WRITTEN)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", default="main")
    parser.add_argument("--run", action="append", default=[], metavar="YAML:STAGE")
    parser.add_argument(
        "--allow-reference", action="append", default=[], metavar="YAML:STAGE"
    )
    parser.add_argument("targets", nargs="*", metavar="DVC.YAML")
    args = parser.parse_args(argv)
    for entry in args.allow_reference:
        if not any(matches(entry, target) for target in args.run):
            parser.error(f"--allow-reference {entry} names no stage given with --run")
    return args


def relock(args: argparse.Namespace, ran: list[str]) -> None:
    """付け直しの段（始めの lock の確認の後）を順に行う。実行したステージを ran に足す（止まったときも、そこまでを報告するため）"""
    run(
        [str(SCRIPTS / "dvc-checkout.py")],
        "outputs did not come back as dvc.lock records them (see MISSING: above)",
    )
    run(
        ["dvc", "repro", "-s", *FORMS],
        "dvc repro -s of the build or a normalized form failed",
        quiet=True,
    )
    status = dvc_status(args.targets)
    run(
        [str(SCRIPTS / "check-references-shell.sh")],
        "references-shell.txt differs from the references devShell; follow the fix printed above",
    )
    for target in args.run:
        allow = [
            x
            for entry in args.allow_reference
            if matches(entry, target)
            for x in ("--allow-reference", entry)
        ]
        run(
            [str(SCRIPTS / "dvc-repro.py"), "-s", *allow, target],
            f"scripts/dvc-repro.py -s {target} stopped or failed",
        )
        ran.append(target)
    if ran:
        status = dvc_status(args.targets)
    check_locks(args.base)
    for key, kinds in left_changed(status):
        print(f"CHANGED: {key}: {', '.join(kinds)}")
    counts = counting.count_changed(status)
    for yaml_path in sorted(counts):
        print(
            f"COUNT: {yaml_path}: "
            + ", ".join(f"{k} {v}" for k, v in counts[yaml_path].items())
        )


def report(ran: list[str], before: dict[str, str]) -> None:
    for target in ran or ["none"]:
        print(f"RAN: {target}")
    after = worktree_files()
    for path in sorted(p for p in after if before.get(p) != after[p] or written(p)):
        print(f"COMMIT: {path}")


def stop_reason(e: Exception) -> str:
    if isinstance(e, subprocess.CalledProcessError):
        return f"{' '.join(e.cmd)} failed: {e.stderr}"
    return str(e)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    sys.stdout.reconfigure(line_buffering=True)
    try:
        root = git("rev-parse", "--show-toplevel").strip()
        # dvc の示すステージの key と --run はリポジトリの直下からのパスで、dvc の目標は cwd からの相対なので、直下でだけ動かす
        if Path.cwd().resolve() != Path(root).resolve():
            raise Stop(
                f"run from the root of the repository ({root}); the stages are paths from there"
            )
        before = worktree_files()
        # 作業ツリーを書き換える前なので、止まってもコミットするファイルは無い（読めない lock そのものを出さない）
        check_locks(args.base)
    except (Stop, subprocess.CalledProcessError) as e:
        print(f"STOP: {stop_reason(e)}", file=sys.stderr)
        return 1
    ran: list[str] = []
    try:
        relock(args, ran)
    except (Stop, subprocess.CalledProcessError) as e:
        report(ran, before)
        print(f"STOP: {stop_reason(e)}", file=sys.stderr)
        print(
            "NEXT: find the cause at the STOP: line (/dvc-relock). Rerun after fixing it; "
            "the COMMIT: files above are not committed yet",
            file=sys.stderr,
        )
        return 1
    report(ran, before)
    print(
        "NEXT: commit the COMMIT: files (git add them by name). Run a CHANGED: stage only if its result is needed, "
        "by naming it with --run"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
