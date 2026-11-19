#!/usr/bin/env python3
"""dvc repro の前に走りうるステージを数え、上限を超えるか、許していない参照解のステージが入っていれば、何も実行せずに止まる（ADR-0257）。

使い方（default の devShell で、リポジトリの直下から）:
    scripts/dvc-repro.py --max-stages <N> [--allow-reference <yaml>:<ステージか foreach の組>]... [--count-only]
                         [--no-run-cache] (-P | <dvc repro の目標>...)

- 走る範囲は `dvc repro --dry -f <目標>` が並べるステージ（目標とその上流）で、cmd もそこから読む
- 走りうるステージは、範囲のうち `dvc status --json --with-deps <目標>` が changed と示すものと、`dvc dag --dot` でその下流に
  あるもの。dvc status と dvc repro --dry は上流が出力を書き換えた後に走る下流を示さないので、下流を足す（上流が同じ出力を
  書けば下流は飛ぶので、数は上限）
- references-shell のステージは always_changed でいつも changed だが、書く references-shell.txt が references の devShell と
  同じなら下流（参照解）は飛ぶ。scripts/check-references-shell.sh が通るときだけ、その下流を足さない
- 参照解のステージは、cmd を references の devShell（`nix develop <flake>#references`）で包んだもの（NJOY と OpenMC）。
  `--allow-reference` に `<yaml>:<ステージ>` か `<yaml>:<組>`（foreach の組の全部）で挙げたものだけを実行してよい
- 一覧（ステージごとに changed の種類か、どの changed の下流か）を出し、通れば `dvc repro <目標>` を同じプロセスで実行する
  （終了コードは dvc repro のもの）。`--count-only` は数えて終わる（通れば 0、止まれば 1）
- `-s` は受けない（上流を飛ばす数秒の再実行は手元で実行する。ADR-0184）。ほかの dvc repro のオプションも、走るステージの
  数え方を変えうるので受けない
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

REFERENCE_CMD = re.compile(r"\bnix develop \S*#references\b")
DRY_STAGE = re.compile(r"Running stage '(.+)':")
DOT_EDGE = re.compile(r'"(.+)" -> "(.+)";')
REFERENCES_SHELL = "references-shell"
CHECK_REFERENCES_SHELL = Path(__file__).resolve().parent / "check-references-shell.sh"


class ReproError(Exception):
    pass


def parse_dry_run(text: str) -> dict[str, list[str]]:
    """`dvc repro --dry -f` の出力から、ステージ → cmd の行（実行の順）"""
    cmds: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if m := DRY_STAGE.fullmatch(line):
            current = m.group(1)
            cmds[current] = []
        elif line.startswith("> "):
            if current is None:
                raise ReproError(
                    f"command without a stage in dvc repro --dry: {line!r}"
                )
            cmds[current].append(line[2:])
    return cmds


def parse_dot(text: str) -> set[tuple[str, str]]:
    """`dvc dag --dot` の出力から、(上流, 下流) の辺"""
    return {
        (m.group(1), m.group(2))
        for line in text.splitlines()
        if (m := DOT_EDGE.fullmatch(line.strip()))
    }


def change_kinds(entries: list) -> list[str]:
    """dvc status --json のステージの項目から、変更の種類（{種類: 中身} の辞書か、中身の無い種類の文字列）"""
    return sorted(
        {kind for e in entries for kind in ([e] if isinstance(e, str) else e)}
    )


def stage_name(key: str) -> str:
    """`<yaml>:<ステージ>` のステージの部分（直下の dvc.yaml のステージはコロンを持たない）"""
    return key.rsplit(":", 1)[-1]


def is_reference(cmd_lines: list[str]) -> bool:
    return any(REFERENCE_CMD.search(line) for line in cmd_lines)


def allows(entry: str, key: str) -> bool:
    """`--allow-reference` の 1 つが、ステージの key（dvc の示す形）か、その foreach の組を指すか"""
    return key == entry or key.split("@", 1)[0] == entry


@dataclass
class Count:
    runs: dict[str, str]  # 走りうるステージ → 理由（実行の順）
    references: list[str]  # runs のうち参照解のステージ
    refused: list[str]  # references のうち許していないもの
    outside: list[
        str
    ]  # dvc status が示したが、dvc repro の範囲に無いもの（.dvc ファイルなど）


def count(
    status: dict,
    edges: set[tuple[str, str]],
    cmds: dict[str, list[str]],
    shell_matches: bool,
    allow: list[str],
) -> Count:
    """走りうるステージを数える。shell_matches は references-shell.txt が references の devShell と同じか"""
    scope = list(cmds)
    in_scope = set(scope)
    for entry in allow:
        if not any(allows(entry, key) for key in scope):
            raise ReproError(
                f"--allow-reference {entry} names no stage of this dvc repro "
                "(write it as <yaml>:<stage> or <yaml>:<foreach group>, as dvc status shows)"
            )
    children: dict[str, set[str]] = {}
    for up, down in edges:
        if up in in_scope and down in in_scope:
            children.setdefault(up, set()).add(down)

    reasons = {}
    sources = []
    for key, entries in status.items():
        if key not in in_scope:
            continue
        kinds = change_kinds(entries)
        reasons[key] = "changed: " + ", ".join(kinds)
        if (
            stage_name(key) == REFERENCES_SHELL
            and kinds == ["always changed"]
            and shell_matches
        ):
            continue
        sources.append(key)
    # 実行の順にたどり、下流の理由には最初に届いた changed のステージを書く
    sources.sort(key=scope.index)
    for source in sources:
        stack = [source]
        while stack:
            for child in sorted(children.get(stack.pop(), ())):
                if child not in reasons:
                    reasons[child] = f"downstream of {source}"
                    stack.append(child)

    runs = {key: reasons[key] for key in scope if key in reasons}
    references = [key for key in runs if is_reference(cmds[key])]
    refused = [
        key for key in references if not any(allows(entry, key) for entry in allow)
    ]
    outside = sorted(key for key in status if key not in in_scope)
    return Count(runs, references, refused, outside)


def dvc(*args: str) -> str:
    start = time.monotonic()
    result = subprocess.run(["dvc", *args], capture_output=True, text=True, check=False)
    print(f"TIME: dvc {' '.join(args)} {time.monotonic() - start:.1f} s", flush=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise ReproError(
            f"dvc {' '.join(args)} failed with exit code {result.returncode}"
        )
    return result.stdout


def references_shell_matches() -> bool:
    """scripts/check-references-shell.sh が通るか（食い違えば、references-shell の下流は再実行される）"""
    result = subprocess.run([str(CHECK_REFERENCES_SHELL)], check=False)
    return result.returncode == 0


def report(c: Count, max_stages: int) -> list[str]:
    """一覧を出し、止まる理由を返す（空なら通る）"""
    for key, reason in c.runs.items():
        mark = " [reference]" if key in c.references else ""
        print(f"RUN: {key}{mark} ({reason})")
    for key in c.outside:
        print(f"NOTE: {key} is changed but is not a stage this dvc repro runs")
    print(
        f"COUNT: stages={len(c.runs)} max={max_stages} references={len(c.references)}"
    )
    problems = []
    if len(c.runs) > max_stages:
        problems.append(
            f"{len(c.runs)} stages may run, more than --max-stages {max_stages}"
        )
    if c.refused:
        problems.append(
            "reference stages not allowed by --allow-reference: " + ", ".join(c.refused)
        )
    return problems


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--max-stages", type=int, required=True)
    parser.add_argument("--allow-reference", action="append", default=[])
    parser.add_argument("--count-only", action="store_true")
    parser.add_argument("--no-run-cache", action="store_true")
    parser.add_argument("-P", "--all-pipelines", action="store_true")
    # -s を受けないことを、argparse の「unrecognized arguments」でなく理由とともに示す
    parser.add_argument(
        "-s", "--single-item", action="store_true", help=argparse.SUPPRESS
    )
    parser.add_argument("targets", nargs="*")
    args = parser.parse_args(argv)
    if args.single_item:
        parser.error(
            "-s is not counted here: run the seconds-long reruns that skip the upstream "
            "locally with dvc repro -s (/dvc-relock, ADR-0184)"
        )
    if args.max_stages < 0:
        parser.error("--max-stages must not be negative")
    if args.all_pipelines == bool(args.targets):
        parser.error("give either -P or the targets of dvc repro")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    # Slurm のログは stdout と stderr を 1 つのファイルにまとめるので、行ごとに書き出して STOP: を一覧の後に置く
    sys.stdout.reconfigure(line_buffering=True)
    selection = ["-P"] if args.all_pipelines else args.targets
    try:
        root = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        # dvc の示すステージの key（<yaml>:<ステージ>）と --allow-reference はリポジトリの直下からのパスで、目標は
        # cwd からの相対なので、直下に移ると別の dvc.yaml を指しうる。移らずに止める
        if Path.cwd().resolve() != Path(root).resolve():
            raise ReproError(
                f"run from the root of the repository ({root}); the targets are paths from there"
            )
        cmds = parse_dry_run(dvc("repro", "--dry", "-f", *selection))
        if not cmds:
            raise ReproError(
                f"dvc repro --dry -f {' '.join(selection)} listed no stage"
            )
        edges = set()
        for target in [None] if args.all_pipelines else args.targets:
            edges |= parse_dot(
                dvc("dag", "--dot", *([] if target is None else [target]))
            )
        status = json.loads(
            dvc(
                "status",
                "--json",
                *([] if args.all_pipelines else ["--with-deps", *args.targets]),
            )
        )
        shell_matches = True
        if any(stage_name(key) == REFERENCES_SHELL and key in cmds for key in status):
            shell_matches = references_shell_matches()
            print(
                f"REFERENCES-SHELL: matches the references devShell: {'yes' if shell_matches else 'no'}"
            )
        c = count(status, edges, cmds, shell_matches, args.allow_reference)
    except (ReproError, subprocess.CalledProcessError, json.JSONDecodeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    problems = report(c, args.max_stages)
    if problems:
        for problem in problems:
            print(f"STOP: {problem}", file=sys.stderr)
        print(
            "NEXT: nothing ran. Find out why each listed stage is changed (a deps change that does not "
            "alter the result, a stale normalized form under python-ast/, toml-values/, binary-forms/ "
            "that /dvc-relock refreshes locally), or narrow the targets. Raise --max-stages or add "
            "--allow-reference only for stages the job is meant to run",
            file=sys.stderr,
        )
        return 1
    if args.count_only:
        return 0
    command = [
        "dvc",
        "repro",
        *(["--no-run-cache"] if args.no_run_cache else []),
        *selection,
    ]
    print(f"RUNNING: {' '.join(command)}", flush=True)
    os.execvp("dvc", command)


if __name__ == "__main__":
    sys.exit(main())
