"""固有値計算を解く順（ADR-0234 のケースの `node-order`）とスレッドの数（`threads`）だけを変えた同じ core の燃焼を比べ、点と中点の k、
ホウ素、燃焼度、出力の割合の差と、色の順の履歴のバイトの一致、燃焼ごとの固有値計算の数、外側の反復の回数、時間を出す。

使い方: python3 compare_core_threads.py --json <要約の出力先> <名前>=<直方体の並びの順の core の入力>,<色の順の core の入力>[,<色の順の core の入力>...]
        python3 compare_core_threads.py --report <要約の JSON>

- 入力から読む履歴とログ、差の取り方は compare_core_pair.py と同じ（差は最初の色の順の燃焼 − 直方体の並びの順の燃焼）。燃焼の順は
  入力の `node-order` で、最初の入力が `lexicographic`、ほかが `red-black`（書かなければ core の既定の `red-black`）でなければ止まる。
  スレッドの数はログの `thread pool threads=` の行から読み、行が無いか、入力の `threads` と違えば止まる
- `--json`: 名前ごとに `points`、`midpoints`（compare_core_pair.py と同じ形）と、`runs` の燃焼ごとの順（`lexicographic` か `red-black`）、
  スレッドの数、最初の色の順の履歴とバイトで同じか（`same_bytes`。直方体の並びの順と最初の色の順は null）、固有値計算（span `nem`）の数、外側の反復の
  回数の和、固有値計算の時間の和（フィールド `seconds`）、core の全体の時間（`core seconds=`）、`per_solve` の固有値計算 1 回の外側の
  反復の回数と時間の [最小, 最大]。`inputs` は読んだファイル（燃焼ごとに入力、履歴、ログの順）のパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import compare_core_pair
import msgpack
from core_history import core_run, setting

POOL = re.compile(r"thread pool threads=(\d+) from=")

ORDERS = {
    "lexicographic": "直方体の並びの順",
    "red-black": "色の順",
}


def threads_of(log: str) -> int:
    """ログのスレッドのプールの行のスレッドの数"""
    pools = POOL.findall(log)
    if len(pools) != 1:
        raise SystemExit(
            f"the log needs one thread pool line (thread pool threads=), found {len(pools)}"
        )
    return int(pools[0])


def run_of(order: str, history: bytes, log: str, first: bytes | None) -> dict:
    core = msgpack.unpackb(history)
    return {
        "order": order,
        "threads": threads_of(log),
        "same_bytes": None if first is None else history == first,
        **compare_core_pair.run_of(core, log),
        "per_solve": compare_core_pair.per_solve(log),
    }


def thread_runs(sources: list[Path]) -> tuple[list[tuple[bytes, str]], list[Path]]:
    """入力ごとの履歴のバイトとログと、読んだファイル。最初の入力の `node-order` が `lexicographic`、ほかが `red-black` でなければ、
    また入力の `threads` がログのスレッドの数と違えば止まる"""
    runs = []
    files = []
    for i, source in enumerate(sources):
        run, history, log_path = core_run(source)
        want = "lexicographic" if i == 0 else "red-black"
        got = setting(run, None, "node-order", "red-black")
        if got != want:
            raise SystemExit(
                f"{source} sets node-order = {got!r}; input {i + 1} must set {want!r}"
            )
        log = log_path.read_text()
        if "threads" in run and run["threads"] != threads_of(log):
            raise SystemExit(
                f"{source} sets threads = {run['threads']}, but {log_path} ran {threads_of(log)}"
            )
        runs.append((history.read_bytes(), log))
        files += [source, history, log_path]
    return runs, files


def compare(
    lexicographic: tuple[bytes, str], red_black: list[tuple[bytes, str]]
) -> dict:
    """最初の色の順の燃焼 − 直方体の並びの順の燃焼の差と、燃焼ごとの順、スレッドの数、バイトの一致、固有値計算の数、反復の回数、時間"""
    if not red_black:
        raise SystemExit("give at least one red-black run")
    base_history, base_log = lexicographic
    first_history, first_log = red_black[0]
    c = compare_core_pair.differences(
        msgpack.unpackb(base_history),
        base_log,
        msgpack.unpackb(first_history),
        first_log,
    )
    runs = [run_of("lexicographic", base_history, base_log, None)]
    # 最初の色の順の燃焼は比べる基準なので、バイトの一致を書かない
    runs.append(run_of("red-black", first_history, first_log, None))
    runs += [run_of("red-black", h, text, first_history) for h, text in red_black[1:]]
    return {**c, "runs": runs}


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        solves = c["points"] + c["midpoints"]
        k = max(abs(s["k"]) for s in solves)
        boron = max(abs(s["boron"]) for s in solves)
        burnup = max(p["burnup"] for p in c["points"])
        power = max(p["power"] for p in c["points"])
        lines.append(
            f"{name}: 点と中点の k の差の絶対値の最大 {k:.2e}、ホウ素 {boron:.2e} ppm、燃焼度の相対差の最大 {burnup:.2e}、"
            f"出力の割合の相対差の最大 {power:.2e}"
        )
        for r in c["runs"]:
            threads = "" if r["threads"] is None else f" {r['threads']} スレッド"
            lines.append(
                f"{name}: {ORDERS[r['order']]}{threads}: 固有値計算 {r['solves']} 回、外側の反復 {r['iterations']:,} 回、"
                f"固有値計算の時間 {r['eigenvalue_seconds']:.1f} s、全体 {r['seconds']:.1f} s"
            )
        compared = [r["same_bytes"] for r in c["runs"] if r["same_bytes"] is not None]
        if compared:
            lines.append(
                f"{name}: 色の順の履歴はスレッドの数によらずバイトで同じ"
                if all(compared)
                else f"{name}: 色の順の履歴にバイトの違うものがある"
            )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    parts = files.split(",")
    paths = [Path(p) for p in parts]
    if not name or len(paths) < 2 or not all(parts):
        raise SystemExit(
            f"give <name>=<lexicographic core input>,<red-black core input>[,...], got {arg}"
        )
    return name, paths


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 3 and args[0] == "--json":
        cases = {}
        files = []
        for name, paths in map(parse, args[2:]):
            if name in cases:
                raise SystemExit(f"a second case named {name}")
            runs, read = thread_runs(paths)
            cases[name] = compare(runs[0], runs[1:])
            files += read
        summary = {
            "format": "xlc-compare-core-threads",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            "cases": cases,
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
