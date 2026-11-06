"""固有値計算を解く順（ADR-0234 のケースの `node-order`）とスレッドの数（`threads`）だけを変えた同じ core の燃焼を比べ、点と中点の k、
ホウ素、燃焼度、出力の割合の差と、色の順の履歴のバイトの一致、燃焼ごとの固有値計算の数、外側の反復の回数、時間を出す。

使い方: python3 compare_core_threads.py --json <要約の出力先> <名前>=<直方体の並びの順の履歴>,<そのログ>,<色の順の履歴>,<そのログ>[,<色の順の履歴>,<そのログ>...]
        python3 compare_core_threads.py --report <要約の JSON>

- 履歴とログ、差の取り方は compare_core_guess.py と同じ（差は最初の色の順の燃焼 − 直方体の並びの順の燃焼）。スレッドの数はログの
  `thread pool threads=` の行から読む。最初の燃焼のログにその行が無ければ、ADR-0234 の前の bin/core（固有値計算は 1 スレッドの
  直方体の並びの順）の燃焼とし、順を `before`、スレッドの数を null にする
- `--json`: 名前ごとに `points`、`midpoints`（compare_core_guess.py と同じ形）と、`runs` の燃焼ごとの順（`lexicographic` か `red-black`）、
  スレッドの数、最初の色の順の履歴とバイトで同じか（`same_bytes`。直方体の並びの順は null）、固有値計算（span `nem`）の数、外側の反復の
  回数の和、固有値計算の時間の和（フィールド `seconds`）、core の全体の時間（`core seconds=`）、`per_solve` の固有値計算 1 回の外側の
  反復の回数と時間の [最小, 最大]。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import compare_core_acceleration
import compare_core_guess
import msgpack

POOL = re.compile(r"thread pool threads=(\d+) from=")

ORDERS = {
    "lexicographic": "直方体の並びの順",
    "red-black": "色の順",
    "before": "ADR-0234 の前の bin/core",
}


def threads_of(log: str) -> int | None:
    """ログのスレッドのプールの行のスレッドの数。行が無ければ None（ADR-0234 の前の bin/core）"""
    pools = POOL.findall(log)
    if len(pools) > 1:
        raise SystemExit("the log has more than one thread pool line")
    return int(pools[0]) if pools else None


def run_of(order: str, history: bytes, log: str, first: bytes | None) -> dict:
    core = msgpack.unpackb(history)
    threads = threads_of(log)
    if threads is None:
        if order != "lexicographic":
            raise SystemExit(
                "a red-black log needs a thread pool line (thread pool threads=)"
            )
        order = "before"
    return {
        "order": order,
        "threads": threads,
        "same_bytes": None if first is None else history == first,
        **compare_core_guess.run_of(core, log),
        "per_solve": compare_core_acceleration.per_solve(log),
    }


def compare(
    lexicographic: tuple[bytes, str], red_black: list[tuple[bytes, str]]
) -> dict:
    """最初の色の順の燃焼 − 直方体の並びの順の燃焼の差と、燃焼ごとの順、スレッドの数、バイトの一致、固有値計算の数、反復の回数、時間"""
    if not red_black:
        raise SystemExit("give at least one red-black run")
    base_history, base_log = lexicographic
    first_history, first_log = red_black[0]
    c = compare_core_guess.compare(
        msgpack.unpackb(base_history),
        base_log,
        msgpack.unpackb(first_history),
        first_log,
    )
    c.pop("runs")
    runs = [run_of("lexicographic", base_history, base_log, None)]
    runs += [run_of("red-black", h, text, first_history) for h, text in red_black]
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
        same = all(r["same_bytes"] for r in c["runs"] if r["same_bytes"] is not None)
        lines.append(
            f"{name}: 色の順の履歴はスレッドの数によらずバイトで同じ"
            if same
            else f"{name}: 色の順の履歴にバイトの違うものがある"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) < 4 or len(paths) % 2 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<lexicographic history>,<log>,<red-black history>,<log>[,...], got {arg}"
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
            pairs = [
                (paths[i].read_bytes(), paths[i + 1].read_text())
                for i in range(0, len(paths), 2)
            ]
            cases[name] = compare(pairs[0], pairs[1:])
            files += paths
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
