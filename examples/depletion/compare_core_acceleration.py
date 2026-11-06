"""固有値計算の加速（ADR-0232 のケースの `acceleration`）だけを変えた同じ core の燃焼を比べ、点と中点の k、ホウ素、燃焼度、出力の割合の差と、
固有値計算の数、外側の反復の回数、時間と、固有値計算 1 回の外側の反復の回数と時間の範囲を出す。

使い方: python3 compare_core_acceleration.py --json <要約の出力先> <名前>=<加速なしの履歴>,<そのログ>,<CMFD の履歴>,<そのログ> ...
        python3 compare_core_acceleration.py --report <要約の JSON>

- 履歴とログ、差の取り方は compare_core_guess.py と同じ（差は CMFD − 加速なし）
- `--json`: 名前ごとに `points`、`midpoints`（compare_core_guess.py と同じ形）と、`runs` の `none` と `cmfd` の固有値計算（span `nem`）の数、
  外側の反復の回数の和、固有値計算の時間の和（フィールド `seconds`）、core の全体の時間（`core seconds=`）、`per_solve` の固有値計算 1 回の
  外側の反復の回数と時間の [最小, 最大]。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import compare_core_guess
import msgpack


def per_solve(log: str) -> dict:
    """固有値計算 1 回の外側の反復の回数と時間の [最小, 最大]"""
    nem = [(int(i), float(s)) for i, s in compare_core_guess.NEM.findall(log)]
    iterations = [i for i, _ in nem]
    seconds = [s for _, s in nem]
    return {
        "iterations": [min(iterations), max(iterations)],
        "seconds": [min(seconds), max(seconds)],
    }


def compare(plain_core: dict, plain_log: str, core: dict, log: str) -> dict:
    """CMFD − 加速なしの差と、2 つの燃焼の固有値計算の数、反復の回数、時間"""
    c = compare_core_guess.compare(plain_core, plain_log, core, log)
    runs = c.pop("runs")
    return {
        **c,
        "runs": {
            "none": {**runs["flat"], "per_solve": per_solve(plain_log)},
            "cmfd": {**runs["default"], "per_solve": per_solve(log)},
        },
    }


def span(values: list, digits: str) -> str:
    low, high = values
    return f"{low:{digits}}〜{high:{digits}}"


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        solves = c["points"] + c["midpoints"]
        k = max(abs(s["k"]) for s in solves)
        boron = max(abs(s["boron"]) for s in solves)
        burnup = max(p["burnup"] for p in c["points"])
        power = max(p["power"] for p in c["points"])
        n, a = c["runs"]["none"], c["runs"]["cmfd"]
        lines.append(
            f"{name}: 点と中点の k の差の絶対値の最大 {k:.2e}、ホウ素 {boron:.2e} ppm、燃焼度の相対差の最大 {burnup:.2e}、"
            f"出力の割合の相対差の最大 {power:.2e}"
        )
        lines.append(
            f"{name}: 固有値計算 {n['solves']} → {a['solves']} 回、外側の反復 {n['iterations']:,} → {a['iterations']:,} 回"
            f"（{a['iterations'] / n['iterations']:.3f} 倍）、固有値計算の時間 {n['eigenvalue_seconds']:.1f} → "
            f"{a['eigenvalue_seconds']:.1f} s（{a['eigenvalue_seconds'] / n['eigenvalue_seconds']:.3f} 倍）、全体 "
            f"{n['seconds']:.1f} → {a['seconds']:.1f} s"
        )
        lines.append(
            f"{name}: 固有値計算 1 回の外側の反復 {span(n['per_solve']['iterations'], ',')} → "
            f"{span(a['per_solve']['iterations'], ',')} 回、時間 {span(n['per_solve']['seconds'], '.3f')} → "
            f"{span(a['per_solve']['seconds'], '.3f')} s"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) != 4 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<unaccelerated history>,<unaccelerated log>,<CMFD history>,<CMFD log>, got {arg}"
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
            plain_history, plain_log, history, log = paths
            cases[name] = compare(
                msgpack.unpackb(plain_history.read_bytes()),
                plain_log.read_text(),
                msgpack.unpackb(history.read_bytes()),
                log.read_text(),
            )
            files += paths
        summary = {
            "format": "xlc-compare-core-acceleration",
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
