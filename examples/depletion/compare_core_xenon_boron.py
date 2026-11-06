"""平衡 Xe と臨界ホウ素の収束のさせ方（ADR-0219 の `xenon-boron`）だけを変えた同じ core の燃焼を比べ、点と中点の k、ホウ素、燃焼度、
出力の割合の差と、固有値計算の数、外側の反復の回数、時間を出す。

使い方: python3 compare_core_xenon_boron.py --json <要約の出力先> <名前>=<nested の履歴>,<そのログ>,<joint の履歴>,<そのログ> ...
        python3 compare_core_xenon_boron.py --report <要約の JSON>

- 履歴とログ、差の取り方は compare_core_guess.py と同じ（差は joint − nested）
- `--json`: 名前ごとに `points`、`midpoints`（compare_core_guess.py と同じ形）と、`runs` の `nested` と `joint` の固有値計算（span `nem`）の数、
  外側の反復の回数の和、固有値計算の時間の和（フィールド `seconds`）、core の全体の時間（`core seconds=`）。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import compare_core_guess
import msgpack


def compare(nested_core: dict, nested_log: str, core: dict, log: str) -> dict:
    """joint − nested の差と、2 つの燃焼の固有値計算の数、反復の回数、時間"""
    c = compare_core_guess.compare(nested_core, nested_log, core, log)
    runs = c.pop("runs")
    return {**c, "runs": {"nested": runs["flat"], "joint": runs["default"]}}


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        solves = c["points"] + c["midpoints"]
        k = max(abs(s["k"]) for s in solves)
        boron = max(abs(s["boron"]) for s in solves)
        burnup = max(p["burnup"] for p in c["points"])
        power = max(p["power"] for p in c["points"])
        n, j = c["runs"]["nested"], c["runs"]["joint"]
        lines.append(
            f"{name}: 点と中点の k の差の絶対値の最大 {k:.2e}、ホウ素 {boron:.2e} ppm、燃焼度の相対差の最大 {burnup:.2e}、"
            f"出力の割合の相対差の最大 {power:.2e}"
        )
        lines.append(
            f"{name}: 固有値計算 {n['solves']} → {j['solves']} 回（{j['solves'] / n['solves']:.3f} 倍）、外側の反復 "
            f"{n['iterations']:,} → {j['iterations']:,} 回（{j['iterations'] / n['iterations']:.3f} 倍）、固有値計算の時間 "
            f"{n['eigenvalue_seconds']:.1f} → {j['eigenvalue_seconds']:.1f} s、全体 {n['seconds']:.1f} → {j['seconds']:.1f} s"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) != 4 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<nested history>,<nested log>,<joint history>,<joint log>, got {arg}"
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
            nested_history, nested_log, history, log = paths
            cases[name] = compare(
                msgpack.unpackb(nested_history.read_bytes()),
                nested_log.read_text(),
                msgpack.unpackb(history.read_bytes()),
                log.read_text(),
            )
            files += paths
        summary = {
            "format": "xlc-compare-core-xenon-boron",
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
