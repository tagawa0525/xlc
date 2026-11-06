"""固有値計算の初期値（ADR-0231）だけを変えた同じ core の燃焼を比べ、点と中点の k、ホウ素、燃焼度、出力の割合の差と、固有値計算の数、
外側の反復の回数、時間を出す。

使い方: python3 compare_core_guess.py --json <要約の出力先> <名前>=<平らな初期値の履歴>,<そのログ>,<既定の初期値の履歴>,<そのログ> ...
        python3 compare_core_guess.py --report <要約の JSON>

- 履歴は core の燃焼の履歴（msgpack）、ログはその標準エラー（tracing）。2 つの燃焼は、同じ燃焼点（EFPD と日数）、状態、定格の出力密度、
  軸方向のメッシュ、集合体を持つこと
- `--json`: 名前ごとに次を書く
  - `points`: 燃焼点ごとの EFPD、日数、状態と、既定 − 平らな初期値の k とホウ素の差、ノードの燃焼度の相対差の絶対値の最大、出力の
    割合の相対差の絶対値の最大（核分裂のないノードを除く）
  - `midpoints`: ログの中点の行（`depletion midpoint`）ごとの、既定 − 平らな初期値の k とホウ素の差
  - `runs`: `flat` と `default` の、固有値計算（span `nem`）の数、外側の反復の回数の和、固有値計算の時間の和（フィールド `seconds`）、
    core の全体の時間（`core seconds=`）、燃焼点と中点の固有値計算の回数の和（履歴の `solves` と中点の行の `solves`）
  `inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack

# span `nem` の閉じる行のフィールド
NEM = re.compile(r"nem\{iterations=(\d+) seconds=([0-9.eE+-]+)\}")
TOTAL = re.compile(r"core seconds=([0-9.eE+-]+)")
MIDPOINT = re.compile(
    r"depletion midpoint k=([0-9.eE+-]+) boron=([0-9.eE+-]+) solves=(\d+)"
)
# 2 つの燃焼で同じであるべき履歴の項目
SHARED = ("specific_power", "states", "axial_mesh", "assemblies", "heavy_metal")


def flat(rows: list[list[float]]) -> list[float]:
    return [v for row in rows for v in row]


def largest_relative(a: list[float], b: list[float]) -> float:
    """b を基準にした相対差の絶対値の最大。b が 0 の値は、a も 0 であること（核分裂のないノード）"""
    if len(a) != len(b):
        raise SystemExit(f"{len(a)} and {len(b)} values")
    worst = 0.0
    for x, y in zip(a, b):
        if y == 0.0:
            if x != 0.0:
                raise SystemExit(f"{x} against 0")
            continue
        worst = max(worst, abs(x - y) / abs(y))
    return worst


def run_of(core: dict, log: str) -> dict:
    """燃焼の固有値計算の数、外側の反復の回数と時間の和、全体の時間、燃焼点と中点の固有値計算の回数の和"""
    nem = [(int(i), float(s)) for i, s in NEM.findall(log)]
    total = TOTAL.findall(log)
    if not nem or len(total) != 1:
        raise SystemExit(
            "the log needs the nem spans and one total time (core seconds=)"
        )
    midpoints = [int(m[2]) for m in MIDPOINT.findall(log)]
    solves = sum(p["solves"] for p in core["points"]) + sum(midpoints)
    if solves != len(nem):
        raise SystemExit(
            f"the history and the midpoint lines count {solves} solves, the log {len(nem)} nem spans"
        )
    return {
        "solves": len(nem),
        "iterations": sum(i for i, _ in nem),
        "eigenvalue_seconds": sum(s for _, s in nem),
        "seconds": float(total[0]),
    }


def compare(flat_core: dict, flat_log: str, core: dict, log: str) -> dict:
    for key in SHARED:
        if flat_core.get(key) != core.get(key):
            raise SystemExit(f"the histories differ in {key}")
    a, b = core["points"], flat_core["points"]
    if [(p["efpd"], p["days"], p["state"]) for p in a] != [
        (p["efpd"], p["days"], p["state"]) for p in b
    ]:
        raise SystemExit("the histories have different points")
    points = [
        {
            "efpd": q["efpd"],
            "days": q["days"],
            "state": q["state"],
            "k": p["k"] - q["k"],
            "boron": p["boron"] - q["boron"],
            "burnup": largest_relative(flat(p["burnups"]), flat(q["burnups"])),
            "power": largest_relative(flat(p["powers"]), flat(q["powers"])),
        }
        for p, q in zip(a, b)
    ]
    mids, flat_mids = MIDPOINT.findall(log), MIDPOINT.findall(flat_log)
    if len(mids) != len(flat_mids):
        raise SystemExit(f"{len(mids)} and {len(flat_mids)} midpoints")
    midpoints = [
        {"k": float(m[0]) - float(f[0]), "boron": float(m[1]) - float(f[1])}
        for m, f in zip(mids, flat_mids)
    ]
    return {
        "points": points,
        "midpoints": midpoints,
        "runs": {"flat": run_of(flat_core, flat_log), "default": run_of(core, log)},
    }


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        solves = c["points"] + c["midpoints"]
        k = max(abs(s["k"]) for s in solves)
        boron = max(abs(s["boron"]) for s in solves)
        burnup = max(p["burnup"] for p in c["points"])
        power = max(p["power"] for p in c["points"])
        f, d = c["runs"]["flat"], c["runs"]["default"]
        lines.append(
            f"{name}: 点と中点の k の差の絶対値の最大 {k:.2e}、ホウ素 {boron:.2e} ppm、燃焼度の相対差の最大 {burnup:.2e}、"
            f"出力の割合の相対差の最大 {power:.2e}"
        )
        lines.append(
            f"{name}: 固有値計算 {f['solves']} → {d['solves']} 回、外側の反復 {f['iterations']:,} → {d['iterations']:,} 回"
            f"（{d['iterations'] / f['iterations']:.3f} 倍）、固有値計算の時間 {f['eigenvalue_seconds']:.1f} → "
            f"{d['eigenvalue_seconds']:.1f} s、全体 {f['seconds']:.1f} → {d['seconds']:.1f} s"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) != 4 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<flat history>,<flat log>,<history>,<log>, got {arg}"
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
            flat_history, flat_log, history, log = paths
            cases[name] = compare(
                msgpack.unpackb(flat_history.read_bytes()),
                flat_log.read_text(),
                msgpack.unpackb(history.read_bytes()),
                log.read_text(),
            )
            files += paths
        summary = {
            "format": "xlc-compare-core-guess",
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
