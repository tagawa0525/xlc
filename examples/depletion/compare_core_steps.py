"""同じサイクルをステップの幅を変えて燃焼させた core の燃焼（ADR-0190、ADR-0218）を比べ、ステップの大きさの誤差と、臨界ホウ素の
探索の中点のホウ素（ADR-0191 の探索と案 S1）の差と、固有値計算の回数を出す。

使い方: python3 compare_core_steps.py --json <要約の出力先> <分割>=<core の履歴（msgpack）>,<その標準エラー（ログ）> ...
        python3 compare_core_steps.py --report <要約の JSON>

- `<分割>` は、基準のサイクルの出力のあるステップ（`efpd`）を等しい EFPD の何個に分けたか（1、2、4 など）。分割 1 の履歴を基準の
  サイクルとし、ほかの履歴は基準の燃焼点をすべて同じ EFPD と日数（日数は 1e-6 日に丸めて）で持つこと。崩壊だけのステップは分けない。
  燃焼の方式（`method_name`: Xe を時間で追う燃焼は中点のホウ素 `search` か `start`、Xe を平衡にした燃焼（#0200）は `equilibrium`）
  ごとに分割 1 があること（分割の組は方式ごとに違ってよい。`splits` は方式 → 分割）
- `--json`: 履歴がどれも同じ定格の出力密度、状態、軸方向のメッシュ、集合体、重金属の質量を持つことを確かめ、次を書く。
  - `steps`: 中点のホウ素ごとに、基準の燃焼点ごとの分割ごとのホウ素濃度、k、ノードの燃焼度と出力の割合と、隣り合う分割（1 と 2、
    2 と 4 …）の差（ホウ素は分割の粗い方 − 細かい方、燃焼度はノードの差の絶対値の最大、出力の割合はノードの相対差の絶対値の最大）と、
    その比（粗い組の差 / 細かい組の差。中点の形が 2 次なら 4 に近い。ホウ素は符号を保ち、負なら隣り合う差の符号が変わった。細かい組の差が
    0 なら null）。分割は倍ずつで、組にした点は基準と同じ状態と出力の割合で解いたものであること
  - `midpoint`: `search` と `start` の両方があれば、両方にある分割ごと、基準の燃焼点ごとの `start` − `search` の差（ホウ素、燃焼度の最大、出力の
    割合の相対差の最大）
  - `solves`: 履歴ごとの燃焼点の数、点の固有値計算の回数（履歴の `solves` の和）、中点の固有値計算の回数（ログの `depletion midpoint`
    の行の `solves` の和）とその和
  - `axial`: 履歴ごとの、すべての燃焼点の軸方向の出力の偏り（炉心の高さの上半分のノードの出力の割合の和 − 下半分の和。割合の炉心の和は
    1。軸方向のメッシュは高さの半分に境界を持つこと）
  `inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from itertools import pairwise
from pathlib import Path

import msgpack

# core のログ（crates/core/src/depletion.rs の tracing）の中点の行の、固有値計算の回数
MIDPOINT_SOLVES = re.compile(r"\bsolves=(\d+)\b")

# 履歴の、比べる燃焼の共通の項目（ミクロ燃焼の追う核種 `nuclides` は、ミクロ燃焼しない履歴には無い）
SHARED = (
    "specific_power",
    "states",
    "axial_mesh",
    "assemblies",
    "heavy_metal",
    "nuclides",
)


def method_name(core: dict) -> str:
    """履歴の燃焼の方式の名前: Xe を時間で追う燃焼（`xenon` が `tracked` か無い）は `midpoint_boron`（`search` か `start`）、Xe を
    平衡にした燃焼は `equilibrium`（中点でも探す燃焼だけ）"""
    midpoint = core.get("midpoint_boron")
    xenon = core.get("xenon", "tracked")
    if midpoint not in ("search", "start") or xenon not in ("tracked", "equilibrium"):
        raise SystemExit(
            f"unexpected midpoint_boron {midpoint} and xenon {xenon} of a history"
        )
    if xenon == "tracked":
        return midpoint
    if midpoint != "search":
        raise SystemExit("equilibrium xenon is compared with the midpoint search only")
    return "equilibrium"


def axial_offset(core: dict, point: dict) -> float:
    """点の軸方向の出力の偏り（上半分のノードの出力の割合の和 − 下半分の和）"""
    mesh = core["axial_mesh"]
    half = 0.5 * (mesh[0] + mesh[-1])
    if half not in mesh:
        raise SystemExit(f"the axial mesh {mesh} has no boundary at half the height")
    upper = mesh.index(half)
    return sum(sum(row[upper:]) - sum(row[:upper]) for row in point["powers"])


def key(point: dict) -> tuple[float, float]:
    """燃焼点を分割の間で組にする鍵（EFPD と、崩壊だけのステップの日数を含む日数。日数は分けた時間の和の丸めを除くため 1e-6 日に丸める）"""
    return (point["efpd"], round(point["days"], 6))


def midpoint_solves(log: str) -> list[int]:
    """ログの中点の行の固有値計算の回数（中点の順）"""
    found = []
    for line in log.splitlines():
        if "depletion midpoint" in line:
            m = MIDPOINT_SOLVES.search(line)
            if m is None:
                raise SystemExit(f"a midpoint line without solves: {line}")
            found.append(int(m[1]))
    return found


def flat(rows: list[list[float]]) -> list[float]:
    return [v for row in rows for v in row]


def differences(coarse: dict, fine: dict) -> dict:
    """同じ燃焼点の、粗い分割 − 細かい分割の差"""
    burnup = max(
        abs(a - b) for a, b in zip(flat(coarse["burnups"]), flat(fine["burnups"]))
    )
    pairs = list(zip(flat(coarse["powers"]), flat(fine["powers"])))
    # 核分裂のないノード（両方 0）だけを除く。片方だけ 0 のノードは相対差を持たない
    if any((a > 0.0) != (b > 0.0) for a, b in pairs):
        raise SystemExit(
            f"the point at {fine['efpd']} EFPD has a node with power in one split only"
        )
    power = max((abs(a - b) / b for a, b in pairs if b > 0.0), default=0.0)
    return {
        "boron": coarse["boron"] - fine["boron"],
        "k": coarse["k"] - fine["k"],
        "burnup": burnup,
        "power": power,
    }


def ratio(coarse: float, fine: float) -> float | None:
    """粗い組の差 / 細かい組の差。符号を保つ（隣り合う差の符号が変われば負で、漸近の範囲の外）"""
    return None if fine == 0.0 else coarse / fine


def check(runs: dict[tuple[str, int], dict]) -> None:
    first = next(iter(runs.values()))
    for (method, split), core in runs.items():
        name = f"{method} {split}"
        if core.get("format") != "xlc-core-history":
            raise SystemExit(f"{name}: unexpected format {core.get('format')}")
        for item in SHARED:
            if core.get(item) != first.get(item):
                raise SystemExit(f"{name}: {item} differs from the other histories")
    for method in {m for m, _ in runs}:
        splits = sorted(s for m, s in runs if m == method)
        if splits[0] != 1:
            raise SystemExit(f"{method}: give the base cycle (split 1)")
        # 比を次数として読む（2 次なら 4）のは、分割が倍ずつのときだけ
        if any(b != 2 * a for a, b in pairwise(splits)):
            raise SystemExit(f"{method}: the splits {splits} must double")


def summarize(
    runs: dict[tuple[str, int], dict], logs: dict[tuple[str, int], str]
) -> dict:
    """分割と中点のホウ素ごとの履歴を比べた要約（モジュールの説明の `--json`）。"""
    check(runs)
    methods = sorted({m for m, _ in runs})
    splits = {m: sorted(s for n, s in runs if n == m) for m in methods}
    base = runs[(methods[0], 1)]["points"]
    # 分割 n の履歴は、燃焼点 0、崩壊だけのステップの点、基準の出力のあるステップの n 倍の点を持つ（分割の名前と履歴を取り違えると、
    # 基準の点はそろっても点の数が合わない）
    base_steps = sum(1 for a, b in pairwise(base) if b["efpd"] != a["efpd"])
    decays = len(base) - 1 - base_steps
    for (method, split), core in runs.items():
        expected = 1 + decays + split * base_steps
        if len(core["points"]) != expected:
            raise SystemExit(
                f"{method} {split}: {len(core['points'])} points, but split {split} of the base cycle has {expected}"
            )
    # 基準の燃焼点を、分割ごとの履歴の点に組にする
    matched: dict[tuple[str, int], list[dict]] = {}
    for (method, split), core in runs.items():
        by_key = {key(p): p for p in core["points"]}
        if len(by_key) != len(core["points"]):
            raise SystemExit(f"{method} {split}: two points share the EFPD and days")
        missing = [key(p) for p in base if key(p) not in by_key]
        if missing:
            raise SystemExit(f"{method} {split}: no points at (EFPD, days) {missing}")
        matched[(method, split)] = [by_key[key(p)] for p in base]
        # 同じ点は同じ状態と出力の割合で解いたもの（分けたケースの状態の書き誤りを、ステップの大きさの誤差に混ぜない）
        for p, q in zip(base, matched[(method, split)]):
            if (q["state"], q.get("power")) != (p["state"], p.get("power")):
                raise SystemExit(
                    f"{method} {split}: the point at {p['efpd']} EFPD is solved in {q['state']} at power "
                    f"{q.get('power')}, the base in {p['state']} at power {p.get('power')}"
                )
    steps = {}
    for method in methods:
        points = []
        for i, p in enumerate(base):
            at = {s: matched[(method, s)][i] for s in splits[method]}
            pairs = [
                {"splits": [a, b], **differences(at[a], at[b])}
                for a, b in pairwise(splits[method])
            ]
            points.append(
                {
                    "efpd": p["efpd"],
                    "days": p["days"],
                    "power": p.get("power"),
                    "state": p["state"],
                    "burnup": p["burnup"],
                    "boron": [at[s]["boron"] for s in splits[method]],
                    "k": [at[s]["k"] for s in splits[method]],
                    "differences": pairs,
                    "ratios": [
                        {
                            "splits": [a["splits"][0], b["splits"][1]],
                            **{
                                q: ratio(a[q], b[q])
                                for q in ("boron", "burnup", "power")
                            },
                        }
                        for a, b in pairwise(pairs)
                    ],
                }
            )
        steps[method] = points
    midpoint = None
    if {"search", "start"} <= set(methods):
        midpoint = {
            str(split): [
                {
                    "efpd": p["efpd"],
                    "days": p["days"],
                    **differences(
                        matched[("start", split)][i], matched[("search", split)][i]
                    ),
                }
                for i, p in enumerate(base)
            ]
            # 両方にある分割
            for split in splits["search"]
            if split in splits["start"]
        }
    solves = []
    for (method, split), core in sorted(runs.items()):
        mids = midpoint_solves(logs[(method, split)])
        run_steps = sum(
            1 for a, b in pairwise(core["points"]) if b["efpd"] != a["efpd"]
        )
        if len(mids) != run_steps:
            raise SystemExit(
                f"{method} {split}: the log has {len(mids)} midpoints for {run_steps} efpd steps"
            )
        point = sum(p["solves"] for p in core["points"])
        solves.append(
            {
                "method": method,
                "split": split,
                "points": len(core["points"]),
                "point_solves": point,
                "midpoint_solves": sum(mids),
                "total": point + sum(mids),
            }
        )
    axial = [
        {
            "method": method,
            "split": split,
            "points": [
                {"efpd": p["efpd"], "days": p["days"], "offset": axial_offset(core, p)}
                for p in core["points"]
            ],
        }
        for (method, split), core in sorted(runs.items())
    ]
    return {
        "specific_power": runs[(methods[0], 1)]["specific_power"],
        "splits": splits,
        "methods": methods,
        "steps": steps,
        "midpoint": midpoint,
        "solves": solves,
        "axial": axial,
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"分割 {summary['splits']}",
    ]
    for method, points in summary["steps"].items():
        lines.append(
            f"[{method}] 燃焼度 (MWd/kg) | 状態 | ホウ素 (ppm) | ホウ素の差 | 燃焼度の差の最大 | 出力の割合の相対差の最大 | 比（ホウ素、燃焼度、出力）"
        )
        for p in points:
            diffs = "; ".join(
                f"{d['splits'][0]}−{d['splits'][1]}: {d['boron']:+.4f} ppm, {d['burnup']:.3e}, {d['power']:.3e}"
                for d in p["differences"]
            )
            ratios = "; ".join(
                ", ".join(
                    "–" if r[q] is None else f"{r[q]:.2f}"
                    for q in ("boron", "burnup", "power")
                )
                for r in p["ratios"]
            )
            lines.append(
                f"{p['burnup']:.4f} | {p['state']} | {p['boron'][-1]:.4f} | {diffs} | {ratios}"
            )
    if summary["midpoint"] is not None:
        for split, points in summary["midpoint"].items():
            boron = max(abs(p["boron"]) for p in points)
            burnup = max(p["burnup"] for p in points)
            power = max(p["power"] for p in points)
            lines.append(
                f"start − search（分割 {split}）: ホウ素の差の絶対値の最大 {boron:.4f} ppm、燃焼度の差の最大 {burnup:.3e} MWd/kg、"
                f"出力の割合の相対差の最大 {power:.3e}"
            )
    for s in summary["solves"]:
        lines.append(
            f"{s['method']} 分割 {s['split']}: 燃焼点 {s['points']}、固有値計算 点 {s['point_solves']} + 中点 {s['midpoint_solves']} = {s['total']}"
        )
    return lines


def parse(arg: str) -> tuple[int, Path, Path]:
    split, _, files = arg.partition("=")
    history, _, log = files.partition(",")
    if not (split.isdigit() and int(split) > 0 and history and log):
        raise SystemExit(f"give <split>=<history>,<log>, got {arg}")
    return int(split), Path(history), Path(log)


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 3 and args[0] == "--json":
        runs = {}
        logs = {}
        files = []
        for split, history, log in map(parse, args[2:]):
            core = msgpack.unpackb(history.read_bytes())
            method = method_name(core)
            if (method, split) in runs:
                raise SystemExit(
                    f"{history}: a second {method} history of split {split}"
                )
            runs[(method, split)] = core
            logs[(method, split)] = log.read_text()
            files += [history, log]
        summary = {
            "format": "xlc-compare-core-steps",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(runs, logs),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
