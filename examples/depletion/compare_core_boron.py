"""集合体 1 体の core の燃焼で、点ごとに探した臨界ホウ素（ADR-0191）を確かめる。

使い方: python3 compare_core_boron.py --json <要約の出力先> <core の履歴（msgpack）> <一様な状態の探索のケース（TOML）> <その標準出力>
        python3 compare_core_boron.py --report <要約の JSON>

- `--json`: core の履歴が集合体 1 体・軸方向のノード 1 つで、状態がホウ素を探す（`search`）ことを確かめ、燃焼点ごとに、点の k と
  目標の差、求めたホウ素濃度と固有値計算の回数を書く。ホウ素濃度が燃焼点の順に下がるか（`boron_falls`）と、|k − 目標| がすべて
  探索の許容の中か（`within_tolerance`）も書く。一様な状態の探索のケースは、点の燃焼度を状態の `burnup` に与え、履歴と同じ
  探索の設定と状態（バンクを含む）で ADR-0145 の探索をするもので、点のノードの燃焼度（解いた値。点の `burnup` は質量の重みの
  平均で、ビットで同じとは限らない）とビットで同じ燃焼度の状態の、標準出力のホウ素濃度（4 桁）、回数、k（8 桁）を並べる。
  点のホウ素濃度を同じ桁に丸めた差（`boron_rounded_difference`）と、回数と k も一致するか（`uniform_agrees`）を書く。
  集合体 1 体では燃焼度が一様なので、点の探索は一様な状態の探索と同じ計算になる。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from itertools import pairwise
from pathlib import Path

import msgpack
import tomllib

# core の状態の出力（crates/core/src/main.rs）の、探索した状態のホウ素濃度と k の行
BORON_LINE = re.compile(r"boron (\S+) ppm \(searched in (\d+) solves\)")
K_LINE = re.compile(r"k-eff (\S+)")


def uniform_searches(case: dict, out: str) -> dict:
    """一様な状態の探索の、状態の名前 → {burnup, boron, solves, k}（標準出力の `state <名前>` に続く行から）"""
    found: dict[str, dict] = {}
    name = None
    for line in out.splitlines():
        if line.startswith("state "):
            name = line.removeprefix("state ")
            found[name] = {}
        elif name is not None and (m := BORON_LINE.fullmatch(line)):
            found[name]["boron"] = float(m[1])
            found[name]["solves"] = int(m[2])
        elif name is not None and (m := K_LINE.fullmatch(line)):
            found[name]["k"] = float(m[1])
    states = case.get("state", {})
    if set(found) != set(states):
        raise SystemExit(
            f"the output has states {sorted(found)}, the case {sorted(states)}"
        )
    for n, s in found.items():
        if set(s) != {"boron", "solves", "k"}:
            raise SystemExit(f"state {n}: the output has no searched boron or k")
        if "burnup" not in states[n]:
            raise SystemExit(f"state {n} gives no burnup")
        s["burnup"] = states[n]["burnup"]
    return found


def summarize(core: dict, case_text: str, out: str) -> dict:
    """履歴、一様な状態の探索のケースと標準出力を比べた要約（モジュールの説明の `--json`）。"""
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"unexpected format {core.get('format')}")
    if len(core["assemblies"]) != 1 or len(core["axial_mesh"]) != 2:
        raise SystemExit("the core must be one assembly of one axial node")
    state = core["state"]
    search = state.get("search")
    if search is None or "boron" in state:
        raise SystemExit("the history state must search the boron")
    case = tomllib.loads(case_text)
    for name, s in case.get("state", {}).items():
        if "search" not in s:
            raise SystemExit(f"state {name} gives no search")
        given = {
            "boron": s["search"]["boron"],
            "k": s["search"]["k"],
            "tolerance": s["search"]["tolerance"],
            "max_solves": s["search"]["max-solves"],
        }
        if (
            given != search
            or any(s[k] != state[k] for k in ("b10", "moderator", "fuel"))
            or s.get("banks", {}) != state["banks"]
        ):
            raise SystemExit(f"state {name} does not search as the history does")
    uniform = uniform_searches(case, out)
    target = search["k"]
    points = []
    for p in core["points"]:
        # 解いたのはノードの燃焼度（点の `burnup` は質量の重みの平均で、ビットで同じとは限らない）
        (burnup,) = p["burnups"][0]
        matched = [n for n, s in uniform.items() if s["burnup"] == burnup]
        if len(matched) != 1:
            raise SystemExit(f"no one uniform state has the burnup {burnup!r}")
        (name,) = matched
        u = uniform[name]
        points.append(
            {
                "efpd": p["efpd"],
                "burnup": burnup,
                "boron": p["boron"],
                "solves": p["solves"],
                "k": p["k"],
                "k_minus_target": p["k"] - target,
                "uniform_state": name,
                "uniform_boron": u["boron"],
                "uniform_solves": u["solves"],
                "uniform_k": u["k"],
                # 一様な状態の出力の桁（ホウ素濃度 4 桁、k 8 桁）に丸めて比べる
                "boron_rounded_difference": round(p["boron"], 4) - u["boron"],
                "k_agrees": round(p["k"], 8) == u["k"],
            }
        )
    largest = max(abs(p["k_minus_target"]) for p in points)
    return {
        "specific_power": core["specific_power"],
        "search": search,
        "points": points,
        "max_abs_k_minus_target": largest,
        "within_tolerance": largest <= search["tolerance"],
        "boron_falls": all(a["boron"] > b["boron"] for a, b in pairwise(points)),
        "uniform_agrees": all(
            p["boron_rounded_difference"] == 0.0
            and p["solves"] == p["uniform_solves"]
            and p["k_agrees"]
            for p in points
        ),
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    points = summary["points"]
    search = summary["search"]
    solves = [p["solves"] for p in points]
    return [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"探索: 範囲 {search['boron']} ppm、目標の k {search['k']}、許容 {search['tolerance']:g}、回数の上限 {search['max_solves']}",
        (
            f"燃焼点 {len(points)}、ホウ素 {points[0]['boron']:.4f}〜{points[-1]['boron']:.4f} ppm、"
            f"点の固有値計算の回数 {min(solves)}〜{max(solves)}（計 {sum(solves)}）"
        ),
        f"|k − 目標| の最大 {summary['max_abs_k_minus_target']:.2e}（許容の中: {summary['within_tolerance']}）",
        f"ホウ素が燃焼点の順に単調に下がる: {summary['boron_falls']}",
        f"一様な状態の探索と一致（4 桁のホウ素、回数、8 桁の k）: {summary['uniform_agrees']}",
        "燃焼度 (MWd/kg) | EFPD | ホウ素 (ppm) | 回数 | k | k − 目標 | 一様のホウ素 (ppm) | 一様の回数 | 丸めた差 (ppm)",
        *(
            f"{p['burnup']:g} | {p['efpd']:.4f} | {p['boron']:.6f} | {p['solves']} | {p['k']:.10f} | "
            f"{p['k_minus_target']:+.2e} | {p['uniform_boron']:.4f} | {p['uniform_solves']} | "
            f"{p['boron_rounded_difference']:+.4f}"
            for p in points
        ),
    ]


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 5 and args[0] == "--json":
        files = [Path(a) for a in args[2:]]
        core = msgpack.unpackb(files[0].read_bytes())
        summary = {
            "format": "xlc-compare-core-boron",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(core, files[1].read_text(), files[2].read_text()),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
