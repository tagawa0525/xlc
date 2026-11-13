"""#0286 の、変種（inputs.py）ごとの時間と、既定（変種 default）からの差を比べる。

ケースと変種ごとに、少数群定数のファイル（`fg-<ケース>-<変種>.msgpack`）とログ（`<ケース>-<変種>.log`。ログの水準は
`info,xlc_lattice::ufg=debug,xlc_lattice::shielding=debug`）から次を読み、既定との差を書く:

- k∞ の反応度の差 1/k(既定) − 1/k(変種)（pcm）
- 2 群の定数（全断面積、吸収、νΣf、散乱行列）の相対の差の最大
- ピン出力の形状関数の相対の差の最大（集合体だけ。テーブルの `pins`）
- 燃料のリングごとの U-238 の、UFG の範囲の中性子束で縮約した 1 群の捕獲断面積の相対の差の最大（ログの `UFG U-238 ring capture`）
- 時間: 自己遮蔽の span の秒（`shielding: close` の和）と、解いたピンの UFG の段の和（スレッド秒。`UFG slowing down` と
  `UFG transfers counted`）を解いたピンの数で割った 1 ピンの値

使い方: python3 compare.py --json <要約の出力先> [<ケース>...]（このディレクトリで。ケースを省くと、ログのあるものすべて）
"""

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


INPUTS = load("inputs", HERE / "inputs.py")
CONSTANTS = ("total", "absorption", "nu_fission", "scatter")
# `UFG slowing down` の段（秒）。slowing は free_gas と sab を含む
STAGES = ("grid", "table", "evaluate", "collision", "slowing", "free_gas", "sab")


def flat(value) -> list[float]:
    """入れ子の数の並びを 1 列に"""
    if isinstance(value, (list, tuple)):
        return [x for v in value for x in flat(v)]
    if value is None:
        return []
    return [float(value)]


def relative(a: list[float], b: list[float]) -> float:
    """`b` に対する `a` の相対の差の最大（`b` が 0 の要素は、`a` も 0 なら差 0、違えば無限大）"""
    if len(a) != len(b):
        sys.exit(f"cannot compare {len(a)} values with {len(b)}")
    worst = 0.0
    for x, y in zip(a, b):
        if y == 0.0:
            worst = max(worst, 0.0 if x == 0.0 else float("inf"))
        else:
            worst = max(worst, abs(x - y) / abs(y))
    return worst


def form(state: dict) -> list[float]:
    """状態のピン出力の形状関数（テーブルの `pins.form`。無ければ空）"""
    return flat(state["pins"]["form"]) if state["pins"] else []


def states(path: Path) -> list[dict]:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)["states"]


FIELD = re.compile(r"(\w+)=(\S+)")
# 自己遮蔽の span（子の span の `shielding:mix:` などは除く）の閉じた行
SPAN = re.compile(
    r" INFO shielding: xlc_lattice::shielding: close time\.busy=([\d.]+)(ms|s|µs)"
)


def fields(line: str) -> dict:
    return dict(FIELD.findall(line))


def log(path: Path) -> dict:
    """ログの時間と、[ピン][リング] の U-238 の捕獲"""
    shielding = 0.0
    stages = dict.fromkeys(STAGES, 0.0)
    transfers = 0.0
    pins = 0
    captures = []
    total = None
    for line in path.read_text().splitlines():
        busy = SPAN.search(line)
        if busy:
            shielding += (
                float(busy.group(1)) * {"s": 1.0, "ms": 1e-3, "µs": 1e-6}[busy.group(2)]
            )
        elif "UFG slowing down" in line:
            f = fields(line)
            pins += 1
            for s in STAGES:
                stages[s] += float(f[s])
        elif "UFG transfers counted" in line:
            transfers += float(fields(line)["seconds"])
        elif "UFG U-238 ring capture" in line:
            f = fields(line)
            captures.append((int(f["cell"]), int(f["ring"]), float(f["capture"])))
        elif "total_seconds=" in line and "lattice state=" in line:
            total = float(fields(line)["total_seconds"])
    return {
        "shielding_seconds": shielding,
        "total_seconds": total,
        "solved": pins,
        "per_pin": {s: v / pins for s, v in stages.items()} if pins else None,
        "transfers_per_pin": transfers / pins if pins else None,
        "captures": sorted(captures),
    }


def difference(case: str, variant: str, base: dict) -> dict:
    fg = states(HERE / f"fg-{case}-{variant}.msgpack")
    timing = log(HERE / f"{case}-{variant}.log")
    row = {"variant": variant, **{k: v for k, v in timing.items() if k != "captures"}}
    if base is None:
        return row | {"k_inf": [s["k_inf"] for s in fg]}
    row["k_inf"] = [s["k_inf"] for s in fg]
    if len(fg) != len(base["fg"]):
        sys.exit(
            f"{case}-{variant}: {len(fg)} states against {len(base['fg'])} of the default"
        )
    row["drho_pcm"] = [
        (1.0 / b["k_inf"] - 1.0 / s["k_inf"]) * 1e5 for s, b in zip(fg, base["fg"])
    ]
    row["constants"] = max(
        relative(flat(s[c]), flat(b[c]))
        for s, b in zip(fg, base["fg"])
        for c in CONSTANTS
    )
    pins = [(form(s), form(b)) for s, b in zip(fg, base["fg"])]
    row["pins"] = (
        max(relative(a, b) for a, b in pins) if any(b for _, b in pins) else None
    )
    ours, theirs = timing["captures"], base["captures"]
    if [c[:2] for c in ours] != [c[:2] for c in theirs]:
        sys.exit(
            f"{case}-{variant}: the ring captures are not of the same pins and rings"
        )
    row["ring_capture"] = (
        relative([c[2] for c in ours], [c[2] for c in theirs]) if ours else None
    )
    return row


def compare(case: str, variants: list[str]) -> list[dict]:
    base_fg = states(HERE / f"fg-{case}-default.msgpack")
    base = {"fg": base_fg, "captures": log(HERE / f"{case}-default.log")["captures"]}
    rows = [difference(case, "default", None)]
    rows += [difference(case, v, base) for v in variants if v != "default"]
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("cases", nargs="*")
    args = parser.parse_args()
    cases = args.cases or [
        c for c in INPUTS.cases() if (HERE / f"{c}-default.log").exists()
    ]
    summary = {}
    for case in cases:
        variants = [v for v in INPUTS.VARIANTS if (HERE / f"{case}-{v}.log").exists()]
        summary[case] = compare(case, variants)
    args.json.write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
