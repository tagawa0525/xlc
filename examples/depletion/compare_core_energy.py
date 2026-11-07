"""ノードの出力の割合を核分裂率 Σf φ V の割合で求めた core の燃焼（ADR-0190 の「## 変更」の前の決定）と、κΣf φ V の割合で求めた燃焼
（今の決定。#0154）を比べ、燃焼点ごとの k、ホウ素、反応度の差、ノードの出力の割合と燃焼度の相対差、集合体の出力の割合と燃焼度を出す。

使い方: python3 compare_core_energy.py --json <要約の出力先> <名前>=<核分裂率の履歴>,<κΣf の履歴> ...
        python3 compare_core_energy.py --report <要約の JSON>

- 履歴は core の燃焼の履歴（msgpack）。2 つの燃焼は、同じ燃焼点（EFPD、日数、状態）、定格の出力密度、状態、軸方向のメッシュ、集合体、
  重金属の質量を持つこと（違えば止める）
- `--json`: 名前ごとに `assemblies`（集合体の ID）と、燃焼点ごとの `points` を書く
  - `efpd`、`days`、`state`
  - `k`、`boron`: [核分裂率, κΣf] の値。`reactivity` は反応度の差 (1/k_核分裂率 − 1/k_κΣf) × 1e5 pcm
  - `power`、`burnup`: ノードの出力の割合と燃焼度の、核分裂率の燃焼を基準にした相対差の絶対値の最大（核分裂率の燃焼で 0 のノードは、κΣf の
    燃焼も 0 であること）
  - `assembly_powers`、`assembly_burnups`: `rate`（核分裂率）と `energy`（κΣf）の、集合体ごとの出力の割合（ノードの和）と燃焼度（ノードの
    重金属の質量の重み。質量のない集合体は 0）
  `inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

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


def assembly_burnups(
    burnups: list[list[float]], masses: list[list[float]]
) -> list[float]:
    """集合体ごとの、ノードの重金属の質量の重みの燃焼度（質量のない集合体は 0）"""
    out = []
    for b, m in zip(burnups, masses):
        total = sum(m)
        out.append(sum(x * y for x, y in zip(b, m)) / total if total > 0.0 else 0.0)
    return out


def compare(rate: dict, energy: dict) -> dict:
    for key in SHARED:
        if rate.get(key) != energy.get(key):
            raise SystemExit(f"the histories differ in {key}")
    key = [(p["efpd"], p["days"], p["state"]) for p in rate["points"]]
    if key != [(p["efpd"], p["days"], p["state"]) for p in energy["points"]]:
        raise SystemExit("the histories have different points")
    masses = rate["heavy_metal"]
    points = []
    for p, q in zip(rate["points"], energy["points"]):
        points.append(
            {
                "efpd": p["efpd"],
                "days": p["days"],
                "state": p["state"],
                "k": [p["k"], q["k"]],
                "boron": [p["boron"], q["boron"]],
                "reactivity": (1.0 / p["k"] - 1.0 / q["k"]) * 1e5,
                "power": largest_relative(flat(q["powers"]), flat(p["powers"])),
                "burnup": largest_relative(flat(q["burnups"]), flat(p["burnups"])),
                "assembly_powers": {
                    "rate": [sum(nodes) for nodes in p["powers"]],
                    "energy": [sum(nodes) for nodes in q["powers"]],
                },
                "assembly_burnups": {
                    "rate": assembly_burnups(p["burnups"], masses),
                    "energy": assembly_burnups(q["burnups"], masses),
                },
            }
        )
    return {"assemblies": rate["assemblies"], "points": points}


def assembly_power_difference(point: dict) -> float:
    """集合体の出力の割合の、核分裂率の燃焼を基準にした相対差の絶対値の最大"""
    shares = point["assembly_powers"]
    return largest_relative(shares["energy"], shares["rate"])


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        points = c["points"]
        boron = max(abs(p["boron"][1] - p["boron"][0]) for p in points)
        reactivity = max(abs(p["reactivity"]) for p in points)
        power = max(p["power"] for p in points)
        assembly = max(assembly_power_difference(p) for p in points)
        burnup = max(p["burnup"] for p in points)
        lines.append(
            f"{name}: ホウ素の差の絶対値の最大 {boron:.4f} ppm、反応度の差の絶対値の最大 {reactivity:.2f} pcm、ノードの出力の割合の"
            f"相対差の絶対値の最大 {power:.2e}、集合体の出力の割合の相対差の絶対値の最大 {assembly:.2e}、ノードの燃焼度の相対差の"
            f"絶対値の最大 {burnup:.2e}"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) != 2 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<fission rate history>,<kappa history>, got {arg}"
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
            rate, energy = paths
            cases[name] = compare(
                msgpack.unpackb(rate.read_bytes()), msgpack.unpackb(energy.read_bytes())
            )
            files += paths
        summary = {
            "format": "xlc-compare-core-energy",
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
