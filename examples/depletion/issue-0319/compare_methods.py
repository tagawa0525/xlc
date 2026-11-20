"""非分離共鳴の負のバックグラウンドの読み方ごとに、同じ評価の遮蔽した群定数を比べる（#0319、ADR-0274）。

使い方: python3 compare_methods.py --json <出力> <範囲の下端 eV> <範囲の上端 eV> <名前>=<xs のファイル>...
        python3 compare_methods.py --report <出力の JSON>

xs のファイルは 1 核種の `xs nuclide` の出力（同じ核種、同じ温度、同じ群と σ0）。最初のファイルを基準にし、ほかのファイルごとに、
温度、σ0、反応（捕獲、核分裂、全断面積）ごとの、非分離共鳴の範囲 [下端, 上端] にかかる群の中で最も大きい相対差
（値/基準 − 1）と、その群、基準の値、値を書く。基準の値が負なら相対差の符号はそのまま（値/基準 − 1）。
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

REACTIONS = ("capture", "fission", "total")


def load(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def differences(reference: dict, other: dict, low: float, high: float) -> list[dict]:
    """温度ごとの、σ0 と反応ごとの最も大きい相対差（群の番号は 1 始まりで、エネルギーの高い順）"""
    bounds = reference["group_bounds"]
    if other["group_bounds"] != bounds or other["sigma0"] != reference["sigma0"]:
        raise SystemExit("the outputs have other groups or sigma zeros")
    a, b = reference["nuclides"][0], other["nuclides"][0]
    temperatures = [t["temperature"] for t in a["temperatures"]]
    if [t["temperature"] for t in b["temperatures"]] != temperatures:
        raise SystemExit("the outputs have other temperatures")
    groups = [
        g for g in range(len(bounds) - 1) if bounds[g + 1] < high and bounds[g] > low
    ]
    out = []
    for ta, tb in zip(a["temperatures"], b["temperatures"]):
        rows = []
        for s0, sa, sb in zip(reference["sigma0"], ta["shielded"], tb["shielded"]):
            row = {"sigma0": s0}
            for reaction in REACTIONS:
                worst = {"relative": 0.0, "group": None, "reference": 0.0, "value": 0.0}
                for g in groups:
                    x, y = sa[reaction][g], sb[reaction][g]
                    if x == 0.0:
                        if y != 0.0:
                            raise SystemExit(
                                f"{reaction} is 0 in group {g + 1} of the reference only"
                            )
                        continue
                    r = y / x - 1.0
                    if worst["group"] is None or abs(r) > abs(worst["relative"]):
                        worst = {
                            "relative": r,
                            "group": g + 1,
                            "reference": x,
                            "value": y,
                        }
                row[reaction] = worst
            rows.append(row)
        out.append({"temperature": ta["temperature"], "sigma0": rows})
    return out


def report_text(summary: dict) -> list[str]:
    lines = [
        f"{summary['nuclide']}: 基準 {summary['reference']}、範囲 {summary['range'][0]:g}〜{summary['range'][1]:g} eV"
    ]
    for name, temperatures in summary["methods"].items():
        for t in temperatures:
            lines.append(f"  {name} / 基準、{t['temperature']:g} K")
            for row in t["sigma0"]:
                cells = []
                for reaction in REACTIONS:
                    d = row[reaction]
                    if d["group"] is None:
                        cells.append(f"{reaction} -")
                        continue
                    cells.append(
                        f"{reaction} {100 * d['relative']:+.3f} %（群 {d['group']}、{d['reference']:.5g} → {d['value']:.5g} b）"
                    )
                lines.append(f"    σ0 {row['sigma0']:g} b: " + "、".join(cells))
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report_text(json.loads(Path(args[1]).read_text()))))
        return
    if len(args) < 6 or args[0] != "--json":
        raise SystemExit(__doc__)
    out, low, high = Path(args[1]), float(args[2]), float(args[3])
    named = [a.split("=", 1) for a in args[4:]]
    if any(len(n) != 2 for n in named):
        raise SystemExit(__doc__)
    data = {name: load(Path(path)) for name, path in named}
    (ref_name, _), *others = named
    reference = data[ref_name]
    nuclide = reference["nuclides"][0]["name"]
    if any(data[n]["nuclides"][0]["name"] != nuclide for n, _ in others):
        raise SystemExit("the outputs are of other nuclides")
    summary = {
        "format": "xlc-unresolved-background-methods",
        "nuclide": nuclide,
        "reference": ref_name,
        "range": [low, high],
        "inputs": [
            {"name": n, "path": p, "md5": hashlib.md5(Path(p).read_bytes()).hexdigest()}
            for n, p in named
        ],
        "methods": {n: differences(reference, data[n], low, high) for n, _ in others},
    }
    out.write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")
    print("\n".join(report_text(summary)))


if __name__ == "__main__":
    main()
