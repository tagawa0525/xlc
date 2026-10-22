"""#0133 の試行: 乱数のシードだけを変えた OpenMC の燃焼（1C、1 リング）の 3 本の、燃焼点ごとの k のばらつき。

使い方: python3 seed_scatter.py <出力の JSON> <msgpack>... （最初が今の参照解 openmc-1C-rings1、残りがシードを変えた実行。
`openmc_depletion.py extract` の出力）

燃焼点ごとに、本数ぶんの k と報告された σ、k の標本標準偏差（本数 - 1 で割る）、標準偏差と報告された σ の平均の比、
最初の実行が残りの実行の平均から外れた量（pcm）とその σ（独立な誤差の合成）を書く。
"""

import json
import math
import statistics
import sys
from pathlib import Path

import msgpack

PCM = 1e5


def load(path: Path) -> dict:
    data = msgpack.unpackb(path.read_bytes())
    if data["format"] != "openmc-depletion":
        raise SystemExit(f"{path}: not an openmc-depletion result")
    return data


def summarize(runs: list[dict]) -> list[dict]:
    burnup = runs[0]["burnup"]
    if any(r["burnup"] != burnup for r in runs):
        raise SystemExit("the runs have different burnup points")
    points = []
    for i, b in enumerate(burnup):
        k = [r["k"][i] for r in runs]
        sigma = [r["k-std"][i] for r in runs]
        point = {
            "burnup": b,
            "k": k,
            "sigma": sigma,
            "mean-sigma": statistics.fmean(sigma),
            "stdev": statistics.stdev(k),
        }
        point["stdev-over-sigma"] = point["stdev"] / point["mean-sigma"]
        rest = k[1:]
        deviation = k[0] - statistics.fmean(rest)
        error = math.sqrt(
            sigma[0] ** 2 + sum(s * s for s in sigma[1:]) / len(rest) ** 2
        )
        point["first-minus-rest-pcm"] = deviation * PCM
        point["first-minus-rest-sigma-pcm"] = error * PCM
        points.append(point)
    return points


def table(points: list[dict]) -> str:
    n = len(points[0]["k"])
    head = (
        "| MWd/kg | "
        + " | ".join(f"k (run {i + 1})" for i in range(n))
        + " | 標準偏差 (pcm) | 報告された σ の平均 (pcm) | 比 | 1 本目 − 残りの平均 (pcm) |"
    )
    rows = [head, "|" + "---|" * (n + 5)]
    for p in points:
        ks = " | ".join(f"{k:.5f} ± {s * PCM:.0f}" for k, s in zip(p["k"], p["sigma"]))
        rows.append(
            f"| {p['burnup']:g} | {ks} | {p['stdev'] * PCM:.0f} | {p['mean-sigma'] * PCM:.0f} | "
            f"{p['stdev-over-sigma']:.2f} | {p['first-minus-rest-pcm']:+.0f} ± {p['first-minus-rest-sigma-pcm']:.0f} |"
        )
    return "\n".join(rows)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    points = summarize([load(Path(p)) for p in sys.argv[2:]])
    Path(sys.argv[1]).write_text(json.dumps(points, indent=1) + "\n")
    print(table(points))
