"""足す核種の、実際の組成で自己遮蔽した 1 群の吸収断面積を無限希釈と比べ、反応度の目安の差に直す（#0292）。

使い方: python3 compare_shielding.py --json <出力> <無限希釈のテーブル> <自己遮蔽のテーブル> <吸収の割合の JSON> <足した核種の一覧>
        python3 compare_shielding.py --report <出力の JSON>

2 つのテーブルは、同じ履歴の組成から同じ状態へのブランチ計算で、冷却計算の 1 群の微視的断面積（ADR-0238）を持つ。無限希釈のテーブルは
今のライブラリ（足す核種は放射化の核種で、無限希釈の群の断面積）、自己遮蔽のテーブルは足す核種を輸送する核種にしたライブラリのコピー
（自己遮蔽した断面積）で解いたもの。足す核種が輸送に入るので、自己遮蔽のテーブルの中性子束のスペクトルも少し変わり、比はその分を含む。

燃焼度の点（吸収の割合の JSON の点）ごとに、核種 j の 1 群の吸収断面積の比 r = σₐ(自己遮蔽) / σₐ(無限希釈) と、自分の吸収で減らした
1 次の目安の差 1e5 fⱼ cⱼ (1 − r) / k∞ (pcm。fⱼ と cⱼ は JSON の割合と数密度の比。正なら無限希釈の目安が吸収を大きく見る) を出し、
燃焼度を通じた差の絶対値の最大の大きい順に並べる。
"""

import argparse
import hashlib
import json
from pathlib import Path

import msgpack
from untransported_absorption import absorption, nominal_states


def read(path):
    return msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)


def sigmas(table, names):
    """燃焼度の値ごとの、核種の 1 群の吸収断面積"""
    burnups, _ = nominal_states(table)
    cooling = table["cooling"]
    missing = [n for n in names if n not in cooling["nuclides"]]
    if missing:
        raise SystemExit(f"the table has no cooling cross sections of {missing}")
    index = {n: i for i, n in enumerate(cooling["nuclides"])}
    return {
        b: {n: absorption(cooling["points"][k][index[n]]) for n in names}
        for k, b in enumerate(burnups)
    }


def compare(infinite, shielded, shares, names):
    a, b = sigmas(infinite, names), sigmas(shielded, names)
    rows = []
    for name in names:
        points = []
        for p in shares["points"]:
            burnup = p["burnup"]
            share = p["untransported_all"].get(name)
            if share is None or burnup not in a or burnup not in b:
                raise SystemExit(f"{name} at {burnup} MWd/kg is missing")
            ratio = b[burnup][name] / a[burnup][name]
            points.append(
                {
                    "burnup": burnup,
                    "ratio": ratio,
                    "shielding_pcm": 1e5
                    * share["fraction"]
                    * share["captured_ratio"]
                    * (1.0 - ratio)
                    / p["k_inf"],
                }
            )
        largest = max(points, key=lambda x: abs(x["shielding_pcm"]))
        rows.append(
            {
                "nuclide": name,
                "max_pcm": largest["shielding_pcm"],
                "burnup": largest["burnup"],
                "points": points,
            }
        )
    return sorted(rows, key=lambda r: -abs(r["max_pcm"]))


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        "核種 | 自己遮蔽による目安の差の最大 (pcm) | 燃焼度 | その点の σₐ の比"
    )
    for r in summary["nuclides"]:
        at = next(p for p in r["points"] if p["burnup"] == r["burnup"])
        lines.append(
            f"{r['nuclide']} | {r['max_pcm']:+.2f} | {r['burnup']:g} | {at['ratio']:.4f}"
        )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not args.json or len(args.files) != 4:
        raise SystemExit(__doc__)
    infinite, shielded = read(args.files[0]), read(args.files[1])
    shares = json.loads(Path(args.files[2]).read_text())
    names = Path(args.files[3]).read_text().split()
    summary = {
        "format": "xlc-compare-shielding",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in args.files
        ],
        "nuclides": compare(infinite, shielded, shares, names),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
