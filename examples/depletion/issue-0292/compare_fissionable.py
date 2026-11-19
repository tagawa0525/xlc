"""核分裂の生成が効く核種も輸送する核種にして、残りを無限希釈の吸収体として数えた燃焼（ADR-0261）の k∞ を、65 核種を輸送して残りを
数えた燃焼と比べる（#0292）。

使い方: python3 compare_fissionable.py --json <出力> [--reference <69 核種の希釈>] <前> <65 核種の希釈> <13 核種の希釈>
          <65 核種の一覧> <65 核種の外の核分裂の核種の一覧> <自己遮蔽する核種の一覧> <核分裂の核種の一覧>
        python3 compare_fissionable.py --report <出力の JSON>

履歴（msgpack）は同じセルの燃焼で、ライブラリだけが違う。前の履歴は吸収体を数えず（輸送する核種 T）、ほかは数える:

- 65 核種の希釈: T と 65 核種 A を輸送し、残り R を数える（compare_dilute.py の「足した希釈」）
- 69 核種の希釈（あれば）: T、A と、核分裂の核種のうち A の外の核種 N を輸送し、R − N を数える（比べる相手）
- 13 核種の希釈: T、自己遮蔽する核種 S と核分裂の核種 F を輸送し、ほかをすべて数える（調べる形）

燃焼度の点ごとの反応度の差（pcm。Δρ(X − Y) = (1/k_Y − 1/k_X) × 1e5）:

- form_from_65_pcm: 13 核種の希釈 − 65 核種の希釈
- form_pcm（69 核種の希釈があるとき）: 13 核種の希釈 − 69 核種の希釈。A ∪ N のうち S と F の外の核種を、自己遮蔽して輸送する代わりに
  無限希釈で数えた差
- rest_fission_pcm（同）: 69 核種の希釈 − 65 核種の希釈。N を輸送する核種にした差（65 核種の希釈が落とした N の核分裂の生成など）。
  form_from_65_pcm = form_pcm + rest_fission_pcm

自己遮蔽の時間は、履歴の `seconds.shielding` の、前の履歴に対する比。
"""

import argparse
import hashlib
import json
from pathlib import Path

from untransported_absorption import read

NAMES = ["before", "added_dilute", "reference", "form"]


def check_sets(histories, added, outside, shielded, fissionable):
    """4 つの履歴の輸送する核種と吸収体が、上の組であること。吸収体の数を返す"""
    base = set(histories["before"]["transportable"])
    added, outside, shielded, fissionable = (
        set(added),
        set(outside),
        set(shielded),
        set(fissionable),
    )
    if outside != fissionable - added:
        raise SystemExit(
            f"the fissionable nuclides outside the added ones are {sorted(fissionable - added)}, "
            f"not {sorted(outside)}"
        )
    if not shielded <= added:
        raise SystemExit(
            f"the shielded nuclides {sorted(shielded - added)} are not added"
        )
    rest = set(histories["added_dilute"].get("dilute_absorbers", []))
    if not outside <= rest:
        raise SystemExit(
            f"{sorted(outside - rest)} are not dilute absorbers of the 65-nuclide burnup"
        )
    everything = added | rest
    expected = {
        "before": (base, set()),
        "added_dilute": (base | added, rest),
        "reference": (base | added | outside, rest - outside),
        "form": (base | shielded | fissionable, everything - shielded - fissionable),
    }
    for name, (transportable, dilute) in expected.items():
        if name not in histories:
            continue
        h = histories[name]
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit(f"{name}: not a lattice depletion history")
        got = (set(h["transportable"]), set(h.get("dilute_absorbers", [])))
        if got != (transportable, dilute):
            raise SystemExit(
                f"{name}: transportable differ in {sorted(got[0] ^ transportable)}, dilute absorbers in "
                f"{sorted(got[1] ^ dilute)}"
            )
    return {
        name: len(histories[name].get("dilute_absorbers", []))
        for name in NAMES
        if name in histories
    }


def compare(histories, added, outside, shielded, fissionable):
    counts = check_sets(histories, added, outside, shielded, fissionable)
    names = [n for n in NAMES if n in histories]
    burnups = [p["burnup"] for p in histories["before"]["points"]]
    for name in names:
        if [p["burnup"] for p in histories[name]["points"]] != burnups:
            raise SystemExit(
                f"{name}: other burnup points than the history without absorbers"
            )
    k = {name: [p["k_inf"] for p in histories[name]["points"]] for name in names}

    def rho(x, y, i):
        return (1.0 / k[y][i] - 1.0 / k[x][i]) * 1e5

    points = []
    for i, b in enumerate(burnups):
        point = {
            "burnup": b,
            **{f"k_{name}": k[name][i] for name in names},
            "form_from_65_pcm": rho("form", "added_dilute", i),
        }
        if "reference" in histories:
            point["form_pcm"] = rho("form", "reference", i)
            point["rest_fission_pcm"] = rho("reference", "added_dilute", i)
        points.append(point)
    before = histories["before"]["seconds"]["shielding"]
    seconds = {
        name: {
            "shielding": histories[name]["seconds"]["shielding"],
            "shieldings": histories[name]["shieldings"],
            "shielding_ratio": histories[name]["seconds"]["shielding"] / before,
        }
        for name in names
    }
    return {
        "added": sorted(added),
        "outside": sorted(outside),
        "shielded": sorted(shielded),
        "fissionable": sorted(fissionable),
        "counts": counts,
        "seconds": seconds,
        "points": points,
    }


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        f"自己遮蔽する核種 {' '.join(summary['shielded'])}、核分裂の核種 {' '.join(summary['fissionable'])}、"
        f"65 核種の外の核分裂の核種 {' '.join(summary['outside'])}、吸収体の数 "
        + ", ".join(f"{n} {c}" for n, c in summary["counts"].items())
    )
    lines.append(
        "自己遮蔽の時間 (s、回数、前に対する比): "
        + ", ".join(
            f"{n} {s['shielding']:.1f}、{s['shieldings']}、{s['shielding_ratio']:.3f}"
            for n, s in summary["seconds"].items()
        )
    )
    lines.append(
        "| 燃焼度 (MWd/kg) | 調べる形 − 65 核種 | 調べる形 − 比べる相手 | 比べる相手 − 65 核種 |"
    )
    lines.append("| --- | --- | --- | --- |")
    for r in summary["points"]:
        rest = (
            f"{r['form_pcm']:+.1f} | {r['rest_fission_pcm']:+.1f}"
            if "form_pcm" in r
            else "− | −"
        )
        lines.append(f"| {r['burnup']:g} | {r['form_from_65_pcm']:+.1f} | {rest} |")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument("--reference")
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not args.json or len(args.files) != 7:
        raise SystemExit(__doc__)
    histories = {
        name: read(f)
        for name, f in zip(["before", "added_dilute", "form"], args.files[:3])
    }
    if args.reference:
        histories["reference"] = read(args.reference)
    lists = [Path(f).read_text().split() for f in args.files[3:]]
    summary = {
        "format": "xlc-compare-fissionable-absorbers",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in [*args.files, *([args.reference] if args.reference else [])]
        ],
        **compare(histories, *lists),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
