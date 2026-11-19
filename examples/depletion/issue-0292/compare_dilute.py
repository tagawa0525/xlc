"""輸送しない核種を無限希釈の吸収体として数えた燃焼（ADR-0261）の k∞ を、65 核種を輸送する核種にした燃焼と比べる（#0292）。

使い方: python3 compare_dilute.py --json <出力> <前> <足した> <希釈> <足した希釈> <遮蔽の希釈> <足した核種の一覧> <自己遮蔽する核種の一覧>
        python3 compare_dilute.py --report <出力の JSON>

履歴（msgpack）は同じセルの燃焼で、ライブラリと入力の `dilute-absorbers` だけが違う:

- 前: p2 のライブラリ、数えない（輸送する核種 T）
- 足した: 足した核種 A を輸送する核種にしたライブラリ、数えない（T + A。比べる相手）
- 希釈: p2 のライブラリ、数える（T、吸収体は A と残り R）
- 足した希釈: 足した核種のライブラリ、数える（T + A、吸収体は R）
- 遮蔽の希釈: 自己遮蔽する核種 S（A の一部）だけを輸送する核種にしたライブラリ、数える（T + S、吸収体は A − S と R）

燃焼度の点ごとの反応度の差（pcm。Δρ(X − Y) = (1/k_Y − 1/k_X) × 1e5）:

- dilute_pcm: 希釈 − 足した。足した核種と残りを無限希釈で数えたときの、比べる相手との差
- shielded_dilute_pcm: 遮蔽の希釈 − 足した
- rest_pcm: 足した希釈 − 足した。残り R の吸収
- dilute_added_pcm: 希釈 − 足した希釈。足した核種 A を、自己遮蔽して輸送する代わりに無限希釈の吸収体として数えた差
  （自己遮蔽、散乱、核分裂の生成、(n,2n) などを落とした分）。dilute_pcm = rest_pcm + dilute_added_pcm
- dilute_unshielded_pcm: 遮蔽の希釈 − 足した希釈。A − S を無限希釈の吸収体として数えた差

自己遮蔽の時間は、履歴の `seconds.shielding`（重ねる吸収の組み立てを含む）の、前の履歴に対する比。
"""

import argparse
import hashlib
import json
from pathlib import Path

from untransported_absorption import read

NAMES = ["before", "added", "dilute", "added_dilute", "shielded_dilute"]


def check_sets(histories, added, shielded):
    """5 つの履歴の輸送する核種と吸収体が、上の組であること。吸収体の数を返す"""
    base = set(histories["before"]["transportable"])
    # 残り R は、希釈の吸収体から足した核種を除いたもの
    rest = set(histories["dilute"].get("dilute_absorbers", [])) - set(added)
    if not set(shielded) <= set(added):
        raise SystemExit(
            f"the shielded nuclides {sorted(set(shielded) - set(added))} are not added"
        )
    expected = {
        "before": (base, set()),
        "added": (base | set(added), set()),
        "dilute": (base, set(added) | rest),
        "added_dilute": (base | set(added), rest),
        "shielded_dilute": (base | set(shielded), (set(added) - set(shielded)) | rest),
    }
    if not rest:
        raise SystemExit(
            "the depletion with the p2 library counts no dilute absorbers besides the added ones"
        )
    for name, (transportable, dilute) in expected.items():
        h = histories[name]
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit(f"{name}: not a lattice depletion history")
        got = (set(h["transportable"]), set(h.get("dilute_absorbers", [])))
        if got != (transportable, dilute):
            raise SystemExit(
                f"{name}: transportable differ in {sorted(got[0] ^ transportable)}, dilute absorbers in "
                f"{sorted(got[1] ^ dilute)}"
            )
    return {name: len(histories[name].get("dilute_absorbers", [])) for name in NAMES}


def compare(histories, added, shielded):
    counts = check_sets(histories, added, shielded)
    ks = {
        name: [(p["burnup"], p["k_inf"]) for p in histories[name]["points"]]
        for name in NAMES
    }
    burnups = [b for b, _ in ks["before"]]
    for name in NAMES:
        if [b for b, _ in ks[name]] != burnups:
            raise SystemExit(
                f"{name}: other burnup points than the history without absorbers"
            )
    k = {name: [v for _, v in ks[name]] for name in NAMES}

    def rho(x, y, i):
        return (1.0 / k[y][i] - 1.0 / k[x][i]) * 1e5

    points = [
        {
            "burnup": b,
            **{f"k_{name}": k[name][i] for name in NAMES},
            "dilute_pcm": rho("dilute", "added", i),
            "shielded_dilute_pcm": rho("shielded_dilute", "added", i),
            "rest_pcm": rho("added_dilute", "added", i),
            "dilute_added_pcm": rho("dilute", "added_dilute", i),
            "dilute_unshielded_pcm": rho("shielded_dilute", "added_dilute", i),
        }
        for i, b in enumerate(burnups)
    ]
    before = histories["before"]["seconds"]["shielding"]
    seconds = {
        name: {
            "shielding": histories[name]["seconds"]["shielding"],
            "shieldings": histories[name]["shieldings"],
            "shielding_ratio": histories[name]["seconds"]["shielding"] / before,
        }
        for name in NAMES
    }
    return {
        "added": list(added),
        "shielded": list(shielded),
        "counts": counts,
        "seconds": seconds,
        "points": points,
    }


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        f"足した核種 {len(summary['added'])}、自己遮蔽する核種 {' '.join(summary['shielded'])}、吸収体の数 "
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
        "| 燃焼度 (MWd/kg) | 全核種を数える | 4 核種を遮蔽 | 残りの吸収 | 65 核種を希釈 | 61 核種を希釈 |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for r in summary["points"]:
        lines.append(
            f"| {r['burnup']:g} | {r['dilute_pcm']:+.1f} | {r['shielded_dilute_pcm']:+.1f} | {r['rest_pcm']:+.1f} | "
            f"{r['dilute_added_pcm']:+.1f} | {r['dilute_unshielded_pcm']:+.1f} |"
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
    if not args.json or len(args.files) != 7:
        raise SystemExit(__doc__)
    histories = {name: read(f) for name, f in zip(NAMES, args.files[:5])}
    added = Path(args.files[5]).read_text().split()
    shielded = Path(args.files[6]).read_text().split()
    summary = {
        "format": "xlc-compare-dilute-absorbers",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in args.files
        ],
        **compare(histories, added, shielded),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
