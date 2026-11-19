"""輸送しない核種を輸送する核種にしたライブラリで燃焼させ直した k∞ の差を、untransported_absorption.py の 1 次の目安と比べる（#0292）。

使い方: python3 compare_added.py --json <出力> <足す前の履歴（msgpack）> <足した後の履歴（msgpack）> <吸収の割合の JSON> <足した核種の一覧>
        python3 compare_added.py --report <出力の JSON>

吸収の割合の JSON は、足す前の履歴と、そのベース計算の組成から基準条件へのブランチ計算のテーブルから untransported_absorption.py が
書いたもの（`untransported_all` に核種ごとの割合と自分の吸収で減らした比を持つ）。足した核種の一覧は空白か改行で区切った名前。

燃焼度の点（吸収の割合の JSON の点）ごとに:

- 直接: ρ(足した後の履歴の k∞) − ρ(足す前の履歴の k∞)。Δρ = (1/k_前 − 1/k_後) × 1e5 pcm
- 目安（崩壊だけ）: −1e5 Σ fⱼ / k∞（足した核種の和。fⱼ は崩壊だけで追った数密度の割合）
- 目安（自分の吸収で減らす）: −1e5 Σ fⱼ rⱼ / k∞（rⱼ は自分の吸収で減らした数密度の比）
"""

import argparse
import hashlib
import json
from pathlib import Path

from untransported_absorption import read


def history_k(history):
    if history.get("format") != "xlc-depletion-history":
        raise SystemExit("the histories must be lattice depletion histories")
    return {p["burnup"]: p["k_inf"] for p in history["points"]}


def compare(before, after, absorption, added):
    k_before, k_after = history_k(before), history_k(after)
    gained = set(after["transportable"]) - set(before["transportable"])
    wrong = sorted(set(added) ^ gained)
    if wrong:
        raise SystemExit(
            f"the added list and the nuclides the second library transports differ in {wrong}"
        )
    if sorted(k_before) != sorted(k_after):
        raise SystemExit("the histories have different burnup points")
    rows = []
    for p in absorption["points"]:
        b = p["burnup"]
        if b not in k_before:
            raise SystemExit(f"burnup {b} is not a point of the histories")
        shares = p["untransported_all"]
        dropped = set(absorption.get("pruned", {}).get("nuclides", []))
        missing = [n for n in added if n not in shares and n not in dropped]
        if missing:
            raise SystemExit(
                f"{missing} are not untransported nuclides of the absorption JSON"
            )
        # JSON から外した核種（最大の目安が床以下。untransported_absorption.py の PRUNE_PCM）は目安が 0
        kept = [n for n in added if n in shares]
        decay = sum(shares[n]["fraction"] for n in kept)
        captured = sum(
            shares[n]["fraction"] * shares[n]["captured_ratio"] for n in kept
        )
        rows.append(
            {
                "burnup": b,
                "k_before": k_before[b],
                "k_after": k_after[b],
                "direct_pcm": (1.0 / k_before[b] - 1.0 / k_after[b]) * 1e5,
                "decay_only_pcm": -1e5 * decay / p["k_inf"],
                "captured_pcm": -1e5 * captured / p["k_inf"],
            }
        )
    return {"added": added, "points": rows}


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(f"足した核種 {len(summary['added'])}")
    lines.append(
        "燃焼度 (MWd/kg) | k∞ 前 | k∞ 後 | 直接 (pcm) | 目安 崩壊だけ (pcm) | 目安 自分の吸収で減らす (pcm)"
    )
    for r in summary["points"]:
        lines.append(
            f"{r['burnup']:g} | {r['k_before']:.6f} | {r['k_after']:.6f} | {r['direct_pcm']:+.1f} | "
            f"{r['decay_only_pcm']:+.1f} | {r['captured_pcm']:+.1f}"
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
    before, after = read(args.files[0]), read(args.files[1])
    absorption = json.loads(Path(args.files[2]).read_text())
    added = Path(args.files[3]).read_text().split()
    summary = {
        "format": "xlc-compare-added-nuclides",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in args.files
        ],
        **compare(before, after, absorption, added),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
