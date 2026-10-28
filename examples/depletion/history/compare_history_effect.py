"""2A の履歴の効き（条件を変えた燃焼と、2A の組成からの瞬時の分岐の差）を、燃焼度の点ごとに出す。

    python3 compare_history_effect.py --json history-effect.json   # このディレクトリの出力から JSON を書く
    python3 compare_history_effect.py --report history-effect.json # 記録の表と文

履歴（`xlc-depletion-history`）の k∞ と数密度、分岐のテーブル（`xlc-table`。燃焼度の軸だけが 2 点以上）の k∞ を読む。反応度の
差は Δρ = (1/k_b − 1/k_a) × 1e5 pcm（a から b への差）。

条件を変えた燃焼（履歴の分岐）ごとに:

- 履歴の効き = ρ(履歴の分岐の k∞) − ρ(2A の組成から、その条件への瞬時の分岐)。同じ燃焼度・同じ条件で、組成の違いだけの差
- 追わない核種の分 = ρ(履歴の分岐) − ρ(2A の組成のうち、ミクロ燃焼で追う元素だけを履歴の分岐の値にした組成の分岐)
- 追う核種の分 = 履歴の効き − 追わない核種の分
- 逆向き = ρ(履歴の分岐の組成から 2A の条件への瞬時の分岐) − ρ(2A の k∞)。2A の条件で、組成だけを履歴の分岐のものにした
  差（条件と組成の交差の項がなければ履歴の効きと等しい）。追わない核種の分は ρ(履歴の分岐の組成から 2A の条件への分岐)
  − ρ(上の置き換えた組成から 2A の条件への分岐)で、追う核種の分は残り
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import msgpack

TRACKED = ("U", "Np", "Pu", "Am", "Xe", "Sm", "Gd")
NUCLIDES = ("U-235", "Pu-239", "Pu-240", "Pu-241", "Xe-135", "Sm-149")

# 条件を変えた燃焼ごとの、履歴と分岐のテーブルのファイル
CASES = {
    "density": {
        "history": "lattice-2A-density-74a7eaa0.msgpack",
        "branch": "table-density-from-2A.msgpack",
        "hybrid": "table-density-hybrid.msgpack",
        "reverse": "table-nominal-from-density.msgpack",
        "reverse_hybrid": "table-density-hybrid-nominal.msgpack",
    },
    "rods": {
        "history": "lattice-2A-rods-db416182.msgpack",
        "branch": "table-rods-from-2A.msgpack",
        "hybrid": "table-rods-hybrid.msgpack",
        "reverse": "table-rods-out-from-rods.msgpack",
        "reverse_hybrid": "table-rods-hybrid-out.msgpack",
    },
}
NOMINAL = "lattice-2A-eb9ca398.msgpack"


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def rho(k_from, k_to):
    """k_from から k_to への反応度の差 (pcm)"""
    return (1 / k_from - 1 / k_to) * 1e5


def table_k(table):
    """燃焼度の点ごとの k∞。燃焼度のほかの軸はどれも 1 点であること"""
    if table["format"] != "xlc-table":
        raise SystemExit("not an xlc-table")
    burnup, *rest = table["axes"]
    if burnup["name"] != "burnup" or any(len(a["points"]) != 1 for a in rest):
        raise SystemExit("the table must vary only the burnup axis")
    if len(table["states"]) != len(burnup["points"]):
        raise SystemExit("the states do not match the burnup points")
    return dict(zip(burnup["points"], (s["k_inf"] for s in table["states"])))


def history_points(history):
    """燃焼度ごとの (k∞, 燃焼の領域の体積で重みを付けた核種の数密度)"""
    if history["format"] != "xlc-depletion-history":
        raise SystemExit("not an xlc-depletion-history")
    volumes = [r["volume"] for r in history["regions"]]
    total = sum(volumes)
    index = {n: i for i, n in enumerate(history["nuclides"])}
    points = {}
    for point in history["points"]:
        mean = {
            n: sum(v * d[index[n]] for v, d in zip(volumes, point["densities"])) / total
            for n in NUCLIDES
        }
        points[point["burnup"]] = (point["k_inf"], mean)
    return points


def compare(nominal, cases):
    """nominal は 2A の履歴、cases は条件ごとの {history, branch, hybrid, reverse, reverse_hybrid}（履歴と 4 つのテーブル）"""
    base = history_points(nominal)
    result = {}
    for name, case in cases.items():
        varied = history_points(case["history"])
        branch, hybrid, reverse, reverse_hybrid = (
            table_k(case[key])
            for key in ("branch", "hybrid", "reverse", "reverse_hybrid")
        )
        burnups = sorted(branch)
        if any(sorted(t) != burnups for t in (hybrid, reverse, reverse_hybrid)):
            raise SystemExit(f"{name}: the tables solve different burnups")
        rows = []
        for u in burnups:
            k_varied, n_varied = varied[u]
            k_base, n_base = base[u]
            effect = rho(branch[u], k_varied)
            untracked = rho(hybrid[u], k_varied)
            reverse_effect = rho(k_base, reverse[u])
            reverse_untracked = rho(reverse_hybrid[u], reverse[u])
            rows.append(
                {
                    "burnup": u,
                    "k_history": k_varied,
                    "k_nominal_history": k_base,
                    "k_branch": branch[u],
                    "k_hybrid": hybrid[u],
                    "k_reverse": reverse[u],
                    "k_reverse_hybrid": reverse_hybrid[u],
                    "between_histories_pcm": rho(k_base, k_varied),
                    "effect_pcm": effect,
                    "untracked_pcm": untracked,
                    "tracked_pcm": effect - untracked,
                    "reverse_pcm": reverse_effect,
                    "reverse_tracked_pcm": reverse_effect - reverse_untracked,
                    "reverse_untracked_pcm": reverse_untracked,
                    "densities_percent": {
                        n: (n_varied[n] / n_base[n] - 1) * 100
                        if n_base[n] > 0
                        else None
                        for n in NUCLIDES
                    },
                }
            )
        result[name] = rows
    return result


def write_json(directory, output):
    d = Path(directory)
    files = [NOMINAL] + [f for case in CASES.values() for f in case.values()]
    cases = {
        name: {key: read(d / f) for key, f in case.items()}
        for name, case in CASES.items()
    }
    summary = {
        "format": "xlc-compare-history-effect",
        "inputs": {f: md5(d / f) for f in files},
        "tracked_elements": list(TRACKED),
        "cases": compare(read(d / NOMINAL), cases),
    }
    Path(output).write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")


def report(path):
    summary = json.loads(Path(path).read_text())
    lines = []
    for name, rows in summary["cases"].items():
        lines.append(f"## {name}")
        lines.append("")
        lines.append(
            "| 燃焼度 (MWd/kg) | 履歴の間 (pcm) | 履歴の効き (pcm) | 追う核種 (pcm) | 追わない核種 (pcm) "
            "| 逆向き (pcm) | 逆向きの追う核種 (pcm) | 逆向きの追わない核種 (pcm) |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for r in rows:
            lines.append(
                f"| {r['burnup']:g} | {r['between_histories_pcm']:+.1f} | {r['effect_pcm']:+.1f} | "
                f"{r['tracked_pcm']:+.1f} | {r['untracked_pcm']:+.1f} | {r['reverse_pcm']:+.1f} | "
                f"{r['reverse_tracked_pcm']:+.1f} | {r['reverse_untracked_pcm']:+.1f} |"
            )
        lines.append("")
        lines.append("| 燃焼度 (MWd/kg) | " + " | ".join(NUCLIDES) + " |")
        lines.append("| --- |" + " --- |" * len(NUCLIDES))
        for r in rows:
            cells = [
                "–" if v is None else f"{v:+.2f}"
                for v in (r["densities_percent"][n] for n in NUCLIDES)
            ]
            lines.append(f"| {r['burnup']:g} | " + " | ".join(cells) + " |")
        lines.append("")
    print("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--json", help="このディレクトリの出力から書く JSON")
    group.add_argument("--report", help="表と文を出す JSON")
    args = parser.parse_args()
    if args.json:
        write_json(Path(__file__).resolve().parent, args.json)
    else:
        report(args.report)


if __name__ == "__main__":
    sys.exit(main())
