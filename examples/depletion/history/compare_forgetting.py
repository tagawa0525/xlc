"""追わない核種の分 U を、指数の忘却の履歴の変数の重みで R(1) から内挿したときの外れを比べる（ADR-0215、#0190）。

    python3 compare_forgetting.py --json forgetting.json     # このディレクトリの出力から JSON を書く
    python3 compare_forgetting.py --report forgetting.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

履歴の変数を燃焼度について指数で忘れる平均

    H_τ(B) = ∫₀ᴮ x(B') e^{−(B − B')/τ} dB' / ∫₀ᴮ e^{−(B − B')/τ} dB'

にすると、途中で条件を変えた燃焼（燃焼度 B₀ まで公称、その後を条件 j）の重みは
w_τ(B) = (1 − e^{−(B − B₀)/τ}) / (1 − e^{−B/τ})。τ → ∞ は燃焼度で重みを付けた平均（ADR-0215 の案 V1。w = (B − B₀)/B）、正の τ は最近の
燃焼に、負の τ は古い燃焼に大きな重みを付ける。一生を条件 j の R(1) は τ によらず w = 1。

compare_split.py と同じテーブルと式で、R(1) と switch の U を作り、switch の参照の点（B₀ より後）ごとに、U(1) を w_τ 倍して足したときの
外れ δ = U_switch − w_τ U(1) の分 ρ(参照の定数 − δ) − ρ(参照の定数)（compare_index.py の外れの分と同じ式）を、τ の列ごとに出す。6 点の
外れの絶対値の最大が最小の τ を `best` にする（τ の列の中から選ぶ。∞ は null）。
"""

import argparse
import json
import math
import re
from pathlib import Path

import compare_split as split
from compare_core_history import same_burnup
from compare_index import gap
from compare_residual import md5, read

FILES = {
    key: split.FILES[key]
    for key in ("nominal", "one", "one_hybrid", "switch", "switch_hybrid")
}
# 比べる τ (MWd/kg)。None は ∞（燃焼度で重みを付けた平均）
TAUS = [
    None,
    *(
        sign * t
        for t in (1000, 500, 200, 100, 60, 40, 30, 20, 15, 10)
        for sign in (1, -1)
    ),
]


def weight(burnup, start, tau):
    """燃焼度 `start` まで公称、その後を条件 j にした履歴の、燃焼度 `burnup` の w_τ（モジュールの説明）。`tau` が None なら ∞"""
    if tau is None:
        return (burnup - start) / burnup
    return -math.expm1(-(burnup - start) / tau) / -math.expm1(-burnup / tau)


def compare(tables, start, taus=TAUS):
    """要約（モジュールの説明）。`start` は switch が条件を変えた燃焼度 B₀"""
    nominal = tables["nominal"]
    _, one = split.split(nominal, tables["one"], tables["one_hybrid"], "R(1)")
    _, switch = split.split(
        nominal, tables["switch"], tables["switch_hybrid"], "switch"
    )
    later = [x for x in switch if x["burnup"] > start]
    pairs = []
    for x in later:
        found = [o for o in one if same_burnup(o["burnup"], x["burnup"])]
        if len(found) != 1:
            raise SystemExit(f"R(1) has no single point at {x['burnup']}")
        pairs.append((x, found[0]))
    rows = []
    for tau in taus:
        points = []
        for x, o in pairs:
            w = weight(x["burnup"], start, tau)
            points.append(
                {
                    "burnup": x["burnup"],
                    "weight": w,
                    "gap_pcm": gap(x["branch"], x["untracked"], o["untracked"], w),
                }
            )
        rows.append(
            {
                "tau": tau,
                "points": points,
                "largest_pcm": max((abs(p["gap_pcm"]) for p in points), default=0.0),
            }
        )
    best = min(rows, key=lambda r: r["largest_pcm"])
    return {
        "format": "xlc-compare-history-forgetting",
        "continued_burnup": start,
        "best": best["tau"],
        "rows": rows,
    }


def report(summary):
    """要約の文章（記録の本文の数値はここからコピーする。ADR-0176）"""

    def name(tau):
        return "∞" if tau is None else f"{tau:g}"

    rows = summary["rows"]
    best = next(r for r in rows if r["tau"] == summary["best"])
    gaps = [p["gap_pcm"] for p in best["points"]]
    lines = [
        f"U の外れの絶対値の最大が最小の τ は {name(best['tau'])} MWd/kg（比べた τ: "
        + "、".join(name(r["tau"]) for r in rows)
        + f"）。そのときの外れは {min(gaps):+.3f}〜{max(gaps):+.3f} pcm（最大 {best['largest_pcm']:.3f} pcm）",
        "τ ごとの外れの絶対値の最大 (pcm): "
        + "、".join(f"{name(r['tau'])} {r['largest_pcm']:.3f}" for r in rows),
    ]
    return re.sub(r"(?<![A-Za-z])-", "−", "\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", help="write the summary to this file")
    parser.add_argument("--report", help="print the report of this summary")
    args = parser.parse_args()
    if args.report:
        with open(args.report, encoding="utf-8") as f:
            print(report(json.load(f)))
        return
    here = Path(__file__).resolve().parent
    tables = {k: read(here / v) for k, v in FILES.items()}
    # switch が条件を変えた燃焼度は、参照の最初の点（格子計算の履歴の `continued` の点。compare_switch.py が照らす）
    start = tables["switch"]["axes"][0]["points"][0]
    result = compare(tables, start)
    result["inputs"] = {v: md5(here / v) for v in FILES.values()}
    text = json.dumps(result, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
