"""集合体 1 体の core のミクロ燃焼（ADR-0216）を、格子計算の燃焼の履歴（ADR-0173）とテーブルの均質化した数密度と比べる。

使い方: python3 compare_core_micro.py --json <要約の出力先> <格子計算の履歴（msgpack）> <テーブル（msgpack）>
          <ミクロ燃焼の core の履歴（msgpack）> <ミクロ燃焼しない core の履歴（msgpack）>
        python3 compare_core_micro.py --report <要約の JSON>

- `--json`: 2 つの core の履歴が同じ集合体 1 体・軸方向のノード 1 つ（反射境界で k∞ を解くケース）、同じ状態と燃焼点で、燃焼点の燃焼度が
  格子計算の履歴の燃焼点のどれかと相対 1e-9 で一致することを確かめて、次を書く。
  - 数密度: テーブルの燃焼度の点ごとに、追う核種の core の数密度の、テーブルの均質化した数密度 N̄⁰（公称の履歴）に対する相対差
    N / N̄⁰ − 1（N̄⁰ が 0 の核種は null）と、核種ごとの絶対値の最大とその燃焼度
  - k: 燃焼点ごとに、ミクロ燃焼の k の履歴の k∞ からの反応度の差、ミクロ燃焼しない k の履歴の k∞ からの差（テーブルの燃焼度の内挿の
    誤差。ADR-0192）、ミクロ燃焼の k のミクロ燃焼しない k からの差（線形の補正の効き）。Δρ = (1/k_基準 − 1/k) × 1e5 pcm
  `inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

# core の燃焼度と履歴とテーブルの燃焼点の一致の許容（core は出力密度 × EFPD で、履歴は燃焼度の点から時間を求めて燃焼させる）
BURNUP_TOLERANCE = 1e-9


def same_burnup(a: float, b: float) -> bool:
    return abs(a - b) <= BURNUP_TOLERANCE * max(1.0, abs(b))


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def check_core(core: dict) -> None:
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"unexpected core format {core.get('format')}")
    if len(core["assemblies"]) != 1 or len(core["axial_mesh"]) != 2:
        raise SystemExit(
            "the core must be one assembly of one axial node (k-infinity of a reflected assembly)"
        )


def summarize(lattice: dict, table: dict, micro: dict, plain: dict) -> dict:
    """4 つのファイルを比べた要約（モジュールの説明の `--json`）。"""
    if lattice.get("format") != "xlc-depletion-history":
        raise SystemExit(f"unexpected history format {lattice.get('format')}")
    if table.get("format") != "xlc-table" or "micro" not in table:
        raise SystemExit("the table must be an xlc-table with micro data")
    check_core(micro)
    check_core(plain)
    if "nuclides" not in micro or "nuclides" in plain:
        raise SystemExit(
            "give the micro-depleted core history first and the one without micro depletion second"
        )
    if (
        micro["state"] != plain["state"]
        or micro["specific_power"] != plain["specific_power"]
    ):
        raise SystemExit("the two core histories differ in their state or power")
    names = micro["nuclides"]
    if names != table["micro"]["nuclides"]:
        raise SystemExit("the core tracks other nuclides than the table")
    axes = [a for a in table["axes"] if a["name"] == "burnup"]
    if len(axes) != 1:
        raise SystemExit("the table has no burnup axis")
    table_burnups = axes[0]["points"]
    densities = table["micro"]["densities"]
    if len(micro["points"]) != len(plain["points"]):
        raise SystemExit("the two core histories have different points")
    points = []
    density_points = []
    for p, q in zip(micro["points"], plain["points"]):
        b = p["burnup"]
        if not same_burnup(q["burnup"], b):
            raise SystemExit(f"the core points differ in burnup ({b}, {q['burnup']})")
        matched = [h for h in lattice["points"] if same_burnup(b, h["burnup"])]
        if len(matched) != 1:
            raise SystemExit(f"core burnup {b} is not a point of the history")
        (h,) = matched
        on_table = [k for k, t in enumerate(table_burnups) if same_burnup(b, t)]
        points.append(
            {
                "efpd": p["efpd"],
                "burnup": h["burnup"],
                "on_table": bool(on_table),
                "history_k": h["k_inf"],
                "micro_k": p["k"],
                "plain_k": q["k"],
                "micro_pcm": rho_pcm(h["k_inf"], p["k"]),
                "plain_pcm": rho_pcm(h["k_inf"], q["k"]),
                "correction_pcm": rho_pcm(q["k"], p["k"]),
            }
        )
        if on_table:
            (k,) = on_table
            ((node,),) = p["densities"]
            density_points.append(
                {
                    "burnup": table_burnups[k],
                    "relative": {
                        name: None if ref == 0.0 else n / ref - 1.0
                        for name, n, ref in zip(names, node, densities[k])
                    },
                }
            )
    largest = {}
    for name in names:
        values = [
            (abs(d["relative"][name]), d["relative"][name], d["burnup"])
            for d in density_points
            if d["relative"][name] is not None
        ]
        _, value, burnup = max(values, default=(None, None, None))
        largest[name] = {"relative": value, "burnup": burnup}
    return {
        "specific_power": micro["specific_power"],
        "nuclides": names,
        "table_burnups": table_burnups,
        "points": points,
        "density_points": density_points,
        "density_largest": largest,
    }


def percent(x: float | None) -> str:
    return "–" if x is None else f"{100.0 * x:+.4f}"


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    points = summary["points"]
    on = [p for p in points if p["on_table"]]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"出力密度 {summary['specific_power']:g} W/g、燃焼点 {len(points)}（テーブルの点 {len(on)}）、追う核種 {len(summary['nuclides'])}",
    ]
    for key, label in [
        ("micro_pcm", "ミクロ燃焼 − 履歴"),
        ("plain_pcm", "ミクロ燃焼なし − 履歴"),
        ("correction_pcm", "ミクロ燃焼 − ミクロ燃焼なし"),
    ]:
        between = [p for p in points if not p["on_table"]]
        for name, subset in [
            ("すべての点", points),
            ("テーブルの点", on),
            ("点の間", between),
        ]:
            if not subset:
                continue
            values = [p[key] for p in subset]
            worst = max(subset, key=lambda p: abs(p[key]))
            lines.append(
                f"{label}（{name}）: {min(values):+.3f}〜{max(values):+.3f} pcm（|Δρ| の最大は {worst['burnup']:g} MWd/kg の {worst[key]:+.3f} pcm）"
            )
    lines.append("核種 | 相対差の絶対値の最大 (%) | 燃焼度 (MWd/kg)")
    for name, w in summary["density_largest"].items():
        lines.append(f"{name} | {percent(w['relative'])} | {w['burnup']}")
    lines.append(
        "燃焼度 (MWd/kg) | EFPD | テーブルの点 | k 履歴 | k ミクロ燃焼 | k なし | ミクロ − 履歴 | なし − 履歴 | ミクロ − なし (pcm)"
    )
    for p in points:
        lines.append(
            f"{p['burnup']:g} | {p['efpd']:.4f} | {'○' if p['on_table'] else ''} | {p['history_k']:.8f} | "
            f"{p['micro_k']:.8f} | {p['plain_k']:.8f} | {p['micro_pcm']:+.3f} | {p['plain_pcm']:+.3f} | "
            f"{p['correction_pcm']:+.3f}"
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 6 and args[0] == "--json":
        files = [Path(a) for a in args[2:]]
        lattice, table, micro, plain = (msgpack.unpackb(f.read_bytes()) for f in files)
        summary = {
            "format": "xlc-compare-core-micro",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(lattice, table, micro, plain),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
