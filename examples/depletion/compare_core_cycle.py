"""サイクルのステップ（ADR-0218: 出力の割合、崩壊だけのステップ、ステップごとの状態）で燃焼させた集合体 1 体の core の燃焼を、
燃焼度の解析値と、同じ EFPD の参照の core の燃焼と比べる。

使い方: python3 compare_core_cycle.py --json <要約の出力先> <サイクルの core の履歴（msgpack）> [<名前>=<参照の core の履歴（msgpack）> ...]
        python3 compare_core_cycle.py --report <要約の JSON>

- `--json`: core の履歴がどれも集合体 1 体で、軸方向のメッシュと定格の出力密度が同じで、サイクルの履歴が I-135、Xe-135、Sm-149 を追う
  ことを確かめて、サイクルの燃焼点ごとに次を書く。
  - 燃焼度: 炉心の平均の燃焼度と、解析値（定格の出力密度 × 累積の EFPD / 1000。EFPD は出力の割合を含むので、割合によらない）、
    その相対差（解析値が 0 の点は 0）、軸方向のノードごとの燃焼度
  - 点の EFPD、燃焼点 0 からの日数、その点で終わるステップの出力の割合（崩壊だけのステップは 0、燃焼点 0 は null）、点を解いた状態の名前、k
  - 数密度: I-135、Xe-135、Sm-149 の軸方向のノードごとの数密度
  - 参照: 参照ごとに、同じ EFPD の点の k と反応度の差 Δρ = (1/k_参照 − 1/k) × 1e5 pcm、その点の数密度。崩壊だけのステップの点と、参照に
    同じ EFPD の点がない点は null
  `inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

# 追う数密度（Xe-135 とその先行核、停止の間に Pm-149 の崩壊で増える Sm-149）
NUCLIDES = ["I-135", "Xe-135", "Sm-149"]


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def check_core(name: str, core: dict, cycle: dict) -> None:
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"{name}: unexpected core format {core.get('format')}")
    if len(core["assemblies"]) != 1:
        raise SystemExit(f"{name}: the core must be one assembly")
    if core["axial_mesh"] != cycle["axial_mesh"]:
        raise SystemExit(f"{name}: the axial mesh differs from the cycle depletion")
    if core["specific_power"] != cycle["specific_power"]:
        raise SystemExit(f"{name}: the specific power differs from the cycle depletion")
    missing = [n for n in NUCLIDES if n not in core.get("nuclides", [])]
    if missing:
        raise SystemExit(f"{name}: the core does not track {missing}")


def densities(core: dict, point: dict) -> dict:
    """点の `NUCLIDES` の、軸方向のノードごとの数密度"""
    index = [core["nuclides"].index(n) for n in NUCLIDES]
    (column,) = point["densities"]
    return {n: [node[i] for node in column] for n, i in zip(NUCLIDES, index)}


def summarize(cycle: dict, references: dict[str, dict]) -> dict:
    """サイクルの燃焼と参照の燃焼を比べた要約（モジュールの説明の `--json`）。"""
    check_core("cycle", cycle, cycle)
    for name, core in references.items():
        check_core(name, core, cycle)
    power = cycle["specific_power"]
    points = []
    for p in cycle["points"]:
        analytic = power * p["efpd"] * 1e-3
        paired = {}
        for name, core in references.items():
            # 崩壊だけのステップの点は EFPD が前の点と同じなので、出力のある点（燃焼点 0 を含む）だけを同じ EFPD の点と比べる
            same = [
                q
                for q in core["points"]
                if q["efpd"] == p["efpd"] and q.get("power") != 0.0
            ]
            if p.get("power") == 0.0 or len(same) != 1:
                paired[name] = None
                continue
            (q,) = same
            paired[name] = {
                "k": q["k"],
                "rho_pcm": rho_pcm(q["k"], p["k"]),
                "densities": densities(core, q),
            }
        points.append(
            {
                "efpd": p["efpd"],
                "days": p["days"],
                "power": p.get("power"),
                "state": p["state"],
                "burnup": p["burnup"],
                "analytic_burnup": analytic,
                "burnup_relative": 0.0
                if analytic == 0.0
                else p["burnup"] / analytic - 1.0,
                "node_burnups": p["burnups"][0],
                "k": p["k"],
                "densities": densities(cycle, p),
                "references": paired,
            }
        )
    return {
        "specific_power": power,
        "axial_mesh": cycle["axial_mesh"],
        "nuclides": NUCLIDES,
        "references": list(references),
        "points": points,
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    points = summary["points"]
    worst = max(points, key=lambda p: abs(p["burnup_relative"]))
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        (
            f"出力密度 {summary['specific_power']:g} W/g、燃焼点 {len(points)}、"
            f"燃焼度の解析値との相対差の絶対値の最大 {abs(worst['burnup_relative']):.1e}（{worst['efpd']:g} EFPD）"
        ),
    ]
    for name in summary["references"]:
        paired = [p for p in points if p["references"][name] is not None]
        if not paired:
            lines.append(f"参照 {name}: 同じ EFPD の点なし")
            continue
        values = [p["references"][name]["rho_pcm"] for p in paired]
        lines.append(
            f"参照 {name}: 同じ EFPD の点 {len(paired)}、Δρ {min(values):+.3f}〜{max(values):+.3f} pcm"
        )
    header = "EFPD | 日数 | 出力 | 状態 | 燃焼度 | k | " + " | ".join(
        f"Δρ {name} (pcm)" for name in summary["references"]
    )
    lines.append(header + " | " + " | ".join(summary["nuclides"]))
    for p in points:
        refs = [
            ""
            if p["references"][n] is None
            else f"{p['references'][n]['rho_pcm']:+.3f}"
            for n in summary["references"]
        ]
        nuclides = [
            ", ".join(f"{v:.6e}" for v in p["densities"][n])
            for n in summary["nuclides"]
        ]
        power = "" if p["power"] is None else f"{p['power']:g}"
        lines.append(
            f"{p['efpd']:.4f} | {p['days']:.4f} | {power} | {p['state']} | {p['burnup']:.6f} | {p['k']:.8f} | "
            + " | ".join(refs + nuclides)
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 3 and args[0] == "--json":
        cycle_path = Path(args[2])
        named = []
        for a in args[3:]:
            name, sep, path = a.partition("=")
            if not sep or not name:
                raise SystemExit(f"give a reference as <name>=<core history>, got {a}")
            named.append((name, Path(path)))
        if len({n for n, _ in named}) != len(named):
            raise SystemExit("the reference names must differ")
        files = [cycle_path, *(p for _, p in named)]
        cycle = msgpack.unpackb(cycle_path.read_bytes())
        references = {n: msgpack.unpackb(p.read_bytes()) for n, p in named}
        summary = {
            "format": "xlc-compare-core-cycle",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(cycle, references),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
