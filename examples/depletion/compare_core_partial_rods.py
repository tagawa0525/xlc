"""部分挿入のノードのミクロ燃焼（ADR-0217）を、挿入の割合を変えた集合体 1 体の core の燃焼で確かめる。

使い方: python3 compare_core_partial_rods.py --json <要約の出力先> <棒のない core の履歴> <棒の入った core の履歴>
          <割合>=<部分挿入の core の履歴> ...
        python3 compare_core_partial_rods.py --report <要約の JSON>

- `--json`: どの core の履歴も、集合体 1 体のミクロ燃焼で、軸方向のメッシュ、追う核種、出力密度、燃焼点（EFPD と燃焼度）と、制御棒の
  バンクの位置を除く状態が同じことを確かめる（違えばエラー）。割合は 0 から 1 へ増え、0 と 1 を
  含むこと。次を書く。
  - 燃焼点ごとの k: 割合ごとの k、割合 0 の k の棒のない燃焼の k からの反応度の差、割合 1 の k の棒の入った燃焼の k からの差
    （Δρ = (1/k_基準 − 1/k) × 1e5 pcm）、k が割合について単調か
  - 追う核種ごと: 割合 0 と棒のない燃焼、割合 1 と棒の入った燃焼の、数密度の相対差 N / N_基準 − 1 の絶対値の最大とその燃焼度と
    軸方向のノードの番号（下から。基準が 0 の点は比べない）、最後の燃焼点のノードごとの割合 1 の割合 0 に対する相対差（割合 0 の
    値が 0 なら null）、数密度が割合について単調でない燃焼点の燃焼度とノード
  単調は、割合の順に並べた値の隣どうしの差が、どれも 0 以上か、どれも 0 以下であること（等しい値は両方に数える）。
  `inputs` は読んだファイルのパスと md5、`banks` は割合ごとの履歴の制御棒のバンクの位置
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import itertools
import json
import sys
from pathlib import Path

import msgpack


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def monotonic(values: list[float]) -> bool:
    steps = [b - a for a, b in itertools.pairwise(values)]
    return all(s >= 0.0 for s in steps) or all(s <= 0.0 for s in steps)


def node_densities(core: dict) -> list[list[list[float]]]:
    """燃焼点ごとの、軸方向のノードごとの追う核種の数密度（[点][ノード][核種]）"""
    out = []
    for p in core["points"]:
        (column,) = p["densities"]
        out.append(column)
    return out


def check_core(name: str, core: dict, first: dict) -> None:
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"{name}: unexpected core format {core.get('format')}")
    if len(core["assemblies"]) != 1:
        raise SystemExit(
            f"{name}: the core must be one assembly (k-infinity of a reflected assembly)"
        )
    if core["axial_mesh"] != first["axial_mesh"]:
        raise SystemExit(f"{name}: the axial mesh differs")
    if "nuclides" not in core:
        raise SystemExit(f"{name}: the core history has no micro depletion")
    if core["nuclides"] != first["nuclides"]:
        raise SystemExit(f"{name}: the tracked nuclides differ")
    if core["specific_power"] != first["specific_power"]:
        raise SystemExit(f"{name}: the specific power differs")
    state = {k: v for k, v in core["state"].items() if k != "banks"}
    if state != {k: v for k, v in first["state"].items() if k != "banks"}:
        raise SystemExit(f"{name}: the state differs beyond the banks")
    ours = [(p["efpd"], p["burnup"]) for p in core["points"]]
    theirs = [(p["efpd"], p["burnup"]) for p in first["points"]]
    if ours != theirs:
        raise SystemExit(f"{name}: the points differ in EFPD or burnup")


def largest(
    names: list[str],
    burnups: list[float],
    values: list[list[list[float]]],
    refs: list[list[list[float]]],
) -> dict:
    """核種ごとの、相対差 N / N_基準 − 1 の絶対値の最大とその燃焼度とノード（基準が 0 の点は比べない）"""
    out = {}
    for i, name in enumerate(names):
        found = [
            (abs(v[i] / r[i] - 1.0), v[i] / r[i] - 1.0, b, z)
            for b, column, reference in zip(burnups, values, refs)
            for z, (v, r) in enumerate(zip(column, reference))
            if r[i] != 0.0
        ]
        _, value, burnup, node = max(found, default=(None, None, None, None))
        out[name] = {"relative": value, "burnup": burnup, "node": node}
    return out


def summarize(out: dict, rodded: dict, partial: list[tuple[float, dict]]) -> dict:
    """棒のない燃焼 `out`、棒の入った燃焼 `rodded`、(割合, 燃焼) の列 `partial` を比べた要約（モジュールの説明の `--json`）"""
    fractions = [f for f, _ in partial]
    if (
        fractions != sorted(set(fractions))
        or fractions[0] != 0.0
        or fractions[-1] != 1.0
    ):
        raise SystemExit(
            f"the fractions {fractions} must increase from 0 to 1 without repeats"
        )
    check_core("out", out, out)
    check_core("rodded", rodded, out)
    for f, core in partial:
        check_core(f"fraction {f}", core, out)
    names = out["nuclides"]
    burnups = [p["burnup"] for p in out["points"]]
    dens = [node_densities(core) for _, core in partial]
    zero, one = dens[0], dens[-1]
    points = []
    for j, p in enumerate(out["points"]):
        ks = [core["points"][j]["k"] for _, core in partial]
        points.append(
            {
                "efpd": p["efpd"],
                "burnup": p["burnup"],
                "out_k": p["k"],
                "rodded_k": rodded["points"][j]["k"],
                "k": ks,
                "zero_pcm": rho_pcm(p["k"], ks[0]),
                "one_pcm": rho_pcm(rodded["points"][j]["k"], ks[-1]),
                "k_monotonic": monotonic(ks),
            }
        )
    zero_limit = largest(names, burnups, zero, node_densities(out))
    one_limit = largest(names, burnups, one, node_densities(rodded))
    nuclides = {}
    nodes = range(len(out["axial_mesh"]) - 1)
    for i, name in enumerate(names):
        nuclides[name] = {
            "zero_limit": zero_limit[name],
            "one_limit": one_limit[name],
            "last_change": [
                None if zero[-1][z][i] == 0.0 else one[-1][z][i] / zero[-1][z][i] - 1.0
                for z in nodes
            ],
            "not_monotonic": [
                {"burnup": b, "node": z}
                for j, b in enumerate(burnups)
                for z in nodes
                if not monotonic([d[j][z][i] for d in dens])
            ],
        }
    return {
        "specific_power": out["specific_power"],
        "fractions": fractions,
        "banks": [core["state"]["banks"] for _, core in partial],
        "axial_mesh": out["axial_mesh"],
        "nuclides": names,
        "points": points,
        "densities": nuclides,
    }


def percent(x: float | None) -> str:
    return "–" if x is None else f"{100.0 * x:+.4f}"


def relative(x: float | None) -> str:
    return "–" if x is None else f"{x:+.3e}"


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    points = summary["points"]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        (
            f"出力密度 {summary['specific_power']:g} W/g、燃焼点 {len(points)}、割合 {summary['fractions']}"
            f"（バンク {summary['banks']}）、追う核種 {len(summary['nuclides'])}"
        ),
    ]
    for key, label in [("zero_pcm", "割合 0 − 棒なし"), ("one_pcm", "割合 1 − 棒あり")]:
        values = [p[key] for p in points]
        lines.append(f"{label}: {min(values):+.6f}〜{max(values):+.6f} pcm")
    bad = [p["burnup"] for p in points if not p["k_monotonic"]]
    lines.append(f"k が割合について単調でない燃焼点: {bad or 'なし'}")
    lines.append(
        "核種 | 割合 0 − 棒なし（相対） | 燃焼度, ノード | 割合 1 − 棒あり（相対） | 燃焼度, ノード | 最後の点の 1 / 0 − 1、ノードごと (%) "
        "| 単調でない (燃焼度, ノード)"
    )
    for name, d in summary["densities"].items():
        z, o = d["zero_limit"], d["one_limit"]
        lines.append(
            f"{name} | {relative(z['relative'])} | {z['burnup']}, {z['node']} | "
            f"{relative(o['relative'])} | {o['burnup']}, {o['node']} | "
            + " ".join(percent(c) for c in d["last_change"])
            + " | "
            + (
                " ".join(f"({m['burnup']:g}, {m['node']})" for m in d["not_monotonic"])
                or "なし"
            )
        )
    lines.append(
        "燃焼度 (MWd/kg) | EFPD | k 棒なし | k 棒あり | k 割合ごと | 0 − 棒なし | 1 − 棒あり (pcm)"
    )
    for p in points:
        lines.append(
            f"{p['burnup']:g} | {p['efpd']:.4f} | {p['out_k']:.8f} | {p['rodded_k']:.8f} | "
            + " ".join(f"{k:.8f}" for k in p["k"])
            + f" | {p['zero_pcm']:+.6f} | {p['one_pcm']:+.6f}"
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 6 and args[0] == "--json":
        pairs = []
        for a in args[4:]:
            fraction, sep, path = a.partition("=")
            if not sep:
                raise SystemExit(f"give <fraction>=<core history>, got {a}")
            pairs.append((float(fraction), Path(path)))
        files = [Path(args[2]), Path(args[3])] + [p for _, p in pairs]
        read = {f: msgpack.unpackb(f.read_bytes()) for f in files}
        summary = {
            "format": "xlc-compare-core-partial-rods",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(
                read[files[0]], read[files[1]], [(f, read[p]) for f, p in pairs]
            ),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
