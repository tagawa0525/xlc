"""4 つの版（ENDF/B-VII.0、VII.1、VIII.0、VIII.1）の核データで解いた、同じピンセルの燃焼の履歴（lattice の `[depletion]` の
出力）を、基準の版と比べる（本番の版を決める比較。ADR-0264 の「確認」、#0303）。

使い方: python3 compare_depletion_versions.py --json <要約の出力先> --base <版> <版>=<履歴（msgpack）>...
        python3 compare_depletion_versions.py --report <要約の JSON>

- `--json`: 4 つの版の履歴が、燃焼度の点、燃焼の領域（物質と体積）、燃焼度 0 の組成（0 でない数密度を核種の名前で）で同じこと
  を確かめてから、版ごとに燃焼点ごとの k∞ と、基準の版との差を書く。差は k の差 (k − k_基準) × 1e5 pcm、反応度の差
  (1/k_基準 − 1/k) × 1e5 pcm、主な核種の燃料のピン全体（領域の体積の重みの平均）の数密度の相対差 N / N_基準 − 1（基準が 0 の点は
  比べない）。核種は名前で引く（版ごとのチェーンで核種の並びが違う）。Problem 9 のサイクル 1 の燃焼度に近い点（PROBLEM9）と最後の
  点の値と、それまでと全体の反応度の差の絶対値の最大も書く。`inputs` は読んだファイルのパス（ステージの cwd から）と md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）

出力はどの版もチェーンの核分裂の Q 値で規格化する（lattice の既定）ので、版ごとの Q 値の差は、同じ燃焼度の点の燃焼の時間の差として
比べる値に入る（#0293 の項目 3）。
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

VERSIONS = ["VII.0", "VII.1", "VIII.0", "VIII.1"]
NUCLIDES = ["U-235", "U-238", "Pu-239", "Pu-240", "Pu-241", "Xe-135", "Sm-149"]
# Problem 9 のサイクル 1 の燃焼度 16.939 GWd/MT（docs/plans/m4-depletion.md）に近い、入力の燃焼度の点
PROBLEM9 = 17.0
# 燃焼度、体積、燃焼度 0 の数密度の一致の許容（どれも同じ入力の値）
TOLERANCE = 1e-9


def close(a: float, b: float) -> bool:
    return abs(a - b) <= TOLERANCE * max(abs(a), abs(b), 1.0)


def labelled(args: list[str]) -> dict[str, Path]:
    """`<版>=<パス>` の並びを版 → パスにする"""
    paths = {}
    for arg in args:
        label, sep, path = arg.partition("=")
        if not sep or not label or not path:
            raise SystemExit(f"expected <version>=<history>, got {arg}")
        if label in paths:
            raise SystemExit(f"version {label} given twice")
        paths[label] = Path(path)
    return paths


def initial(history: dict) -> list[dict[str, float]]:
    """燃焼度 0 の点の、領域ごとの 0 でない数密度（名前 → 数密度）"""
    names = history["nuclides"]
    return [
        {n: d for n, d in zip(names, region) if d != 0.0}
        for region in history["points"][0]["densities"]
    ]


def check_same(histories: dict[str, dict], base: str) -> list[float]:
    """4 つの版の履歴が同じ点、領域、初めの組成を持つことを確かめ、燃焼度の点を返す"""
    if sorted(histories) != sorted(VERSIONS):
        raise SystemExit(f"expected the versions {VERSIONS}, got {sorted(histories)}")
    if base not in VERSIONS:
        raise SystemExit(f"base {base} is not one of {VERSIONS}")
    for v, h in histories.items():
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit(f"{v}: not an xlc depletion history")
    ref = histories[base]
    burnup = [p["burnup"] for p in ref["points"]]
    regions = [(r["material"], r["volume"]) for r in ref["regions"]]
    composition = initial(ref)
    for v, h in histories.items():
        other = [p["burnup"] for p in h["points"]]
        if len(other) != len(burnup) or not all(map(close, other, burnup)):
            raise SystemExit(f"{v}: burnup points {other}, {base} has {burnup}")
        mine = [(r["material"], r["volume"]) for r in h["regions"]]
        if len(mine) != len(regions) or any(
            m != n or not close(a, b) for (m, a), (n, b) in zip(mine, regions)
        ):
            raise SystemExit(f"{v}: regions {mine}, {base} has {regions}")
        for r, (a, b) in enumerate(zip(initial(h), composition)):
            if a.keys() != b.keys() or not all(close(a[n], b[n]) for n in b):
                raise SystemExit(
                    f"{v}: region {r} starts from another composition than {base}"
                )
    if not any(close(b, PROBLEM9) for b in burnup):
        raise SystemExit(f"no burnup point at {PROBLEM9} MWd/kg: {burnup}")
    return burnup


def pin_densities(history: dict, name: str) -> list[float]:
    """燃焼点ごとの、燃料のピン全体（領域の体積の重みの平均）の数密度"""
    try:
        i = history["nuclides"].index(name)
    except ValueError:
        raise SystemExit(f"{name} is not in the chain of a history") from None
    volumes = [r["volume"] for r in history["regions"]]
    total = sum(volumes)
    return [
        sum(region[i] * v for region, v in zip(p["densities"], volumes)) / total
        for p in history["points"]
    ]


def largest(burnup: list[float], values: list[float], upto: float) -> dict:
    """燃焼度 `upto` までの点のうち、絶対値の最大の値とその燃焼度"""
    b, x = max(
        ((b, x) for b, x in zip(burnup, values) if b <= upto + TOLERANCE),
        key=lambda bx: abs(bx[1]),
    )
    return {"upto": upto, "burnup": b, "pcm": x}


def summarize(histories: dict[str, dict], base: str) -> dict:
    """4 つの版の履歴を基準の版と比べた要約（モジュールの説明の `--json`）"""
    burnup = check_same(histories, base)
    ref = histories[base]
    k_ref = [p["k_inf"] for p in ref["points"]]
    ref_density = {n: pin_densities(ref, n) for n in NUCLIDES}
    versions = {}
    for v in VERSIONS:
        k = [p["k_inf"] for p in histories[v]["points"]]
        entry = {"k_inf": k}
        if v != base:
            reactivity = [(1 / a - 1 / b) * 1e5 for a, b in zip(k_ref, k)]
            densities = {}
            for n in NUCLIDES:
                ours = pin_densities(histories[v], n)
                densities[n] = [
                    None if r == 0.0 else x / r - 1.0
                    for x, r in zip(ours, ref_density[n])
                ]
            entry.update(
                {
                    "k_pcm": [(b - a) * 1e5 for a, b in zip(k_ref, k)],
                    "reactivity_pcm": reactivity,
                    "densities": densities,
                    "largest_reactivity_pcm": {
                        "problem9": largest(burnup, reactivity, PROBLEM9),
                        "all": largest(burnup, reactivity, burnup[-1]),
                    },
                }
            )
        versions[v] = entry
    marks = []
    for target in (PROBLEM9, burnup[-1]):
        i = next(i for i, b in enumerate(burnup) if close(b, target))
        marks.append(
            {
                "burnup": burnup[i],
                "index": i,
                "reactivity_pcm": {
                    v: e["reactivity_pcm"][i]
                    for v, e in versions.items()
                    if "reactivity_pcm" in e
                },
            }
        )
    return {
        "base": base,
        "nuclides": NUCLIDES,
        "burnup": burnup,
        "versions": versions,
        "marks": marks,
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the stage")
    base = summary["base"]
    others = [v for v in summary["versions"] if v != base]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"基準 {base}。反応度の差は (1/k_{base} − 1/k) × 1e5 pcm",
    ]
    for m in summary["marks"]:
        i = m["index"]
        ks = ", ".join(
            f"{v} {e['k_inf'][i]:.5f}" for v, e in summary["versions"].items()
        )
        diffs = ", ".join(
            f"{v} − {base} {m['reactivity_pcm'][v]:+.0f} pcm" for v in others
        )
        lines.append(f"{m['burnup']:g} MWd/kg: k∞ {ks}; 反応度の差 {diffs}")
        dens = "; ".join(
            f"{v} − {base} "
            + ", ".join(
                f"{n} "
                + (
                    "–"
                    if summary["versions"][v]["densities"][n][i] is None
                    else f"{summary['versions'][v]['densities'][n][i] * 100:+.2f}%"
                )
                for n in summary["nuclides"]
            )
            for v in others
        )
        lines.append(f"{m['burnup']:g} MWd/kg の数密度の相対差: {dens}")
    for v in others:
        big = summary["versions"][v]["largest_reactivity_pcm"]
        lines.append(
            f"{v} − {base} の反応度の差の絶対値の最大: {big['problem9']['upto']:g} MWd/kg までで "
            f"{big['problem9']['pcm']:+.0f} pcm（{big['problem9']['burnup']:g} MWd/kg）、全体で "
            f"{big['all']['pcm']:+.0f} pcm（{big['all']['burnup']:g} MWd/kg）"
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 4 and args[0] == "--json" and args[2] == "--base":
        paths = labelled(args[4:])
        histories = {v: msgpack.unpackb(p.read_bytes()) for v, p in paths.items()}
        summary = {
            "format": "xlc-compare-depletion-versions",
            # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
            "inputs": [
                {"path": str(p), "md5": hashlib.md5(p.read_bytes()).hexdigest()}
                for p in paths.values()
            ],
            **summarize(histories, args[3]),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
