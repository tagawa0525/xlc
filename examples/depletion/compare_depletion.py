"""xlc の格子計算の燃焼の履歴（lattice の `[depletion]` の出力、ADR-0173）と、OpenMC の燃焼（openmc_depletion.py extract）を比べる。

使い方: python3 compare_depletion.py --json <要約の出力先> <xlc の履歴（msgpack）> <OpenMC の結果（msgpack）>
        python3 compare_depletion.py --report <要約の JSON>

- `--json`: 燃料で輸送する核種の集合が同じこと、燃焼度の点とリングの体積が同じことを確かめてから、燃焼点ごとの k∞ の差
  (k_xlc − k_OpenMC) × 1e5 pcm（OpenMC の標準偏差つき）とその 0 の点からの変化、主要な核種の数密度の相対差
  (N_xlc − N_OpenMC) / N_OpenMC を、燃料のピン全体（リングの体積の重みの平均）と、最後の点のリングごとに書く。xlc の燃焼の
  領域と OpenMC のリングは、内側から順に対応させる。`inputs` は読んだファイルのパス（段階の cwd から）と md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack

NUCLIDES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"]
# 燃焼度とリングの体積の一致の許容（燃焼度は入力の同じ値、体積は xlc の格子の円の面積と OpenMC の面積）
BURNUP_TOLERANCE = 1e-9
VOLUME_TOLERANCE = 1e-6


def openmc_name(name: str) -> str:
    """xlc の核種名（`Am-242m`）を OpenMC の名前（`Am242_m1`）にする。"""
    m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m(\d*))?", name)
    if m is None:
        raise SystemExit(f"unexpected nuclide name {name}")
    symbol, mass, isomer, level = m.groups()
    return f"{symbol}{mass}" + (f"_m{level or 1}" if isomer else "")


def relative(ours: float, theirs: float) -> float | None:
    """相対差。OpenMC が 0 なら比べない（None）"""
    return None if theirs == 0.0 else ours / theirs - 1.0


def summarize(xlc: dict, omc: dict) -> dict:
    """2 つの履歴を比べた要約（モジュールの説明の `--json`）。"""
    if xlc["format"] != "xlc-depletion-history" or omc["format"] != "openmc-depletion":
        raise SystemExit("unexpected formats")
    ours = {openmc_name(n) for n in xlc["transportable"]}
    theirs = set(omc["transportable"])
    if ours != theirs:
        raise SystemExit(
            f"transported nuclides differ: xlc only {sorted(ours - theirs)}, "
            f"OpenMC only {sorted(theirs - ours)}"
        )
    burnup = [p["burnup"] for p in xlc["points"]]
    if len(burnup) != len(omc["burnup"]) or any(
        abs(a - b) > BURNUP_TOLERANCE * max(b, 1.0)
        for a, b in zip(burnup, omc["burnup"])
    ):
        raise SystemExit(f"burnup points differ: {burnup} vs {omc['burnup']}")
    regions, rings = xlc["regions"], omc["rings"]
    if len(regions) != len(rings):
        raise SystemExit(f"{len(regions)} xlc regions for {len(rings)} OpenMC rings")
    for r, ring in zip(regions, rings):
        if abs(r["volume"] - ring["volume"]) > VOLUME_TOLERANCE * ring["volume"]:
            raise SystemExit(
                f"{r['material']}: volume {r['volume']} vs {ring['volume']}"
            )
    index = {name: i for i, name in enumerate(xlc["nuclides"])}
    total = sum(r["volume"] for r in regions)

    def pin(values) -> float:
        return sum(v * r["volume"] for v, r in zip(values, regions)) / total

    first = (xlc["points"][0]["k_inf"] - omc["k"][0]) * 1e5
    k = []
    for i, point in enumerate(xlc["points"]):
        diff = (point["k_inf"] - omc["k"][i]) * 1e5
        k.append(
            {
                "xlc": point["k_inf"],
                "openmc": omc["k"][i],
                "openmc_std_pcm": omc["k-std"][i] * 1e5,
                "diff_pcm": diff,
                "change_pcm": diff - first,
            }
        )
    density = {}
    for name in NUCLIDES:
        series = []
        for i, point in enumerate(xlc["points"]):
            ours_pin = pin(
                [point["densities"][r][index[name]] for r in range(len(regions))]
            )
            theirs_pin = pin(
                [ring["densities"][openmc_name(name)][i] for ring in rings]
            )
            series.append(relative(ours_pin, theirs_pin))
        density[name] = series
    last = len(burnup) - 1
    by_ring = [
        {
            name: relative(
                xlc["points"][last]["densities"][r][index[name]],
                ring["densities"][openmc_name(name)][last],
            )
            for name in NUCLIDES
        }
        for r, ring in enumerate(rings)
    ]
    return {
        "burnup": burnup,
        "transported": len(ours),
        "k": k,
        "pin": density,
        "rings": by_ring,
    }


def percent(x: float | None) -> str:
    return "–" if x is None else f"{x * 100:+.2f}"


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the stage")
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"輸送する核種: 両方で同じ {summary['transported']} 核種",
        "燃焼度 (MWd/kg) | k xlc | k OpenMC | 差 (pcm) | 0 からの変化 (pcm)",
    ]
    for b, k in zip(summary["burnup"], summary["k"]):
        lines.append(
            f"{b:g} | {k['xlc']:.6f} | {k['openmc']:.6f} ± {k['openmc_std_pcm'] / 1e5:.6f} | "
            f"{k['diff_pcm']:+.0f} ± {k['openmc_std_pcm']:.0f} | {k['change_pcm']:+.0f}"
        )
    lines.append("ピン全体の数密度の相対差 (%): 燃焼度 | " + " | ".join(NUCLIDES))
    for i, b in enumerate(summary["burnup"]):
        lines.append(
            f"{b:g} | " + " | ".join(percent(summary["pin"][n][i]) for n in NUCLIDES)
        )
    if len(summary["rings"]) > 1:
        lines.append(
            "最後の点のリングごとの相対差 (%): リング | " + " | ".join(NUCLIDES)
        )
        for r, ring in enumerate(summary["rings"]):
            lines.append(f"{r} | " + " | ".join(percent(ring[n]) for n in NUCLIDES))
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 4 and args[0] == "--json":
        files = [Path(args[2]), Path(args[3])]
        xlc, omc = (msgpack.unpackb(f.read_bytes()) for f in files)
        summary = {
            "format": "xlc-compare-depletion",
            # 読んだファイルのパス（段階の cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(xlc, omc),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
