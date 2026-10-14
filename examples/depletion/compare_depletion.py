"""xlc の格子計算の燃焼の履歴（lattice の `[depletion]` の出力、ADR-0173）と、OpenMC の燃焼（openmc_depletion.py extract）を比べる。

使い方: python3 compare_depletion.py <xlc の履歴（msgpack）> <OpenMC の結果（msgpack）>

燃料で輸送する核種の集合が同じことを確かめてから、燃焼点ごとに、k∞ の差 (k_xlc − k_OpenMC) × 1e5 pcm（± は OpenMC の標準偏差）と、その 0 の点からの変化、主要な核種の数密度の
相対差 (N_xlc − N_OpenMC) / N_OpenMC を、燃料のピン全体（リングの体積の重みの平均）とリングごとに書く。xlc の燃焼の領域と
OpenMC のリングは、内側から順に対応させ、体積が合うことを確かめる。
"""

import re
import sys
from pathlib import Path

import msgpack

NUCLIDES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"]
# 燃焼度と体積の一致の許容（燃焼度は入力の同じ値、体積は xlc の格子の円の面積と OpenMC の面積）
TOLERANCE = 1e-9


def openmc_name(name: str) -> str:
    """xlc の核種名（`Am-242m`）を OpenMC の名前（`Am242_m1`）にする。"""
    m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m(\d*))?", name)
    if m is None:
        raise SystemExit(f"unexpected nuclide name {name}")
    symbol, mass, isomer, level = m.groups()
    return f"{symbol}{mass}" + (f"_m{level or 1}" if isomer else "")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    xlc = msgpack.unpackb(Path(sys.argv[1]).read_bytes())
    omc = msgpack.unpackb(Path(sys.argv[2]).read_bytes())
    if xlc["format"] != "xlc-depletion-history" or omc["format"] != "openmc-depletion":
        raise SystemExit("unexpected formats")
    burnup = [p["burnup"] for p in xlc["points"]]
    if len(burnup) != len(omc["burnup"]) or any(
        abs(a - b) > TOLERANCE * max(b, 1.0) for a, b in zip(burnup, omc["burnup"])
    ):
        raise SystemExit(f"burnup points differ: {burnup} vs {omc['burnup']}")
    regions, rings = xlc["regions"], omc["rings"]
    if len(regions) != len(rings):
        raise SystemExit(f"{len(regions)} xlc regions for {len(rings)} OpenMC rings")
    for r, ring in zip(regions, rings):
        if abs(r["volume"] - ring["volume"]) > 1e-6 * ring["volume"]:
            raise SystemExit(
                f"{r['material']}: volume {r['volume']} vs {ring['volume']}"
            )
    index = {name: i for i, name in enumerate(xlc["nuclides"])}
    # 燃料で輸送する核種の集合が同じこと（ADR-0173）
    ours = {openmc_name(n) for n in xlc["transportable"]}
    theirs = set(omc["transportable"])
    if ours != theirs:
        raise SystemExit(
            f"transported nuclides differ: xlc only {sorted(ours - theirs)}, OpenMC only {sorted(theirs - ours)}"
        )
    print(f"輸送する核種: 両方で同じ {len(ours)} 核種\n")

    print("| 燃焼度 (MWd/kg) | k xlc | k OpenMC | 差 (pcm) | 0 からの変化 (pcm) |")
    print("| --- | --- | --- | --- | --- |")
    first = None
    for i, point in enumerate(xlc["points"]):
        k, sigma = omc["k"][i], omc["k-std"][i]
        diff = (point["k_inf"] - k) * 1e5
        first = diff if first is None else first
        print(
            f"| {burnup[i]:g} | {point['k_inf']:.6f} | {k:.6f} ± {sigma:.6f} | {diff:+.0f} ± {sigma * 1e5:.0f} | "
            f"{diff - first:+.0f} |"
        )

    def pin(values: list[float]) -> float:
        total = sum(r["volume"] for r in regions)
        return sum(v * r["volume"] for v, r in zip(values, regions)) / total

    print("\n燃料のピン全体の数密度の相対差（%）\n")
    print("| 燃焼度 (MWd/kg) | " + " | ".join(NUCLIDES) + " |")
    print("| --- |" + " --- |" * len(NUCLIDES))
    for i, point in enumerate(xlc["points"]):
        cells = []
        for name in NUCLIDES:
            ours = pin(
                [point["densities"][r][index[name]] for r in range(len(regions))]
            )
            theirs = pin([ring["densities"][openmc_name(name)][i] for ring in rings])
            cells.append(
                "–" if theirs == 0.0 else f"{(ours - theirs) / theirs * 100:+.2f}"
            )
        print(f"| {burnup[i]:g} | " + " | ".join(cells) + " |")

    if len(regions) > 1:
        print("\n最後の点のリングごとの数密度の相対差（%）\n")
        print("| リング | " + " | ".join(NUCLIDES) + " |")
        print("| --- |" + " --- |" * len(NUCLIDES))
        last = len(burnup) - 1
        for r, ring in enumerate(rings):
            cells = []
            for name in NUCLIDES:
                ours = xlc["points"][last]["densities"][r][index[name]]
                theirs = ring["densities"][openmc_name(name)][last]
                cells.append(
                    "–" if theirs == 0.0 else f"{(ours - theirs) / theirs * 100:+.2f}"
                )
            print(f"| {r} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main()
