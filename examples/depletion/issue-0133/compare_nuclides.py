"""#0133 の試行: xlc の燃焼の履歴（複数）と OpenMC の結果の、k∞ の差と核種のピン全体の数密度の相対差 (%) を並べる。

使い方: python3 compare_nuclides.py <OpenMC の msgpack> <名前>=<xlc の履歴> ...
"""

import re
import sys
from pathlib import Path

import msgpack

NUCLIDES = [
    "U-235",
    "Pu-239",
    "Pu-241",
    "Xe-135",
    "I-135",
    "Cs-133",
    "Rh-103",
    "Nd-149",
    "Pm-149",
    "Sm-149",
    "Pm-147",
    "Pm-148",
    "Pm-148m",
    "Sm-150",
    "Sm-151",
    "Sm-152",
    "Eu-153",
    "Nd-143",
    "Nd-145",
    "Tc-99",
    "Ru-101",
]


def openmc_name(name):
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m?)", name).groups()
    return f"{symbol}{mass}" + ("_m1" if isomer else "")


def main():
    omc = msgpack.unpackb(Path(sys.argv[1]).read_bytes())
    rings = omc["rings"]
    total = sum(r["volume"] for r in rings)
    runs = []
    for arg in sys.argv[2:]:
        label, path = arg.split("=", 1)
        runs.append((label, msgpack.unpackb(Path(path).read_bytes())))
    burnup = omc["burnup"]
    print("k∞ の差 (pcm)、± は OpenMC の標準偏差")
    print("燃焼度 | " + " | ".join(label for label, _ in runs))
    for i, b in enumerate(burnup):
        diffs = [(x["points"][i]["k_inf"] - omc["k"][i]) * 1e5 for _, x in runs]
        print(
            f"{b:g} | "
            + " | ".join(f"{d:+.0f}" for d in diffs)
            + f" | ±{omc['k-std'][i] * 1e5:.0f}"
        )
    for name in NUCLIDES:
        theirs = [
            sum(r["densities"][openmc_name(name)][i] * r["volume"] for r in rings)
            / total
            for i in range(len(burnup))
        ]
        print(f"{name}: 燃焼度 | " + " | ".join(label for label, _ in runs))
        for i, b in enumerate(burnup):
            if theirs[i] == 0.0:
                continue
            cells = []
            for _, x in runs:
                k = x["nuclides"].index(name)
                regions = x["regions"]
                ours = sum(
                    x["points"][i]["densities"][r][k] * reg["volume"]
                    for r, reg in enumerate(regions)
                ) / sum(reg["volume"] for reg in regions)
                cells.append(f"{(ours / theirs[i] - 1) * 100:+.2f}")
            print(f"  {b:g} | " + " | ".join(cells))


main()
