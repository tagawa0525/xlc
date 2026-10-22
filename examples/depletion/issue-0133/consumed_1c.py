"""#0133 の試行: 1C の燃焼の燃焼点ごとに、U-235 と U-238 の消えた量（ピン全体）と、安定な核分裂生成物（Ru-101、Mo-95。核分裂の
数の積分）の数密度を xlc/OpenMC − 1（%）で比べる。ピン全体はリングの体積の重みの平均で、2 つの履歴の対応（形式、輸送する核種、
燃焼点、リングの名前と体積）は `compare_depletion.py` の `summarize` で確かめる。

使い方: python3 consumed_1c.py <lattice-1C-rings*.msgpack> <openmc-1C-rings*.msgpack>
"""

import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_depletion import openmc_name, summarize


def main():
    xlc = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
    omc = msgpack.unpackb(Path(sys.argv[2]).read_bytes(), strict_map_key=False)
    summarize(xlc, omc)
    index = {n: i for i, n in enumerate(xlc["nuclides"])}
    regions, rings = xlc["regions"], omc["rings"]
    total = sum(r["volume"] for r in regions)

    def pin(values):
        return sum(v * r["volume"] for v, r in zip(values, regions)) / total

    for name in ("U-235", "U-238", "Ru-101", "Mo-95"):
        x = [
            pin([p["densities"][r][index[name]] for r in range(len(regions))])
            for p in xlc["points"]
        ]
        o = [
            pin([ring["densities"][openmc_name(name)][k] for ring in rings])
            for k in range(len(omc["burnup"]))
        ]
        if name.startswith("U-"):
            values = [((x[0] - a) / (o[0] - b) - 1) * 100 for a, b in zip(x[1:], o[1:])]
            label = "consumed"
        else:
            values = [(a / b - 1) * 100 for a, b in zip(x[1:], o[1:])]
            label = "density"
        print(
            f"{name} {label} xlc/OpenMC - 1 (%):", " ".join(f"{v:+.3f}" for v in values)
        )


main()
