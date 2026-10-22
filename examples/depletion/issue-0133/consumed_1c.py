"""#0133 の試行: 1C の燃焼の燃焼点ごとに、U-235 と U-238 の消えた量（ピン全体）と、安定な核分裂生成物（Ru-101、Mo-95。核分裂の
数の積分）の数密度を xlc/OpenMC − 1（%）で比べる。ピン全体はリングの体積の重みの平均（`compare_depletion.py` と同じ）。

使い方: python3 consumed_1c.py <lattice-1C-rings*.msgpack> <openmc-1C-rings*.msgpack>
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_depletion import openmc_name


def main():
    xlc = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
    omc = msgpack.unpackb(Path(sys.argv[2]).read_bytes(), strict_map_key=False)
    index = {n: i for i, n in enumerate(xlc["nuclides"])}
    volumes = np.array([r["volume"] for r in xlc["regions"]])
    weights = volumes / volumes.sum()
    rings = omc["rings"]
    if len(rings) != len(weights):
        raise SystemExit(f"{len(weights)} xlc regions, {len(rings)} OpenMC rings")
    if [p["burnup"] for p in xlc["points"]] != list(omc["burnup"]):
        raise SystemExit("the burnup points differ")
    # リングの体積は xlc の領域と同じ（`compare_depletion.py` が確かめる）なので、OpenMC も xlc の重みで平均する
    for name in ("U-235", "U-238", "Ru-101", "Mo-95"):
        x = [
            sum(w * p["densities"][r][index[name]] for r, w in enumerate(weights))
            for p in xlc["points"]
        ]
        o = [
            sum(
                w * ring["densities"][openmc_name(name)][k]
                for ring, w in zip(rings, weights)
            )
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
