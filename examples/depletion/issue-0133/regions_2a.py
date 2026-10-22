"""#0133 の試行: 2A の対称操作の組ごとに、U-235 の消えた量と U-238 の消えた量の xlc/OpenMC − 1（%）を燃焼点ごとに出す。

使い方: python3 regions_2a.py <lattice-2A.msgpack> <openmc-2A.msgpack> <燃焼点の番号 ...>
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_depletion_assembly import check_same_run, openmc_name

x = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
o = msgpack.unpackb(Path(sys.argv[2]).read_bytes(), strict_map_key=False)
# 形式、輸送する核種、燃焼点、領域の順（対称操作の組と体積）が同じこと
check_same_run(x, o)
points = [int(a) for a in sys.argv[3:]]
idx = {n: i for i, n in enumerate(x["nuclides"])}


def consumed(n, r, i):
    k = openmc_name(n)
    nx = [p["densities"][r][idx[n]] for p in x["points"]]
    no = o["regions"][r]["densities"][k]
    return ((nx[0] - nx[i]) / (no[0] - no[i]) - 1) * 100


print(
    "組 (行, 列の中央からの距離) | "
    + " | ".join(f"{x['points'][i]['burnup']:g}: U-235 / U-238" for i in points)
)
rows = []
for r, reg in enumerate(o["regions"]):
    vals = [(consumed("U-235", r, i), consumed("U-238", r, i)) for i in points]
    rows.append(vals)
    print(
        f"{tuple(reg['orbit'])} | "
        + " | ".join(f"{a:+.3f} / {b:+.3f}" for a, b in vals)
    )
a = np.array(rows)
for k, i in enumerate(points):
    print(
        f"{x['points'][i]['burnup']:g}: U-235 の範囲 {a[:, k, 0].min():+.3f}〜{a[:, k, 0].max():+.3f}、"
        f"U-238 の範囲 {a[:, k, 1].min():+.3f}〜{a[:, k, 1].max():+.3f}、相関 {np.corrcoef(a[:, k, 0], a[:, k, 1])[0, 1]:+.2f}"
    )
