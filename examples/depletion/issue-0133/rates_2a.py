"""#0133 の試行: 2A の xlc と OpenMC の燃焼の履歴の数密度から、ステップごとの 1 原子あたりの消滅率 ln(N_i/N_{i+1})/Δt を
核種ごとに出し、xlc/OpenMC − 1 を比べる（生成の項は無視するので、U-235、U-238 のように生成の小さい核種で読む）。

使い方: python3 rates_2a.py <lattice-2A.msgpack> <openmc-2A.msgpack> [核種 ...]
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
names = sys.argv[3:] or [
    "U-235",
    "U-238",
    "Pu-239",
    "Pu-240",
    "Pu-241",
    "Pu-242",
    "U-236",
    "Am-241",
    "Np-237",
]
idx = {n: i for i, n in enumerate(x["nuclides"])}
w = np.array([r["volume"] for r in x["regions"]])
w = w / w.sum()
t = np.array([p["time"] for p in x["points"]])
bu = np.array([p["burnup"] for p in x["points"]])


def dens(n):
    nx = np.array(
        [[p["densities"][r][idx[n]] for r in range(len(w))] for p in x["points"]]
    )
    no = np.array(
        [o["regions"][r]["densities"][openmc_name(n)] for r in range(len(w))]
    ).T
    return nx, no


res = {}
for n in names:
    nx, no = dens(n)
    ax, ao = nx @ w, no @ w
    res[n] = (ax, ao)
print("燃焼度 | " + " | ".join(f"{n} N差% / 率差%" for n in names))
for i in range(len(t) - 1):
    row = [f"{bu[i]:g}-{bu[i + 1]:g}"]
    for n in names:
        ax, ao = res[n]
        dn = (ax[i + 1] / ao[i + 1] - 1) * 100 if ao[i + 1] > 0 else float("nan")
        if ax[i] > 0 and ao[i] > 0:
            rx = np.log(ax[i] / ax[i + 1])
            ro = np.log(ao[i] / ao[i + 1])
            dr = (rx / ro - 1) * 100
        else:
            dr = float("nan")
        row.append(f"{dn:+.3f} / {dr:+.3f}")
    print(" | ".join(row))
