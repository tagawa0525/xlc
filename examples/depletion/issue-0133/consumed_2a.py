"""#0133 の試行: 2A の燃焼点ごとに、U-235 の消えた量、U-236 の増えた量（U-235 の捕獲の積分。U-236 自身の吸収を含む）、
U-238 の消えた量、安定な核分裂生成物（Ru-101、Mo-95。核分裂の数の積分）を xlc/OpenMC − 1（%）で比べる。

使い方: python3 consumed_2a.py <lattice-2A.msgpack> <openmc-2A.msgpack>
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
idx = {n: i for i, n in enumerate(x["nuclides"])}
w = np.array([r["volume"] for r in x["regions"]])
w /= w.sum()
bu = [p["burnup"] for p in x["points"]]


def pair(n):
    k = openmc_name(n)
    ax = (
        np.array(
            [[p["densities"][r][idx[n]] for r in range(len(w))] for p in x["points"]]
        )
        @ w
    )
    ao = np.array([o["regions"][r]["densities"][k] for r in range(len(w))]).T @ w
    return ax, ao


u5, u6, u8, ru, mo = (pair(n) for n in ["U-235", "U-236", "U-238", "Ru-101", "Mo-95"])
print(
    "燃焼度 | U-235 消費 | U-236 増 | U-238 消費 | Ru-101 | Mo-95 | U-235 消費/OpenMC の U-235 初期 (%)"
)


def f(a, i):
    return (a[0][i] - a[0][0]) / (a[1][i] - a[1][0]) * 100 - 100


def g(a, i):
    return (a[0][i] / a[1][i] - 1) * 100


for i in range(1, len(bu)):
    print(
        f"{bu[i]:g} | {f(u5, i):+.3f} | {f(u6, i):+.3f} | {f(u8, i):+.3f} | {g(ru, i):+.3f} | {g(mo, i):+.3f} | {(u5[1][0] - u5[1][i]) / u5[1][0] * 100:.1f}"
    )
