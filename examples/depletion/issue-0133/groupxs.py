"""#0133 の試行: 燃料の核種の群ごとの吸収断面積（全 − 弾性）を xlc と OpenMC（absorption / 中性子束）で比べる。
使い方: python3 groupxs.py <作業ディレクトリ> <核種> [群の下限 (eV)] [群の上限 (eV)]"""

import sys
from pathlib import Path

import msgpack
import numpy as np

work = Path(sys.argv[1])
name = sys.argv[2]
lo = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
hi = float(sys.argv[4]) if len(sys.argv) > 4 else 2e7
omc = msgpack.unpackb((work / "openmc-rates.msgpack").read_bytes())
xs = msgpack.unpackb((work / "xlc-fuel-xs.msgpack").read_bytes())["nuclides"][
    name.replace("-", "")
]
b = np.array(omc["bounds"])
P = sum(np.sum(m["nu-fission"]) for m in omc["materials"].values())
phi = np.array(omc["flux"])
n = omc["densities"][name]
r = omc["rates"][name]
ao = np.array(r["absorption"]["mean"])
ae = np.array(r["absorption"]["std"])
co = np.array(r["(n,gamma)"]["mean"])
sx = np.array(xs["total"]) - np.array(xs["elastic"])
print("群 | 上端 eV | 下端 eV | OpenMC σa (b) | ± % | xlc σa | 差 % | 吸収/生成 ×1e5")
for g in range(len(b) - 1):
    if b[g + 1] < lo or b[g] > hi or ao[g] <= 0:
        continue
    so = ao[g] / phi[g] / n
    print(
        f"{g + 1} | {b[g]:.4g} | {b[g + 1]:.4g} | {so:.5g} | {ae[g] / ao[g] * 100:.2f} | {sx[g]:.5g} | {(sx[g] / so - 1) * 100:+.2f} | {ao[g] / P * 1e5:.3f}"
    )
