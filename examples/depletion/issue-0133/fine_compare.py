"""#0133 の試行: UFG の等価なセルの燃料の細かい中性子束（fine_flux の例）を、OpenMC の燃料の細かい tally と、群ごとに規格化して比べる。
使い方: python3 fine_compare.py <作業ディレクトリ（xlc-fine.msgpack、openmc-rates.msgpack）> <群の下端 eV> <群の上端 eV> [区間の数]"""

import itertools
import sys
from pathlib import Path

import msgpack
import numpy as np

w = Path(sys.argv[1])
lo, hi = float(sys.argv[2]), float(sys.argv[3])
nsub = int(sys.argv[4]) if len(sys.argv) > 4 else 10
x = msgpack.unpackb((w / "xlc-fine.msgpack").read_bytes())
o = msgpack.unpackb((w / "openmc-rates.msgpack").read_bytes())["fine"]
e = np.array(x["energies"])[::-1]
F = np.array(x["fuel"])[::-1]
u = np.log(e)
b = np.array(o["bounds"])
phi = np.array(o["flux"])


def xlc_integral(a, c):
    # レサジーについて台形で ∫F du（点の間は線形）
    grid = np.concatenate(
        [[np.log(a)], u[(u > np.log(a)) & (u < np.log(c))], [np.log(c)]]
    )
    return np.trapezoid(np.interp(grid, u, F), grid)


inside = (b[:-1] >= lo * 0.9999) & (b[1:] <= hi * 1.0001)
idx = np.where(inside)[0]
px = np.array([xlc_integral(b[i], b[i + 1]) for i in idx])
po = phi[idx]
px /= px.sum()
po /= po.sum()
print(f"{lo}-{hi} eV: {len(idx)} bins")
edges = np.linspace(0, len(idx), nsub + 1).astype(int)
for a, c in itertools.pairwise(edges):
    print(
        f"  {b[idx[a]]:.4g}-{b[idx[c - 1] + 1]:.4g} eV: xlc {px[a:c].sum():.5f} OpenMC {po[a:c].sum():.5f} 比 {px[a:c].sum() / po[a:c].sum() - 1:+.4f}"
    )
for name, cap in o["capture"].items():
    c = np.array(cap)[idx]
    if c.sum() <= 0:
        continue
    s = c / phi[idx]  # 区間の σ × N（OpenMC の重み）
    print(
        f"  {name}: 区間の断面積を xlc の中性子束で重みづけた群の捕獲 / OpenMC − 1 = {(s * px).sum() / (s * po).sum() - 1:+.4f}"
    )
