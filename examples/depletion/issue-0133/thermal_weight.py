"""#0133 の試行: 熱群の群定数の群の中の重みを確かめる。核種の点の捕獲断面積（OpenMC の HDF5）を、xlc の重み関数
（NJOY iwt=4 の熱の部分 E exp(−E/0.0253 eV)。0.1 eV より上は 1/E）と OpenMC の燃料の細かい中性子束で群平均し、その比を、xlc と
OpenMC の群の捕獲断面積の比と並べる。
使い方: python3 thermal_weight.py <作業ディレクトリ> <cross_sections.xml> <温度 K> <核種 ...>"""

import re
import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc.data


def openmc_name(name):
    """xlc の核種の名前（Pm-148m）を OpenMC の名前（Pm148_m1）にする"""
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m)?", name).groups()
    return f"{symbol}{mass}" + ("_m1" if isomer else "")


w = Path(sys.argv[1])
libs = openmc.data.DataLibrary.from_xml(sys.argv[2])
T = float(sys.argv[3])
o = msgpack.unpackb((w / "openmc-rates.msgpack").read_bytes())
xs = msgpack.unpackb((w / "xlc-fuel-xs.msgpack").read_bytes())["nuclides"]
fb = np.array(o["fine"]["bounds"])
fphi = np.array(o["fine"]["flux"])
b = np.array(o["bounds"])


def weight(e):
    return np.where(
        e <= 0.1, e * np.exp(-e / 0.0253) / (0.1 * np.exp(-0.1 / 0.0253)) / 0.1, 1 / e
    )


for name in sys.argv[4:]:
    d = openmc.data.IncidentNeutron.from_hdf5(
        libs.get_by_material(openmc_name(name))["path"]
    )
    key = next(k for k in d.temperatures if abs(float(k.rstrip("K")) - T) < 1)
    cap = d[102].xs[key]
    r = o["rates"][name]
    n = o["densities"][name]
    phi = np.array(o["flux"])
    sx = (
        np.array(xs[name.replace("-", "")]["total"])
        - np.array(xs[name.replace("-", "")]["elastic"])
        - np.array(xs[name.replace("-", "")]["fission"])
    )
    so = np.array(r["(n,gamma)"]["mean"]) / phi / n
    print(
        f"{name}: 群 | 下端 eV | xlc/OpenMC − 1 (%) | 重み関数/細かい中性子束 − 1 (%)"
    )
    for g in range(len(b) - 1):
        lo, hi = b[g + 1], b[g]
        if hi > 0.63 or lo < 1.0e-3:
            continue
        m = (fb[:-1] >= lo * 0.9999) & (fb[1:] <= hi * 1.0001)
        # 重み関数の平均は、細かい区間の覆う範囲（群の端の、群をまたぐ区間を除いた範囲）でとる
        e = np.geomspace(fb[:-1][m].min(), fb[1:][m].max(), 2001)
        u = np.log(e)
        sw = np.trapezoid(cap(e) * weight(e) * e, u) / np.trapezoid(weight(e) * e, u)
        # 細かい区間の中では断面積を区間の幾何平均の点の値とする
        mid = np.sqrt(fb[:-1][m] * fb[1:][m])
        sf = (cap(mid) * fphi[m]).sum() / fphi[m].sum()
        print(
            f"  {g + 1} | {lo:.4g} | {(sx[g] / so[g] - 1) * 100:+.3f} | {(sw / sf - 1) * 100:+.3f}"
        )
