"""#0133 の試行: 核種の 1/E の重みの群の捕獲断面積を、xlc のライブラリのポイントワイズと OpenMC の HDF5 で比べる。
使い方: python3 dilute.py <xlc のライブラリ> <OpenMC の cross_sections.xml> <温度> <下端 eV> <上端 eV> <核種 ...>"""

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


lib = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
xml = Path(sys.argv[2])
T = float(sys.argv[3])
lo, hi = float(sys.argv[4]), float(sys.argv[5])
libs = openmc.data.DataLibrary.from_xml(xml)


def avg(e, s):
    m = (e >= lo) & (e <= hi)
    e2 = np.concatenate([[lo], e[m], [hi]])
    s2 = np.interp(e2, e, s)
    # 1/E の重み、点の間はレサジーで台形
    u = np.log(e2)
    return np.trapezoid(s2, u) / (u[-1] - u[0])


for name in sys.argv[6:]:
    nd = next(n for n in lib["nuclides"] if n["name"] == name)
    td = next(t for t in nd["temperatures"] if t["temperature"] == T)
    pw = td["pointwise"]
    e = np.array(pw["energy"])
    c = np.array(pw["capture"])
    xs_x = avg(e, c)
    h5 = libs.get_by_material(openmc_name(name))["path"]
    o = openmc.data.IncidentNeutron.from_hdf5(h5)
    key = next(k for k in o.temperatures if abs(float(k.rstrip("K")) - T) < 1)
    eo = o.energy[key]
    co = o[102].xs[key](eo)
    xs_o = avg(eo, co)
    print(
        f"{name} {T} K {lo}-{hi} eV: xlc pointwise {xs_x:.5g} b, OpenMC {xs_o:.5g} b, 差 {(xs_x / xs_o - 1) * 100:+.3f}%"
    )
