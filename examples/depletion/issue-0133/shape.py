"""#0133 の試行: ポイントワイズの形の比較。NR の重み 1/(σt + σ0)/E の捕獲の群平均を xlc のライブラリと OpenMC の HDF5 で比べる。
使い方: python3 shape.py <xlc のライブラリ> <cross_sections.xml> <温度> <σ0> <下端> <上端> <核種 ...>"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc.data

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_chain import openmc_name

lib = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), strict_map_key=False)
libs = openmc.data.DataLibrary.from_xml(Path(sys.argv[2]))
T, s0, lo, hi = map(float, sys.argv[3:7])


def avg(e, c, t):
    grid = np.geomspace(lo, hi, 400001)
    w = 1.0 / (np.interp(grid, e, t) + s0)
    u = np.log(grid)
    return np.trapezoid(np.interp(grid, e, c) * w, u) / np.trapezoid(w, u)


for name in sys.argv[7:]:
    nd = next(n for n in lib["nuclides"] if n["name"] == name)
    pw = next(t for t in nd["temperatures"] if t["temperature"] == T)["pointwise"]
    x = avg(np.array(pw["energy"]), np.array(pw["capture"]), np.array(pw["total"]))
    o = openmc.data.IncidentNeutron.from_hdf5(
        libs.get_by_material(openmc_name(name))["path"]
    )
    key = next(k for k in o.temperatures if abs(float(k.rstrip("K")) - T) < 1)
    eo = o.energy[key]
    y = avg(eo, o[102].xs[key](eo), o[1].xs[key](eo))
    print(
        f"{name} σ0 {s0} {lo}-{hi}: xlc {x:.6g} OpenMC {y:.6g} 差 {(x / y - 1) * 100:+.3f}%"
    )
