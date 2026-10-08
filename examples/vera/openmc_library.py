"""VERA Problem 1B〜1D, 2B〜2D と MOX ピンセル（examples/mox）の参照解に使う OpenMC の連続エネルギーのデータを、ENDF/B-VIII.0 から
手元の NJOY で作る（xs と同じ評価済みデータから作り、核データの違いを比較から除くため）。

使い方: python3 openmc_library.py <ENDF/B-VIII.0 のディレクトリ> <出力ディレクトリ> [並列数]

`njoy` は PATH から呼ぶ。出力ディレクトリに核種ごとの HDF5 と cross_sections.xml を書く。
燃料の核種は 600, 900, 1200 K、ほかの核種（MOX の Pu と Am を含む）と H in H2O の S(α,β) は 600 K だけを作る。
出力ディレクトリにすでにある HDF5 は、求める温度をすべて含めば作り直さない（核種を足したときに、足した分だけを
作る）。HDF5 は一時ファイルに書いてから名前を変えるので、途中で止まっても書きかけのファイルは残らない。
"""

import os
import sys
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import openmc.data

FUEL_TEMPERATURES = [600.0, 900.0, 1200.0]
OTHER_TEMPERATURES = [600.0]

# (名前, ENDF ファイル名)。Hf と He は ADR-0054 のとおり除く
FUEL = [
    ("U234", "n-092_U_234"),
    ("U235", "n-092_U_235"),
    ("U236", "n-092_U_236"),
    ("U238", "n-092_U_238"),
    ("O16", "n-008_O_016"),
]
OTHERS = (
    [("H1", "n-001_H_001"), ("B10", "n-005_B_010"), ("B11", "n-005_B_011")]
    + [(f"Zr{a}", f"n-040_Zr_{a:03d}") for a in (90, 91, 92, 94, 96)]
    + [
        (f"Sn{a}", f"n-050_Sn_{a:03d}")
        for a in (112, 114, 115, 116, 117, 118, 119, 120, 122, 124)
    ]
    + [(f"Fe{a}", f"n-026_Fe_{a:03d}") for a in (54, 56, 57, 58)]
    + [(f"Cr{a}", f"n-024_Cr_{a:03d}") for a in (50, 52, 53, 54)]
    + [(f"Pu{a}", f"n-094_Pu_{a:03d}") for a in (238, 239, 240, 241, 242)]
    + [("Am241", "n-095_Am_241")]
)


def has_temperatures(path: Path, temperatures: list[float]) -> bool:
    """`path` の HDF5 が、求める温度（OpenMC は "600K" のような名前で持つ）をすべて含むか。"""
    if not path.exists():
        return False
    with h5py.File(path, "r") as f:
        (name,) = f.keys()
        stored = set(f[name]["kTs"].keys())
    return {f"{round(t)}K" for t in temperatures} <= stored


def export(path: Path, write: Callable[[str], None]) -> None:
    """`write`（HDF5 の出力先を受け取る）で一時ファイルに書き、成功したら `path` に置き換える。"""
    temporary = path.with_suffix(".h5.tmp")
    write(str(temporary))
    os.replace(temporary, path)


def neutron(args: tuple[str, Path, list[float], Path]) -> Path:
    name, endf, temperatures, out = args
    path = out / f"{name}.h5"
    if has_temperatures(path, temperatures):
        return path
    data = openmc.data.IncidentNeutron.from_njoy(str(endf), temperatures=temperatures)
    export(path, lambda p: data.export_to_hdf5(p, "w"))
    print(f"{name}: {temperatures}", flush=True)
    return path


def main() -> None:
    endf = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    out.mkdir(parents=True, exist_ok=True)
    neutrons = endf / "ENDF-B-VIII.0_neutrons"
    jobs = [(n, neutrons / f"{f}.endf", FUEL_TEMPERATURES, out) for n, f in FUEL] + [
        (n, neutrons / f"{f}.endf", OTHER_TEMPERATURES, out) for n, f in OTHERS
    ]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        paths = list(pool.map(neutron, jobs))

    thermal_path = out / "c_H_in_H2O.h5"
    if not has_temperatures(thermal_path, OTHER_TEMPERATURES):
        thermal = openmc.data.ThermalScattering.from_njoy(
            str(neutrons / "n-001_H_001.endf"),
            str(endf / "ENDF-B-VIII.0_thermal_scatt" / "tsl-HinH2O.endf"),
            temperatures=OTHER_TEMPERATURES,
        )
        export(thermal_path, lambda p: thermal.export_to_hdf5(p, "w"))
        print(f"{thermal.name}: {OTHER_TEMPERATURES}", flush=True)

    library = openmc.data.DataLibrary()
    for path in [*paths, thermal_path]:
        library.register_file(str(path))
    library.export_to_xml(str(out / "cross_sections.xml"))


if __name__ == "__main__":
    main()
