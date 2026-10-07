"""VERA Problem 1B〜1D, 2B〜2D の参照解に使う OpenMC の連続エネルギーのデータを、ENDF/B-VIII.0 から
手元の NJOY で作る（xs と同じ評価済みデータから作り、核データの違いを比較から除くため）。

使い方: python3 openmc_library.py <ENDF/B-VIII.0 のディレクトリ> <出力ディレクトリ> [並列数]

`njoy` は PATH から呼ぶ。出力ディレクトリに核種ごとの HDF5 と cross_sections.xml を書く。
燃料の核種は 600, 900, 1200 K、ほかの核種と H in H2O の S(α,β) は 600 K だけを作る。
"""

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import openmc.data

FUEL_TEMPERATURES = [600.0, 900.0, 1200.0]
OTHER_TEMPERATURES = [600.0]

# (名前, ENDF ファイル名)。Hf と He は計画のとおり除く
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
)


def neutron(args: tuple[str, Path, list[float], Path]) -> Path:
    name, endf, temperatures, out = args
    path = out / f"{name}.h5"
    data = openmc.data.IncidentNeutron.from_njoy(str(endf), temperatures=temperatures)
    data.export_to_hdf5(str(path), "w")
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

    thermal = openmc.data.ThermalScattering.from_njoy(
        str(neutrons / "n-001_H_001.endf"),
        str(endf / "ENDF-B-VIII.0_thermal_scatt" / "tsl-HinH2O.endf"),
        temperatures=OTHER_TEMPERATURES,
    )
    thermal_path = out / "c_H_in_H2O.h5"
    thermal.export_to_hdf5(str(thermal_path), "w")
    print(f"{thermal.name}: {OTHER_TEMPERATURES}", flush=True)

    library = openmc.data.DataLibrary()
    for path in [*paths, thermal_path]:
        library.register_file(str(path))
    library.export_to_xml(str(out / "cross_sections.xml"))


if __name__ == "__main__":
    main()
