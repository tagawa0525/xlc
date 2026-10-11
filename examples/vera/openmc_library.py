"""VERA Problem 1B〜1D, 2A〜2K, 2Q と MOX ピンセル（examples/mox）の参照解に使う OpenMC の連続エネルギーのデータを、
ENDF/B-VIII.0 から手元の NJOY で作る（xs と同じ評価済みデータから作り、核データの違いを比較から除くため）。

使い方: python3 openmc_library.py <ENDF/B-VIII.0 のディレクトリ> <出力ディレクトリ> <並列数>
            --fuel <K>... --coolant-and-clad <K>... --others <K>... --thermal <K>...

`njoy` は PATH から呼ぶ。出力ディレクトリに核種ごとの HDF5 と cross_sections.xml を書く。温度は核種の組ごとに渡す:
`--fuel` は燃料の核種と O-16、`--coolant-and-clad` は冷却材と被覆管（Zircaloy-4）の核種（H, B, Zr, Sn, Fe, Cr）、
`--others` はそのほか（Pyrex、SS304、AIC と MOX の Pu と Am）、`--thermal` は H in H2O の S(α,β)。H-1 は
`--coolant-and-clad` と `--thermal` の両方の温度で作る（OpenMC は内挿のモードで物質の温度そのものの H-1 を使い、
S(α,β) の温度の H-1 は、xlc が内挿する H-1 との違いを見積もるのに使う。compare_h1_temperatures.py）。
OpenMC の温度の内挿のモードは、温度を 1 つしか持たない核種が 1 つでもあると全体を最も近い温度のモードに戻す
（OpenMC 0.15.2 の src/nuclide.cpp）ので、内挿のモードで解く問題の核種には 2 つ以上の温度を渡す。

出力ディレクトリにすでにある HDF5 は、求める温度をすべて含めば作り直さない（核種を足したときに、足した分だけを
作る）。含まなければ、求める温度だけで作り直す（持っていた温度でも、渡さなければ消える）ので、温度を足すときは
持っている温度も並べて渡す。HDF5 は一時ファイルに書いてから名前を変えるので、途中で止まっても書きかけのファイルは残らない。
"""

import argparse
import os
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import openmc.data

# (名前, ENDF ファイル名)。Hf と He は ADR-0054 と ADR-0120 のとおり除く
FUEL = [
    ("U234", "n-092_U_234"),
    ("U235", "n-092_U_235"),
    ("U236", "n-092_U_236"),
    ("U238", "n-092_U_238"),
    ("O16", "n-008_O_016"),
]
HYDROGEN = ("H1", "n-001_H_001")
# 冷却材（H は HYDROGEN）と被覆管の核種
COOLANT_AND_CLAD = (
    [("B10", "n-005_B_010"), ("B11", "n-005_B_011")]
    + [(f"Zr{a}", f"n-040_Zr_{a:03d}") for a in (90, 91, 92, 94, 96)]
    + [
        (f"Sn{a}", f"n-050_Sn_{a:03d}")
        for a in (112, 114, 115, 116, 117, 118, 119, 120, 122, 124)
    ]
    + [(f"Fe{a}", f"n-026_Fe_{a:03d}") for a in (54, 56, 57, 58)]
    + [(f"Cr{a}", f"n-024_Cr_{a:03d}") for a in (50, 52, 53, 54)]
)
OTHERS = (
    [("C12", "n-006_C_012"), ("C13", "n-006_C_013"), ("P31", "n-015_P_031")]
    + [(f"Si{a}", f"n-014_Si_{a:03d}") for a in (28, 29, 30)]
    + [("Mn55", "n-025_Mn_055")]
    + [(f"Ni{a}", f"n-028_Ni_{a:03d}") for a in (58, 60, 61, 62, 64)]
    + [(f"Pu{a}", f"n-094_Pu_{a:03d}") for a in (238, 239, 240, 241, 242)]
    + [("Am241", "n-095_Am_241")]
    + [(f"Ag{a}", f"n-047_Ag_{a:03d}") for a in (107, 109)]
    + [(f"In{a}", f"n-049_In_{a:03d}") for a in (113, 115)]
    + [
        (f"Cd{a}", f"n-048_Cd_{a:03d}")
        for a in (106, 108, 110, 111, 112, 113, 114, 116)
    ]
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


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="OpenMC data for the VERA references")
    parser.add_argument("endf", type=Path, help="ENDF/B-VIII.0 のディレクトリ")
    parser.add_argument("out", type=Path, help="出力ディレクトリ")
    parser.add_argument("workers", type=int, help="並列数")
    for kind in ("fuel", "coolant-and-clad", "others", "thermal"):
        parser.add_argument(f"--{kind}", type=float, nargs="+", required=True)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    endf, out = args.endf.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    neutrons = endf / "ENDF-B-VIII.0_neutrons"
    hydrogen = sorted({*args.coolant_and_clad, *args.thermal})
    kinds = [
        (FUEL, args.fuel),
        ([HYDROGEN], hydrogen),
        (COOLANT_AND_CLAD, args.coolant_and_clad),
        (OTHERS, args.others),
    ]
    jobs = [
        (n, neutrons / f"{f}.endf", sorted(temperatures), out)
        for nuclides, temperatures in kinds
        for n, f in nuclides
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        paths = list(pool.map(neutron, jobs))

    thermal_path = out / "c_H_in_H2O.h5"
    thermal_temperatures = sorted(args.thermal)
    if not has_temperatures(thermal_path, thermal_temperatures):
        thermal = openmc.data.ThermalScattering.from_njoy(
            str(neutrons / "n-001_H_001.endf"),
            str(endf / "ENDF-B-VIII.0_thermal_scatt" / "tsl-HinH2O.endf"),
            temperatures=thermal_temperatures,
        )
        export(thermal_path, lambda p: thermal.export_to_hdf5(p, "w"))
        print(f"{thermal.name}: {thermal_temperatures}", flush=True)

    library = openmc.data.DataLibrary()
    for path in [*paths, thermal_path]:
        library.register_file(str(path))
    library.export_to_xml(str(out / "cross_sections.xml"))


if __name__ == "__main__":
    main()
