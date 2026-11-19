"""#0293 の項目 3: 出力の規格化に使う核分裂あたりのエネルギーの定義を、評価済みデータから比べる。

2 つの使い方がある:

- `heating <cross_sections.xml> <核種>...`: OpenMC のライブラリ（HDF5）から、0.0253 eV の値を出す（核分裂しない核種はライブラリの温度ごと。核分裂の値は温度に依らない）。核分裂する核種は、チェーンの Q に
  当たる MT18 の QI、`heating-local` の核分裂の反応の分（openmc/data/neutron.py が MT18 の KERMA を置き換える EFR + EGP + EGD + EB）、
  即発中性子の運動エネルギー ENP、MF1 MT458 の ER と、ER / QI − 1、核分裂の分 / QI − 1。核分裂しない核種は heating-local / 捕獲の断面積
  （散乱の反跳と、(n,α) などの捕獲のほかの反応の分も分子に入る）
- `kappa <ENDF の中性子のディレクトリ> <casl-kappa.toml> <核種>...`: VERA の燃焼のベンチマークの仕様書の Kappa と Kappa0（J）を、その
  評価の MF3 MT18 の QI と、MT458 の ET（ニュートリノを含む全量。0.0253 eV）と比べる

使い方: python3 fission_energy.py heating <cross_sections.xml> <核種>...
        python3 fission_energy.py kappa <ENDF のディレクトリ> <casl-kappa.toml> <核種>...
"""

import argparse
import re
import sys
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path

import tomllib

# CODATA 2018 の電子ボルト（J）。仕様書の表は J で書く
EV_IN_JOULES = 1.602176634e-19
THERMAL = 0.0253


def fission_row(qi: float, er: float, local: float) -> dict:
    """QI、ER、heating-local の核分裂の反応の分（eV）の比。"""
    return {"er_over_qi": er / qi - 1, "local_over_qi": local / qi - 1}


def kappa_row(qi: float, et: float, kappa: float, kappa0: float) -> dict:
    """QI と ET（eV）と、仕様書の Kappa と Kappa0（J）の比。"""
    return {
        "kappa0_over_qi": kappa0 / EV_IN_JOULES / qi - 1,
        "kappa_over_et": kappa / EV_IN_JOULES / et - 1,
        "kappa_over_kappa0": kappa / kappa0,
    }


def heating(cross_sections: Path, nuclides: list[str]) -> None:
    import openmc.data

    root = ET.parse(cross_sections).getroot()
    paths = {
        m.get("materials"): cross_sections.parent / m.get("path")
        for m in root
        if m.get("type") == "neutron"
    }
    for name in nuclides:
        if name not in paths:
            sys.exit(f"{name} is not in {cross_sections}")
        data = openmc.data.IncidentNeutron.from_hdf5(paths[name])
        if 18 not in data.reactions:
            # `temperatures` は HDF5 の群の名前の順（"1200K" が "600K" より前）なので、温度の値で並べる
            for t in sorted(data.temperatures, key=lambda t: float(t.rstrip("K"))):
                per_capture = data[901].xs[t](THERMAL) / data[102].xs[t](THERMAL)
                print(
                    f"{name} {t}: heating-local / capture {per_capture / 1e6:.4f} MeV"
                )
            continue
        f = data.fission_energy
        if f is None:
            sys.exit(f"{name} has no MF1 MT458 fission energy release")
        qi = data[18].q_value
        er = f.recoverable(THERMAL)
        part = (
            f.fragments(THERMAL)
            + f.prompt_photons(THERMAL)
            + f.delayed_photons(THERMAL)
            + f.betas(THERMAL)
        )
        row = fission_row(qi, er, part)
        print(
            f"{name}: QI {qi / 1e6:.4f} MeV, fission heating {part / 1e6:.4f} MeV,"
            f" ENP {f.prompt_neutrons(THERMAL) / 1e6:.4f} MeV, ER {er / 1e6:.4f} MeV,"
            f" ER / QI - 1 {row['er_over_qi'] * 100:+.4f}%,"
            f" fission heating / QI - 1 {row['local_over_qi'] * 100:+.4f}%"
        )


def endf_file(directory: Path, name: str) -> Path:
    """OpenMC の名前（U235）の ENDF の中性子のファイル（n-092_U_235.endf）。"""
    match = re.fullmatch(r"([A-Z][a-z]?)(\d+)", name)
    if match is None:
        sys.exit(f"{name} is not an OpenMC ground-state nuclide name")
    element, mass = match.groups()
    import openmc.data

    z = openmc.data.ATOMIC_NUMBER[element]
    path = directory / f"n-{z:03d}_{element}_{int(mass):03d}.endf"
    if not path.exists():
        sys.exit(f"{path} does not exist")
    return path


def kappa(directory: Path, table: Path, nuclides: list[str]) -> None:
    import openmc.data
    from openmc.data.endf import get_cont_record, get_head_record

    values = tomllib.loads(table.read_text())
    for name in nuclides:
        if name not in values:
            sys.exit(f"{name} is not in {table}")
        evaluation = openmc.data.Evaluation(endf_file(directory, name))
        section = StringIO(evaluation.section[3, 18])
        get_head_record(section)
        qi = get_cont_record(section)[1]
        if (1, 458) not in evaluation.section:
            sys.exit(f"{name} has no MF1 MT458 in {directory}")
        release = openmc.data.FissionEnergyRelease.from_endf(
            evaluation, openmc.data.IncidentNeutron.from_endf(evaluation)
        )
        et = release.q_total(THERMAL)
        k = values[name]
        row = kappa_row(qi, et, k["kappa"], k["kappa0"])
        print(
            f"{name}: QI {qi / 1e6:.4f} MeV, ET {et / 1e6:.4f} MeV,"
            f" kappa0 {k['kappa0'] / EV_IN_JOULES / 1e6:.4f} MeV, kappa {k['kappa'] / EV_IN_JOULES / 1e6:.4f} MeV,"
            f" kappa0 / QI - 1 {row['kappa0_over_qi'] * 100:+.4f}%,"
            f" kappa / ET - 1 {row['kappa_over_et'] * 100:+.4f}%,"
            f" kappa / kappa0 {row['kappa_over_kappa0']:.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    h = sub.add_parser("heating")
    h.add_argument("cross_sections", type=Path)
    h.add_argument("nuclides", nargs="+")
    k = sub.add_parser("kappa")
    k.add_argument("directory", type=Path)
    k.add_argument("table", type=Path)
    k.add_argument("nuclides", nargs="+")
    args = parser.parse_args()
    if args.command == "heating":
        heating(args.cross_sections, args.nuclides)
    else:
        kappa(args.directory, args.table, args.nuclides)


if __name__ == "__main__":
    main()
