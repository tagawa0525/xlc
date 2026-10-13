"""VERA Problem 1B〜1D, 2A〜2K, 2Q と MOX ピンセル（examples/mox）の参照解に使う OpenMC の連続エネルギーのデータを、
ENDF/B-VIII.0 から手元の NJOY で作る（xs と同じ評価済みデータから作り、核データの違いを比較から除くため）。

使い方: python3 openmc_library.py neutron <入力> <核種> <出力の HDF5>
        python3 openmc_library.py thermal <入力> <熱中性子散乱則> <出力の HDF5>
        python3 openmc_library.py index <入力> <ライブラリのディレクトリ>

<入力> は openmc-library.toml の形で、中性子の核種ごと（`[neutron.<OpenMC の名前>]`）に ENDF のファイル（`endf`）と
温度（`temperatures`、K）を、熱中性子散乱則ごと（`[thermal.<OpenMC の名前>]`）に中性子と熱の ENDF のファイル
（`neutron`、`thermal`）と温度を書く。パスは入力のファイルからの相対。neutron と thermal は 1 つの核種（散乱則）の
HDF5 を作り、index は <ライブラリのディレクトリ> の `neutron/<名前>.h5` と `thermal/<名前>.h5` を、入力にあるものだけ
並べた `cross_sections.xml` を同じディレクトリに書く。核種ごとに分けるのは、DVC の段階を核種ごとにして、足した核種
だけを作るため（ADR-0130）。`njoy` は PATH から呼ぶ。
"""

import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import openmc.data

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_input import check_keys, input_path, read_input

KINDS = {
    "neutron": {"endf", "temperatures"},
    "thermal": {"neutron", "thermal", "temperatures"},
}


def entry(source: Path, kind: str, name: str) -> dict[str, Any]:
    """入力の `[<kind>.<name>]` の項目。キーを確かめ、温度を float にする。"""
    library = read_input(source, set(), set(KINDS))
    if name not in library.get(kind, {}):
        raise SystemExit(f"{source}: no [{kind}.{name}]")
    item = library[kind][name]
    check_keys(item, f"{source}: [{kind}.{name}]", KINDS[kind], set())
    return {**item, "temperatures": [float(t) for t in item["temperatures"]]}


def check_name(kind: str, expected: str, actual: str) -> None:
    """NJOY から読んだデータの名前が、入力に書いた名前（cross_sections.xml と OpenMC の物質が使う名前）と同じか。"""
    if actual != expected:
        raise SystemExit(f"[{kind}.{expected}] makes data named {actual}")


def neutron(source: Path, name: str, out: Path) -> None:
    item = entry(source, "neutron", name)
    data = openmc.data.IncidentNeutron.from_njoy(
        str(input_path(source, item["endf"])), temperatures=item["temperatures"]
    )
    check_name("neutron", name, data.name)
    out.parent.mkdir(parents=True, exist_ok=True)
    data.export_to_hdf5(str(out), "w")
    print(f"{name}: {item['temperatures']}", flush=True)


def thermal(source: Path, name: str, out: Path) -> None:
    item = entry(source, "thermal", name)
    data = openmc.data.ThermalScattering.from_njoy(
        str(input_path(source, item["neutron"])),
        str(input_path(source, item["thermal"])),
        temperatures=item["temperatures"],
    )
    check_name("thermal", name, data.name)
    out.parent.mkdir(parents=True, exist_ok=True)
    data.export_to_hdf5(str(out), "w")
    print(f"{name}: {item['temperatures']}", flush=True)


def index(source: Path, directory: Path) -> None:
    """`cross_sections.xml` を書く。パスは `directory` からの相対にする（OpenMC は directory の要素がなければ
    cross_sections.xml のあるディレクトリから読む）。`DataLibrary.export_to_xml` は、HDF5 の共通のディレクトリが
    cross_sections.xml のディレクトリと違うと、その絶対パスを書くので使わない（worktree のパスが出力に入る）。"""
    library = read_input(source, set(), set(KINDS))
    root = ET.Element("cross_sections")
    for kind in KINDS:
        for name in library.get(kind, {}):
            path = directory / kind / f"{name}.h5"
            data = openmc.data.DataLibrary()
            data.register_file(path)
            (registered,) = data
            if registered["materials"] != [name]:
                raise SystemExit(f"{path} holds {registered['materials']}, not {name}")
            ET.SubElement(
                root,
                "library",
                materials=name,
                path=f"{kind}/{name}.h5",
                type=registered["type"],
            )
    ET.indent(root)
    ET.ElementTree(root).write(
        directory / "cross_sections.xml", xml_declaration=True, encoding="utf-8"
    )


def main() -> None:
    match sys.argv[1:]:
        case ["neutron", source, name, out]:
            neutron(Path(source), name, Path(out))
        case ["thermal", source, name, out]:
            thermal(Path(source), name, Path(out))
        case ["index", source, directory]:
            index(Path(source), Path(directory))
        case _:
            raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
