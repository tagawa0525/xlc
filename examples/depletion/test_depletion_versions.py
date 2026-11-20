"""版ごとのピンセルの燃焼の比較（#0303）の入力が、本番の版（ENDF/B-VIII.0）の入力と、版の違いのほかで同じであること。

実行: python3 -m unittest discover -s examples/depletion -p 'test_depletion_versions.py'
"""

import copy
import sys
import unittest
from pathlib import Path

import tomllib
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from process_survey import DECLARATION_KEYS

# 版のディレクトリ → ENDF の版の名前。VIII.0 のライブラリとチェーンは HERE の xs.toml と chain.toml
VERSIONS = {
    "endfb70": "ENDF-B-VII.0",
    "endfb71": "ENDF-B-VII.1",
    "endfb80": "ENDF-B-VIII.0",
    "endfb81": "ENDF-B-VIII.1",
}
# tsl-HinH2O.endf に 573.6 K がない版（docs/environment.md）
WITHOUT_5736 = {"endfb70", "endfb71"}


def read(path: Path) -> dict:
    return tomllib.loads(path.read_text())


class Libraries(unittest.TestCase):
    def test_each_library_is_the_production_one_with_its_version_files(self):
        production = read(HERE / "xs.toml")
        for key, name in VERSIONS.items():
            if key == "endfb80":
                continue
            with self.subTest(key):
                ours = read(HERE / key / "xs.toml")
                self.assertEqual(
                    ours.pop("endf-directory"),
                    f"../../../references/endf/{name}_neutrons",
                )
                self.assertEqual(
                    ours.pop("thermal-directory"),
                    f"../../../references/endf/{name}_thermal_scatt",
                )
                expected = copy.deepcopy(production)
                del expected["endf-directory"], expected["thermal-directory"]
                if key in WITHOUT_5736:
                    expected["nuclide"]["H-1"]["temperatures"] = [550, 600]
                # 評価の欠陥を読み替える宣言は版の評価ごとに違う（あることは test_input_declarations.py が確かめる）
                for entry in [*ours["nuclide"].values(), *expected["nuclide"].values()]:
                    for k in DECLARATION_KEYS:
                        entry.pop(k, None)
                self.assertEqual(ours, expected)

    def test_each_chain_reads_the_neutrons_of_its_library(self):
        for key in VERSIONS:
            with self.subTest(key):
                if key == "endfb80":
                    chain, xs = HERE / "chain.toml", HERE / "xs.toml"
                else:
                    chain, xs = HERE / key / "chain.toml", HERE / key / "xs.toml"
                self.assertEqual(
                    (chain.parent / read(chain)["neutron-directory"]).resolve(),
                    (xs.parent / read(xs)["endf-directory"]).resolve(),
                )


class Stages(unittest.TestCase):
    """examples/depletion/dvc.yaml の foreach の組が、入力の読むパスと同じこと（違うと、入力が読むファイルが変わってもステージが走らない）"""

    @classmethod
    def setUpClass(cls):
        cls.stages = yaml.safe_load((HERE / "dvc.yaml").read_text())["stages"]

    def test_depletion_deps_are_the_library_and_chain_of_each_input(self):
        items = self.stages["lattice-1C-version"]["foreach"]
        self.assertEqual(sorted(items), sorted(VERSIONS))
        for key, item in items.items():
            with self.subTest(key):
                d = read(HERE / key / "lattice-1C-60.toml")
                for dep, path in (
                    (item["library"], d["library"]),
                    (item["chain"], d["depletion"]["chain"]),
                ):
                    self.assertEqual(
                        (HERE / dep).resolve(), (HERE / key / path).resolve()
                    )

    def test_check_deps_are_the_neutrons_of_each_chain(self):
        items = self.stages["natural-isotopes"]["foreach"]
        self.assertEqual(sorted(items), sorted(VERSIONS))
        for key, item in items.items():
            with self.subTest(key):
                chain = HERE / f"{item['chain']}.toml"
                self.assertEqual(
                    (chain.parent / read(chain)["neutron-directory"]).resolve(),
                    (
                        HERE / f"../../references/endf/{item['neutrons']}_neutrons"
                    ).resolve(),
                )
                self.assertEqual(item["neutrons"], VERSIONS[key])


class Depletions(unittest.TestCase):
    def test_the_depletions_differ_only_in_the_library_and_chain(self):
        inputs = {k: read(HERE / k / "lattice-1C-60.toml") for k in VERSIONS}
        for key, d in inputs.items():
            with self.subTest(key):
                up = "../" if key == "endfb80" else ""
                self.assertEqual(d.pop("library"), f"{up}depletion.msgpack")
                self.assertEqual(d["depletion"].pop("chain"), f"{up}chain.msgpack")
        first = inputs["endfb70"]
        for key, d in inputs.items():
            with self.subTest(key):
                self.assertEqual(d, first)

    def test_the_depletion_is_the_1c_pin_cell_burned_to_60(self):
        ours = read(HERE / "endfb80" / "lattice-1C-60.toml")
        pin = read(HERE / "lattice-1C-rings1.toml")
        burnup = ours["depletion"].pop("burnup")
        short = pin["depletion"].pop("burnup")
        self.assertEqual(burnup[: len(short)], short)
        self.assertEqual(burnup, sorted(burnup))
        self.assertEqual(burnup[-1], 60)
        self.assertIn(17, burnup)
        self.assertEqual(ours["depletion"].pop("output"), "lattice-1C-60.msgpack")
        pin["depletion"].pop("output")
        self.assertEqual(ours.pop("library"), "../depletion.msgpack")
        pin.pop("library")
        self.assertEqual(ours["depletion"].pop("chain"), "../chain.msgpack")
        pin["depletion"].pop("chain")
        self.assertEqual(ours, pin)


if __name__ == "__main__":
    unittest.main()
