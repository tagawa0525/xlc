"""check_natural_isotopes.py（燃焼チェーンと多群ライブラリの入力が、天然の元素の評価を当てる核種で揃っているか。ADR-0262）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_check_natural_isotopes.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_natural_isotopes as cni


def endf_float(x: float) -> str:
    """ENDF の 11 桁の数（指数の e を省く形）"""
    mantissa, exponent = f"{x:.6e}".split("e")
    return f"{mantissa}{int(exponent):+d}".rjust(11)


def write_endf(path: Path, za: float, awr: float = 1.0) -> None:
    """TPID と MF1/MT451 の HEAD だけの評価"""
    head = endf_float(za) + endf_float(awr) + "0".rjust(11) * 4
    path.write_text(
        "a tape identification".ljust(66)
        + "   1 0  0    0\n"
        + head
        + " 125 1451    1\n"
    )


def chain_nuclide(name, z, a):
    return {"name": name, "z": z, "a": a}


class Check(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.neutrons = root / "neutrons"
        self.neutrons.mkdir()
        write_endf(self.neutrons / "n-006_C_000.endf", 6000.0)
        write_endf(self.neutrons / "n-092_U_235.endf", 92235.0)
        write_endf(self.neutrons / "n-095_Am_242m1.endf", 95242.0)
        (root / "chain").mkdir()
        (root / "xs").mkdir()
        self.chain_dir = root / "chain"
        self.xs_dir = root / "xs"
        self.chain_input = {
            "neutron-directory": "../neutrons",
            "natural-isotopes": ["C-12", "C-13"],
        }
        self.chain = {
            "format": "xlc-chain",
            "nuclides": [
                chain_nuclide("C-12", 6, 12),
                chain_nuclide("C-13", 6, 13),
                chain_nuclide("U-235", 92, 235),
                chain_nuclide("Am-242m", 95, 242),
            ],
        }
        self.xs_input = {
            "endf-directory": "../neutrons",
            "nuclide": {
                "U-235": {"endf": "n-092_U_235.endf"},
                "C-12": {"endf": "n-006_C_000.endf"},
                "C-13": {"endf": "n-006_C_000.endf"},
                "Am-242m": {"endf": "n-095_Am_242m1.endf"},
            },
        }

    def tearDown(self):
        self.tmp.cleanup()

    def check(self):
        return cni.check(
            self.chain_input, self.chain_dir, self.chain, self.xs_input, self.xs_dir
        )

    def test_matching_inputs_list_the_natural_evaluations(self):
        result = self.check()
        self.assertEqual(result["errors"], [])
        self.assertEqual(
            result["natural"],
            [["C-12", "n-006_C_000.endf"], ["C-13", "n-006_C_000.endf"]],
        )
        self.assertEqual(result["chain_only"], [])

    def test_natural_isotopes_outside_the_library_are_listed(self):
        del self.xs_input["nuclide"]["C-12"]
        del self.xs_input["nuclide"]["C-13"]
        result = self.check()
        self.assertEqual(result["errors"], [])
        self.assertEqual(result["natural"], [])
        self.assertEqual(result["chain_only"], ["C-12", "C-13"])

    def test_a_natural_evaluation_on_an_isotope_the_chain_does_not_list_is_an_error(
        self,
    ):
        self.chain_input["natural-isotopes"] = ["C-12"]
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_a_listed_isotope_with_its_own_evaluation_in_the_library_is_an_error(self):
        write_endf(self.neutrons / "n-006_C_012.endf", 6012.0)
        self.xs_input["nuclide"]["C-12"] = {"endf": "n-006_C_012.endf"}
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_a_natural_evaluation_of_another_element_is_an_error(self):
        write_endf(self.neutrons / "n-023_V_000.endf", 23000.0)
        self.xs_input["nuclide"]["C-13"] = {"endf": "n-023_V_000.endf"}
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_an_isotope_evaluation_of_another_nuclide_is_an_error(self):
        self.xs_input["nuclide"]["Am-242m"] = {"endf": "n-092_U_235.endf"}
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_a_library_nuclide_outside_the_chain_is_an_error(self):
        self.chain["nuclides"] = self.chain["nuclides"][:-1]
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_different_neutron_directories_are_an_error(self):
        other = Path(self.tmp.name) / "other"
        other.mkdir()
        for f in self.neutrons.iterdir():
            (other / f.name).write_bytes(f.read_bytes())
        self.xs_input["endf-directory"] = "../other"
        self.assertEqual(len(self.check()["errors"]), 1)

    def test_a_tape_without_the_mf1_head_is_an_error(self):
        (self.neutrons / "n-092_U_235.endf").write_text("no records\n")
        with self.assertRaises(SystemExit):
            self.check()


class EndfFloat(unittest.TestCase):
    def test_numbers_without_the_exponent_letter(self):
        self.assertEqual(cni.endf_number(" 9.223500+4"), 92235.0)
        self.assertEqual(cni.endf_number("-1.500000-3"), -1.5e-3)
        self.assertEqual(cni.endf_number(" 6.000000e3"), 6000.0)


if __name__ == "__main__":
    unittest.main()
