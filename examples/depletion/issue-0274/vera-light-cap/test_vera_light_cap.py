"""vera_light_cap.py の比較と、置き場の入力（VERA の 1B〜1D、2A〜2K、2Q のコピー。ADR-0246）。

実行: cd examples/depletion/issue-0274/vera-light-cap && python3 -m unittest test_vera_light_cap
"""

import sys
import unittest
from pathlib import Path
from typing import ClassVar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vera_light_cap as vlc

VERA = HERE.parents[2] / "vera"


def result(diff: float, rms: float | None = None, std: float = 2e-5) -> dict:
    pins = (
        None if rms is None else {"positions": 39, "rms": rms, "max": rms, "at": [1, 1]}
    )
    return {"k": {"diff_pcm": diff, "keno_std": std}, "pins": pins}


def results_with(diffs: dict) -> tuple[dict, dict]:
    """問題ごとに (記録, 400-heavy400, 400, e1000) の差を与えた要約。ほかの問題は差 0"""
    results, recorded = {}, {}
    for p in vlc.PINCELLS + vlc.ASSEMBLIES:
        rec, both, light, default = diffs.get(p, (0.0, 0.0, 0.0, 0.0))
        rms = None if p in vlc.PINCELLS else 0.001
        recorded[p] = result(rec, rms)
        results[p, "400-heavy400"] = result(both, rms)
        results[p, "400"] = result(light, rms)
        results[p, "e1000"] = result(default, rms)
    return results, recorded


class CompareTest(unittest.TestCase):
    def test_the_shift_and_the_reproduction_are_differences_of_the_diffs(self):
        results, recorded = results_with({"2A": (185.2, 185.4, 170.0, 150.0)})
        row = next(
            r for r in vlc.problem_rows(results, recorded) if r["problem"] == "2A"
        )
        self.assertAlmostEqual(row["reproduction_pcm"], 0.2)
        self.assertAlmostEqual(row["shift_pcm"], -20.0)
        self.assertAlmostEqual(row["keno_std_pcm"], 2.0)
        self.assertAlmostEqual(row["pins_rms_percent"]["e1000"], 0.1)

    def test_pincells_have_no_pin_powers(self):
        results, recorded = results_with({})
        row = next(
            r for r in vlc.problem_rows(results, recorded) if r["problem"] == "1B"
        )
        self.assertIsNone(row["record_pins_rms_percent"])
        self.assertIsNone(row["pins_rms_percent"]["e1000"])

    def test_totals_are_taken_over_pincells_and_assemblies_apart(self):
        results, recorded = results_with(
            {"1B": (0, 0, 30.0, 0), "1C": (0, 0, -30.0, 0), "2A": (0, 0, 120.0, 0)}
        )
        totals = vlc.aggregate(vlc.problem_rows(results, recorded))
        self.assertAlmostEqual(totals["pincells"]["400"]["mean_pcm"], 0.0)
        self.assertAlmostEqual(totals["pincells"]["400"]["rms_pcm"], (1800 / 3) ** 0.5)
        self.assertAlmostEqual(totals["assemblies"]["400"]["mean_pcm"], 10.0)
        self.assertAlmostEqual(
            totals["assemblies"]["400"]["rms_pcm"], (14400 / 12) ** 0.5
        )


class InputsTest(unittest.TestCase):
    """コピーは元の入力と、先頭の注釈、ライブラリ、ピン出力の名前、上限の行だけが違う"""

    CAPS: ClassVar[dict[str, list[str]]] = {
        "e1000": [],
        "400": ["ufg-free-gas = { upper = 400 }"],
        "400-heavy400": ["ufg-free-gas = { upper = 400, heavy = 400 }"],
    }

    def test_each_copy_differs_from_the_original_only_in_the_library_pins_and_cap(self):
        for p in vlc.PINCELLS + vlc.ASSEMBLIES:
            rel = (
                f"p1/pincell-{p}.toml" if p in vlc.PINCELLS else f"p2/assembly-{p}.toml"
            )
            original = (VERA / rel).read_text().splitlines()
            for v in vlc.VARIANTS:
                with self.subTest(problem=p, variant=v):
                    copy = (HERE / f"{p}-{v}.toml").read_text().splitlines()
                    self.assertTrue(
                        copy[0].startswith(
                            f"# #0274 の VERA の確認: ../../../vera/{rel} のコピーで"
                        )
                    )
                    expected = []
                    for line in original:
                        if line == 'library = "../vera.msgpack"':
                            expected.append(
                                'library = "../../../vera/endfb70/vera.msgpack"'
                            )
                            expected.extend(self.CAPS[v])
                        elif line == f'pin-powers = "pins-{p}.msgpack"':
                            expected.append(f'pin-powers = "pins-{p}-{v}.msgpack"')
                        else:
                            expected.append(line)
                    self.assertEqual(copy[1:], expected)
                    self.assertNotIn('library = "../vera.msgpack"', copy)
                    if p in vlc.ASSEMBLIES:
                        self.assertIn(f'pin-powers = "pins-{p}-{v}.msgpack"', copy)

    def test_the_originals_leave_the_free_gas_to_the_defaults(self):
        # コピーの e1000 が既定（ADR-0246）で解くのは、元の入力が自由ガスの条件を書かないため
        for p in vlc.PINCELLS + vlc.ASSEMBLIES:
            rel = (
                f"p1/pincell-{p}.toml" if p in vlc.PINCELLS else f"p2/assembly-{p}.toml"
            )
            self.assertNotIn("ufg-", (VERA / rel).read_text(), rel)


if __name__ == "__main__":
    unittest.main()
