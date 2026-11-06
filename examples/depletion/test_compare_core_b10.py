"""compare_core_b10.py の要約（B-10 の割合を変えた状態の core の k と、直接の格子計算の k∞。ADR-0237）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_b10.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_b10


def table(b10: float, ppms: list[float], k) -> dict:
    """燃焼度 0、20 とホウ素 `ppms` の、k∞ が k(ppm, 燃焼度) のテーブル（減速材と燃料は 565 K の 1 点）"""
    burnups = [0.0, 20.0]
    return {
        "axes": [
            {"name": "burnup", "points": burnups},
            {"name": "boron", "points": ppms, "b10": b10},
            {"name": "moderator", "points": [565.0], "density": [0.743]},
            {"name": "fuel", "points": [565.0]},
        ],
        "points": [[u, p, 0, 0] for u in range(2) for p in range(len(ppms))],
        "states": [{"k_inf": k(ppm, u)} for u in burnups for ppm in ppms],
    }


def nominal_k(ppm: float, u: float) -> float:
    return 1.2 - 1e-4 * ppm - 1e-3 * u


def case_and_out(k_of) -> tuple[dict, str]:
    """18 at% の 1300 ppm の状態と、テーブルの割合の同じ ppm の状態のケースと、`k_of(名前, 状態)` の k を書いた標準出力"""
    states = {}
    for u in [0, 20]:
        states[f"b180-p1300-u{u}"] = {"boron": 1300.0, "b10": 0.18, "burnup": u}
        states[f"nominal-p1300-u{u}"] = {"boron": 1300.0, "b10": 0.199, "burnup": u}
    for s in states.values():
        s.update({"moderator": 565, "fuel": 565})
    out = "".join(
        f"state {n}\nk-eff {k_of(n, s):.8f}\n" for n, s in sorted(states.items())
    )
    return {"state": states}, out


class Summary(unittest.TestCase):
    def test_differences_are_reactivities_from_the_direct_k(self):
        base = table(0.199, [600.0, 1300.0], nominal_k)
        direct = table(0.18, [1300.0], lambda p, u: nominal_k(p, u) + 0.01)

        def k_of(name: str, s: dict) -> float:
            k = nominal_k(1300.0, s["burnup"])
            return k if name.startswith("nominal") else k + 0.0098

        case, out = case_and_out(k_of)
        summary = compare_core_b10.summarize(base, case, out, [direct])
        self.assertEqual(summary["table_b10"], 0.199)
        self.assertEqual(
            [(p["b10"], p["boron"], p["burnup"]) for p in summary["points"]],
            [(0.18, 1300.0, 0.0), (0.18, 1300.0, 20.0)],
        )
        for p in summary["points"]:
            k = nominal_k(1300.0, p["burnup"])
            self.assertEqual(
                set(p), {"b10", "boron", "burnup", "direct", "effect", "converted"}
            )
            self.assertAlmostEqual(p["direct"], k + 0.01, places=12)
            self.assertAlmostEqual(
                p["effect"], (1 / k - 1 / (k + 0.01)) * 1e5, places=3
            )
            self.assertAlmostEqual(
                p["converted"], (1 / (k + 0.01) - 1 / (k + 0.0098)) * 1e5, places=3
            )

    def test_states_that_do_not_match_their_names_are_errors(self):
        base = table(0.199, [600.0, 1300.0], nominal_k)
        direct = table(0.18, [1300.0], nominal_k)

        def k_of(n: str, s: dict) -> float:
            return nominal_k(1300.0, s["burnup"])

        for change in [
            lambda c: c["state"]["b180-p1300-u0"].update(b10=0.199),
            lambda c: c["state"]["b180-p1300-u0"].update(boron=1174.0),
            lambda c: c["state"]["nominal-p1300-u20"].update(burnup=0),
            lambda c: c["state"]["nominal-p1300-u20"].update(moderator=580),
        ]:
            case, out = case_and_out(k_of)
            change(case)
            with self.assertRaises(SystemExit):
                compare_core_b10.summarize(base, case, out, [direct])

    def test_nominal_states_must_reproduce_the_table(self):
        base = table(0.199, [600.0, 1300.0], nominal_k)
        direct = table(0.18, [1300.0], nominal_k)
        case, out = case_and_out(
            lambda n, s: (
                nominal_k(1300.0, s["burnup"])
                + (1e-3 if n.startswith("nominal") else 0)
            )
        )
        with self.assertRaises(SystemExit):
            compare_core_b10.summarize(base, case, out, [direct])


if __name__ == "__main__":
    unittest.main()
