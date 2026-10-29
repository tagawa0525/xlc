"""compare_core_restart.py の要約（restart した燃焼と切らない燃焼の差、装荷替えの 1 サイクルの点ごとの k、ホウ素、集合体の燃焼度）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_restart.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_restart

NUCLIDES = ["U-235", "Xe-135"]


def core(points, assemblies=("A01", "A02"), auxiliary=1.0e-9):
    """軸方向のノード 2 つのミクロ燃焼の core の履歴。`points` は (EFPD, 日数, 状態, k, ホウ素, [集合体ごとの燃焼度])。ノードの燃焼度は
    集合体の燃焼度の 0.9 倍と 1.1 倍、数密度は U-235 が 7e-4 − 1e-5 × 燃焼度、Xe-135 が 2e-9"""
    n = len(assemblies)
    return {
        "format": "xlc-core-history",
        "specific_power": 38.41,
        "states": {
            "aro": {
                "search": {},
                "b10": 0.199,
                "moderator": 565.0,
                "fuel": 565.0,
                "banks": {},
            }
        },
        "axial_mesh": [0.0, 10.0, 20.0],
        "assemblies": list(assemblies),
        "heavy_metal": [[4000.0, 4000.0]] * n,
        "nuclides": list(NUCLIDES),
        "midpoint_boron": "search",
        "xenon": "equilibrium",
        "auxiliary": ["I-135"],
        "auxiliary_densities": [[[auxiliary], [auxiliary]]] * n,
        "points": [
            {
                "efpd": e,
                "days": d,
                **({} if k == 0 else {"power": 1.0}),
                "state": s,
                "burnup": sum(b) / n,
                "k": kk,
                "boron": ppm,
                "solves": 5,
                "assembly_burnups": list(b),
                "burnups": [[0.9 * x, 1.1 * x] for x in b],
                "powers": [[0.5 / n, 0.5 / n]] * n,
                "densities": [
                    [[7e-4 - 1e-5 * y, 2e-9] for y in (0.9 * x, 1.1 * x)] for x in b
                ],
            }
            for k, (e, d, s, kk, ppm, b) in enumerate(points)
        ],
    }


FULL = [
    (0.0, 0.0, "aro", 1.071, 1288.9, [0.0, 0.0]),
    (10.0, 10.0, "aro", 1.071, 900.0, [0.4, 0.36]),
    (20.0, 20.0, "aro", 1.071, 800.0, [0.8, 0.72]),
    (20.0, 50.0, "aro", 1.0710001, 1100.0, [0.8, 0.72]),
]


def halves(k_shift=0.0):
    """FULL を EFPD 10 の点で切った前と後の半分。後の半分の EFPD と日数は 0 から"""
    first = core(FULL[:2])
    second = core(
        [
            (e - 10.0, d - 10.0, s, k + k_shift, ppm, b)
            for e, d, s, k, ppm, b in FULL[1:]
        ]
    )
    return first, second


class Restarts(unittest.TestCase):
    def test_second_halves_pair_with_the_uncut_points_after_the_cut(self):
        first, second = halves(k_shift=1.0e-8)
        (cut,) = compare_core_restart.summarize(
            core(FULL), {"a": (first, second)}, None
        )["cuts"]
        self.assertEqual(cut["name"], "a")
        self.assertEqual(cut["cut"], {"efpd": 10.0, "days": 10.0})
        self.assertEqual([p["efpd"] for p in cut["points"]], [10.0, 20.0, 20.0])
        self.assertEqual([p["days"] for p in cut["points"]], [10.0, 20.0, 50.0])
        p = cut["points"][1]
        self.assertAlmostEqual(
            p["rho_pcm"], (1 / 1.071 - 1 / (1.071 + 1e-8)) * 1e5, delta=1e-9
        )
        self.assertEqual(p["boron_diff"], 0.0)
        self.assertEqual(p["burnup_max_diff"], 0.0)
        self.assertEqual(p["density_max_relative"], 0.0)
        self.assertEqual(cut["auxiliary_max_relative"], 0.0)
        # 前の半分は、切らない燃焼の初めの点と同じ
        self.assertEqual(cut["first_k_max_diff"], 0.0)

    def test_density_differences_are_relative_to_the_largest_density_of_each_nuclide(
        self,
    ):
        first, second = halves()
        second["points"][1]["densities"][0][0][1] += 1.0e-12
        second["auxiliary_densities"] = [[[1.1e-9], [1.0e-9]]] * 2
        (cut,) = compare_core_restart.summarize(
            core(FULL), {"a": (first, second)}, None
        )["cuts"]
        self.assertAlmostEqual(cut["points"][1]["density_max_relative"], 1e-12 / 2e-9)
        # 補助の核種も、切らない燃焼の最大の数密度（1e-9）に対する比
        self.assertAlmostEqual(cut["auxiliary_max_relative"], 0.1)

    def test_halves_that_do_not_join_the_uncut_depletion_are_errors(self):
        for edit in [
            lambda f, s: s.update(axial_mesh=[0.0, 5.0, 20.0]),
            lambda f, s: s.update(assemblies=["A01", "A03"]),
            lambda f, s: s.update(heavy_metal=[[4000.0, 3999.0]] * 2),
            lambda f, s: s["points"].pop(),
            lambda f, s: s["points"][1].update(efpd=11.0),
            lambda f, s: f["points"][1].update(burnups=[[0.0, 0.0]] * 2),
            lambda f, s: s.update(xenon="tracked"),
            # 数密度がない、形が違う（比べずに差 0 と書かない）
            lambda f, s: [p.pop("densities") for p in s["points"][1:]],
            lambda f, s: s["points"][2]["densities"][0].pop(),
            lambda f, s: s.update(auxiliary_densities=None),
            lambda f, s: s["auxiliary_densities"][0][0].clear(),
        ]:
            first, second = halves()
            edit(first, second)
            with self.assertRaises(SystemExit):
                compare_core_restart.summarize(core(FULL), {"a": (first, second)}, None)


RELOAD = [
    (0.0, 0.0, "aro", 1.071, 1200.0, [0.0, 0.8, 0.72]),
    (10.0, 10.0, "aro", 1.071, 950.0, [0.45, 1.15, 1.1]),
]


class Reloads(unittest.TestCase):
    def test_points_give_the_k_boron_and_assembly_burnups_by_id(self):
        # 燃焼点 0 の A01 と A02 は前のサイクルの最後の点の状態、新しい A05 は燃焼度 0
        previous = core(FULL)
        reload = core(RELOAD, assemblies=("A01", "A02", "A05"))
        last = previous["points"][-1]
        start = reload["points"][0]
        start["burnups"] = [*last["burnups"], [0.0, 0.0]]
        start["densities"] = [*last["densities"], [[7e-4, 0.0], [7e-4, 0.0]]]
        s = compare_core_restart.summarize(previous, {}, reload)["reload"]
        self.assertEqual(s["assemblies"], ["A01", "A02", "A05"])
        self.assertEqual(s["inherited"], ["A01", "A02"])
        self.assertEqual(s["fresh"], ["A05"])
        self.assertEqual(s["discharged"], [])
        self.assertEqual(s["previous"], {"efpd": 20.0, "days": 50.0})
        self.assertEqual([p["k"] for p in s["points"]], [1.071, 1.071])
        self.assertEqual([p["boron"] for p in s["points"]], [1200.0, 950.0])
        self.assertEqual(s["points"][1]["assembly_burnups"], [0.45, 1.15, 1.1])

    def test_reloads_must_start_from_the_inherited_states_and_fresh_fuel(self):
        previous = core(FULL)
        reload = core(RELOAD, assemblies=("A01", "A02", "A05"))
        # 点 0 の A01 が前のサイクルの最後の点の A01 と違う
        with self.assertRaises(SystemExit):
            compare_core_restart.summarize(previous, {}, reload)

    def test_reloads_of_other_nuclides_masses_or_without_densities_are_errors(self):
        def started():
            previous = core(FULL)
            reload = core(RELOAD, assemblies=("A01", "A02", "A05"))
            last, start = previous["points"][-1], reload["points"][0]
            start["burnups"] = [*last["burnups"], [0.0, 0.0]]
            start["densities"] = [*last["densities"], [[7e-4, 0.0], [7e-4, 0.0]]]
            return previous, reload

        # 直さない組は止まらない
        previous, reload = started()
        compare_core_restart.summarize(previous, {}, reload)
        for edit in [
            lambda r: r["points"][0].pop("densities"),
            lambda r: r.update(nuclides=["U-235", "Xe-136"]),
            lambda r: r.update(auxiliary=["Te-135"]),
            lambda r: r["heavy_metal"].__setitem__(0, [4000.0, 3999.0]),
        ]:
            previous, reload = started()
            edit(reload)
            with self.assertRaises(SystemExit):
                compare_core_restart.summarize(previous, {}, reload)

    def test_discharged_assemblies_are_listed(self):
        previous = core(FULL, assemblies=("A01", "A02", "A03"))
        reload = core(
            [(0.0, 0.0, "aro", 1.071, 1200.0, [0.8, 0.0])], assemblies=("A01", "A05")
        )
        reload["points"][0]["burnups"][0] = previous["points"][-1]["burnups"][0]
        reload["points"][0]["densities"][0] = previous["points"][-1]["densities"][0]
        reload["points"][0]["densities"][1] = [[7e-4, 0.0], [7e-4, 0.0]]
        s = compare_core_restart.summarize(previous, {}, reload)["reload"]
        self.assertEqual(s["discharged"], ["A02", "A03"])
        self.assertEqual(s["fresh"], ["A05"])


class Report(unittest.TestCase):
    def test_reports_give_the_largest_differences_and_the_reload_points(self):
        first, second = halves(k_shift=1.0e-8)
        summary = {
            "inputs": [{"path": "a", "md5": "0" * 32}],
            **compare_core_restart.summarize(core(FULL), {"a": (first, second)}, None),
        }
        text = "\n".join(compare_core_restart.report(summary))
        self.assertIn("切る点 a", text)
        self.assertIn("Δρ", text)


if __name__ == "__main__":
    unittest.main()
