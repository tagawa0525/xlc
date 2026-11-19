"""compare_background.py（同じ組成で、輸送しない核種の扱いだけを変えた状態の k の差を、核種、効果、エネルギーの範囲に分ける）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_background.py'
"""

import copy
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_background

# 3 群: 0 は UFG の範囲より上、1 は UFG の範囲、2 は UFG の範囲より下
BOUNDS = [1e7, 1e6, 10.0, 1e-5]
RANGE = [10.0, 1e6]
TEMPERATURE = 565.0
# 燃焼の領域の体積 / 格子の体積
VOLUME_RATIO = 0.5
# 燃料の領域の中性子束の、格子で均質化した中性子束に対する比（O-16 の実効断面積 / 無限希釈）
FLUX_RATIO = np.array([1.0, 0.9, 0.7])
O16_TOTAL = np.array([2.0, 3.8, 3.9])
CHI = np.array([0.9, 0.1, 0.0])


def transported(n, total, scatter, nu_fission=(0.0, 0.0, 0.0)):
    return {
        "density": n,
        "total": np.array(total, dtype=float),
        "scatter": np.array(scatter, dtype=float),
        "nu_fission": np.array(nu_fission, dtype=float),
    }


def base_nuclides():
    return {
        "U-238": transported(
            0.02,
            [7.0, 20.0, 10.0],
            [[5.0, 1.0, 0.0], [0.0, 9.0, 0.5], [0.0, 0.0, 8.0]],
            [0.6, 0.0, 0.0],
        ),
        "U-235": transported(
            0.0008,
            [6.0, 40.0, 500.0],
            [[4.0, 0.5, 0.0], [0.0, 10.0, 0.2], [0.0, 0.0, 10.0]],
            [3.0, 50.0, 900.0],
        ),
        "O-16": transported(
            0.045,
            O16_TOTAL * FLUX_RATIO,
            np.diag(O16_TOTAL * FLUX_RATIO) * 0.999,
        ),
    }


# 減速材などの燃焼させない物質の、格子で均質化した巨視的断面積
REST_TOTAL = np.array([0.3, 1.0, 1.5])
REST_SCATTER = np.array([[0.1, 0.19, 0.0], [0.0, 0.5, 0.49], [0.0, 0.0, 1.48]])


def build(nuclides, dilute, cooling_scale=1.0, regions=1):
    """`nuclides` を輸送し（名前 → transported）、`dilute`（名前 → (燃焼の領域の数密度, 群ごとの吸収)）を無限希釈で数える
    状態の、ブランチ計算の出力、履歴、ライブラリ"""
    names = list(nuclides)
    dilute_names = list(dilute)
    total = REST_TOTAL.copy()
    scatter = REST_SCATTER.copy()
    nu_fission = np.zeros(3)
    for x in nuclides.values():
        total += x["density"] * x["total"]
        scatter += x["density"] * x["scatter"]
        nu_fission += x["density"] * x["nu_fission"]
    for n, sigma in dilute.values():
        total += n * VOLUME_RATIO * np.array(sigma) * FLUX_RATIO
    removal = np.diag(total) - scatter.T
    phi = np.linalg.solve(removal, CHI)
    k = float(nu_fission @ phi)
    cooling = []
    for n, sigma in dilute.values():
        one = float((np.array(sigma) * FLUX_RATIO) @ phi / phi.sum()) * cooling_scale
        cooling.append(
            {
                "fission": 0.0,
                "thermal_fission": None,
                "reactions": [{"mt": 102, "channels": [one]}],
            }
        )
    branch = {
        "format": "xlc-table",
        "bounds": BOUNDS,
        "axes": [
            {"name": "burnup", "points": [60.0]},
            {"name": "fuel", "points": [TEMPERATURE]},
        ],
        "states": [
            {
                "k_inf": k,
                "total": total.tolist(),
                "scatter": scatter.tolist(),
                "nu_fission": nu_fission.tolist(),
                "chi": CHI.tolist(),
            }
        ],
        "micro": {
            "nuclides": names,
            "densities": [[nuclides[n]["density"] for n in names]],
            "states": [
                [
                    {
                        "total": nuclides[n]["total"].tolist(),
                        "scatter": nuclides[n]["scatter"].tolist(),
                        "nu_fission": nuclides[n]["nu_fission"].tolist(),
                    }
                    for n in names
                ]
            ],
        },
        "cooling": {"nuclides": dilute_names, "points": [cooling]},
    }
    every = sorted(set(names) | set(dilute_names) | {"Xe-999"})
    region = {n: nuclides[n]["density"] / VOLUME_RATIO for n in names}
    region.update({n: d[0] for n, d in dilute.items()})
    history = {
        "format": "xlc-depletion-history",
        "nuclides": every,
        "regions": [{"volume": 1.0}] * regions,
        "transportable": names,
        "dilute_absorbers": dilute_names,
        "points": [
            {
                "burnup": 60.0,
                "k_inf": k,
                "densities": [[region.get(n, 0.0) for n in every]] * regions,
            }
        ],
    }
    library = {
        "group_bounds": BOUNDS,
        "pointwise_range": RANGE,
        "nuclides": [
            {
                "name": "O-16",
                "temperatures": [
                    {
                        "temperature": TEMPERATURE,
                        "shielded": [{"total": O16_TOTAL.tolist()}],
                    }
                ],
            }
        ],
        "activation": [
            {
                "name": n,
                "temperatures": [
                    {
                        "temperature": TEMPERATURE,
                        "reactions": [
                            {"mt": 102, "channels": [{"cross_section": list(sigma)}]},
                            {
                                "mt": 16,
                                "channels": [{"cross_section": [9.0, 0.0, 0.0]}],
                            },
                        ],
                    }
                ],
            }
            for n, (_, sigma) in dilute.items()
        ],
    }
    return branch, history, library


# 輸送する形と無限希釈で数える形を比べる核種。散乱は群の中だけ（減速しない）
FP_SCATTER_IN_GROUP = [[5.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 4.0]]
FP_ABSORPTION = [1.0, 25.0, 8.0]
# 輸送したときの吸収（全断面積 − 散乱）が、無限希釈の吸収を燃料の領域の中性子束で均質化したものと同じになる全断面積
FP_TOTAL = np.sum(FP_SCATTER_IN_GROUP, axis=1) + np.array(FP_ABSORPTION) * FLUX_RATIO


def case(nuclides, dilute, **kw):
    return compare_background.case(*build(nuclides, dilute, **kw), 60.0)


def rho(a, b):
    return (1.0 / a["k_matrix"] - 1.0 / b["k_matrix"]) * 1e5


class Step(unittest.TestCase):
    def test_parts_sum_to_the_reactivity_difference(self):
        reference = base_nuclides()
        reference["Fp"] = transported(
            0.002, FP_TOTAL, [[4.0, 0.8, 0.0], [0.0, 9.0, 0.4], [0.0, 0.0, 4.0]]
        )
        form = base_nuclides()
        form["U-238"]["total"] = form["U-238"]["total"] - np.array([0.0, 2.0, 0.0])
        a = case(reference, {"Ab": (0.001, [0.5, 3.0, 2.0])})
        b = case(form, {"Ab": (0.001, [0.5, 3.0, 2.0]), "Fp": (0.004, FP_ABSORPTION)})
        step = compare_background.step(a, b)
        self.assertAlmostEqual(step["total_pcm"], rho(a, b), places=9)
        self.assertAlmostEqual(sum(step["parts"].values()), step["total_pcm"], places=9)
        self.assertLess(abs(step["closure_pcm"]), 1e-9)

    def test_a_change_of_one_resonance_cross_section_falls_in_its_own_part(self):
        form = base_nuclides()
        form["U-238"]["total"] = form["U-238"]["total"] - np.array([0.0, 2.0, 0.0])
        step = compare_background.step(case(base_nuclides(), {}), case(form, {}))
        self.assertGreater(step["total_pcm"], 0.0)
        self.assertAlmostEqual(
            step["parts"]["u238_resonance"], step["total_pcm"], places=9
        )
        others = {k: v for k, v in step["parts"].items() if k != "u238_resonance"}
        for name, value in others.items():
            self.assertAlmostEqual(value, 0.0, places=9, msg=name)

    def test_in_group_scattering_of_a_moved_nuclide_does_not_count(self):
        reference = base_nuclides()
        reference["Fp"] = transported(0.002, FP_TOTAL, FP_SCATTER_IN_GROUP)
        a = case(reference, {})
        b = case(base_nuclides(), {"Fp": (0.004, FP_ABSORPTION)})
        step = compare_background.step(a, b)
        self.assertAlmostEqual(step["parts"]["moved_scattering_fast"], 0.0, places=9)
        self.assertAlmostEqual(step["parts"]["moved_scattering_below"], 0.0, places=9)
        # 輸送の吸収（全断面積 − 散乱）と同じ無限希釈の吸収なので、吸収の分も 0
        self.assertAlmostEqual(step["parts"]["moved_absorption"], 0.0, places=9)
        fp = step["nuclides"]["Fp"]
        self.assertEqual(fp["kind"], "to_dilute")

    def test_dropping_the_down_scattering_of_a_moved_nuclide_counts_in_its_range(self):
        reference = base_nuclides()
        # 群 0 から 1 への減速を持つ（群の中の散乱を 1 減らし、その分を移す）
        reference["Fp"] = transported(
            0.002, FP_TOTAL, [[4.0, 1.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 4.0]]
        )
        a = case(reference, {})
        b = case(base_nuclides(), {"Fp": (0.004, FP_ABSORPTION)})
        step = compare_background.step(a, b)
        self.assertNotAlmostEqual(step["parts"]["moved_scattering_fast"], 0.0, places=3)
        self.assertAlmostEqual(step["parts"]["moved_scattering_below"], 0.0, places=9)
        self.assertAlmostEqual(step["parts"]["moved_absorption"], 0.0, places=9)
        ranges = step["nuclides"]["Fp"]["scattering"]
        self.assertAlmostEqual(
            ranges["fast"], step["parts"]["moved_scattering_fast"], places=9
        )

    def test_a_nuclide_counted_in_both_states_is_rest(self):
        a = case(base_nuclides(), {"Ab": (0.001, [0.5, 3.0, 2.0])})
        b = case(base_nuclides(), {"Ab": (0.002, [0.5, 3.0, 2.0])})
        step = compare_background.step(a, b)
        self.assertLess(step["total_pcm"], 0.0)
        self.assertAlmostEqual(step["parts"]["rest"], step["total_pcm"], places=9)
        self.assertEqual(step["nuclides"]["Ab"]["kind"], "dilute")


class Case(unittest.TestCase):
    def test_the_dilute_absorption_matches_the_cooling_cross_sections(self):
        c = case(base_nuclides(), {"Fp": (0.004, FP_ABSORPTION)})
        self.assertLess(abs(c["dilute_check"]), 1e-12)

    def test_a_dilute_absorption_off_the_cooling_cross_sections_is_an_error(self):
        with self.assertRaisesRegex(SystemExit, "cooling"):
            case(base_nuclides(), {"Fp": (0.004, FP_ABSORPTION)}, cooling_scale=1.01)

    def test_micro_cross_sections_of_every_transported_nuclide_are_needed(self):
        branch, history, library = build(base_nuclides(), {})
        history["transportable"] = history["transportable"] + ["Xe-999"]
        with self.assertRaisesRegex(SystemExit, "micro"):
            compare_background.case(branch, history, library, 60.0)

    def test_several_depletion_regions_are_an_error(self):
        with self.assertRaisesRegex(SystemExit, "one depletion region"):
            compare_background.case(*build(base_nuclides(), {}, regions=2), 60.0)

    def test_few_groups_other_than_the_library_groups_are_an_error(self):
        branch, history, library = build(base_nuclides(), {})
        library = copy.deepcopy(library)
        library["group_bounds"] = [1e7, 1e5, 10.0, 1e-5]
        with self.assertRaisesRegex(SystemExit, "library groups"):
            compare_background.case(branch, history, library, 60.0)


class Report(unittest.TestCase):
    def test_report_has_a_line_per_step(self):
        a = case(base_nuclides(), {})
        form = base_nuclides()
        form["U-238"]["total"] = form["U-238"]["total"] - np.array([0.0, 2.0, 0.0])
        b = case(form, {})
        summary = compare_background.summarize([("ref", a), ("form", b)])
        lines = compare_background.report({**summary, "inputs": []})
        self.assertTrue(any(line.startswith("ref → form") for line in lines))


if __name__ == "__main__":
    unittest.main()
