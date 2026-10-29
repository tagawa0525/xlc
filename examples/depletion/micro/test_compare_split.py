"""compare_split.py（ミクロ燃焼の線形の補正の誤差を、自己遮蔽、空間の分布、少数群の重みに分ける）。

実行: python3 -m unittest discover -s examples/depletion/micro -p 'test_compare_split.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_linearization as cl
import compare_split as cs

TOTAL = [0.5, 1.2]
SCATTER = [[0.45, 0.02], [0.0, 1.1]]
NU_FISSION = [0.006, 0.14]
CHI = [1.0, 0.0]
# 追う核種 1 つの微視的断面積 (barn)
SIGMA = {
    "total": [10.0, 400.0],
    "scatter": [[9.0, 0.0], [0.0, 10.0]],
    "nu_fission": [1.0, 700.0],
}
DENSITIES = [[1e-4], [0.9e-4]]


def constants(total, scatter, nu_fission, chi=CHI):
    out = {"total": total, "scatter": scatter, "nu_fission": nu_fission, "chi": chi}
    out["k_inf"] = cl.k_inf(out)
    return out


def sigma_plus(sigma, dt, dnf):
    """熱群の total と nu_fission を足した微視的断面積（散乱はそのまま）"""
    return {
        "total": [sigma["total"][0], sigma["total"][1] + dt],
        "scatter": sigma["scatter"],
        "nu_fission": [sigma["nu_fission"][0], sigma["nu_fission"][1] + dnf],
    }


def with_nuclide(state, density, sigma_b, sigma_a, chi=None):
    """状態に追う核種の巨視的断面積の差 N (σ_b − σ_a) を足し、k∞ を求め直したもの"""
    out = cs.add(state, cs.nuclide_delta(density, sigma_b, sigma_a))
    if chi is not None:
        out["chi"] = chi
    return constants(out["total"], out["scatter"], out["nu_fission"], out["chi"])


def table(point1, sigma1, densities=DENSITIES):
    """燃焼度 2 点、状態は点ごとに 1 つ。点 0 は公称の定数"""
    base = constants(TOTAL, SCATTER, NU_FISSION)
    return {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "points": [[0], [1]],
        "states": [base, point1],
        "micro": {
            "nuclides": ["U-235"],
            "densities": densities,
            "states": [[SIGMA], [sigma1]],
        },
    }


def tables():
    """少数群の重み、空間の分布、自己遮蔽の分を、それぞれ追う核種の熱群の σ の差として持つ 5 つのテーブル。少数群の重みは χ も変える"""
    base = constants(TOTAL, SCATTER, NU_FISSION)
    n1 = DENSITIES[1][0]
    predicted = cl.corrected(base, [SIGMA], [n1 - DENSITIES[0][0]])
    s_spread = sigma_plus(SIGMA, -3.0, -5.0)
    spread = with_nuclide(predicted, n1, s_spread, SIGMA, chi=[0.98, 0.02])
    s_fixed = sigma_plus(s_spread, 1.0, 2.0)
    fixed = with_nuclide(spread, n1, s_fixed, s_spread)
    s_shifted = sigma_plus(s_fixed, -6.0, -4.0)
    shifted = with_nuclide(fixed, n1, s_shifted, s_fixed)
    s_live = sigma_plus(s_shifted, 0.5, 0.5)
    live = with_nuclide(shifted, n1, s_live, s_shifted)
    return {
        "nominal": table(constants(TOTAL, SCATTER, NU_FISSION), SIGMA),
        "shifted": table(shifted, s_shifted),
        "fixed": table(fixed, s_fixed),
        "spread": table(spread, s_spread),
        "spread_live": table(live, s_live),
    }


class Split(unittest.TestCase):
    def test_the_error_is_the_sum_of_shielding_spatial_and_weight(self):
        (row,) = cs.compare(tables())["rows"]
        self.assertAlmostEqual(
            row["error_pcm"],
            row["shielding_pcm"] + row["spatial_pcm"] + row["weight_pcm"],
            delta=1e-9,
        )
        # 順を入れ替えても、2 つの候補の和は同じ
        self.assertAlmostEqual(
            row["spatial_live_pcm"] + row["shielding_spread_pcm"],
            row["spatial_pcm"] + row["shielding_pcm"],
            delta=1e-9,
        )
        for key in ("shielding_pcm", "spatial_pcm", "weight_pcm"):
            self.assertGreater(abs(row[key]), 1.0, key)

    def test_a_change_of_one_tracked_cross_section_falls_on_that_nuclide(self):
        (row,) = cs.compare(tables())["rows"]
        for part, key in (("shielding", "shielding_pcm"), ("spatial", "spatial_pcm")):
            s = row["splits"][part]
            # 2 群の k と k∞ は同じ定数から求めたので、輸送の差と内訳の全体は同じ
            self.assertAlmostEqual(s["whole_pcm"], row[key], delta=1e-9)
            self.assertAlmostEqual(s["nuclides"]["U-235"], s["whole_pcm"], delta=1e-9)
            for rest in ("others_pcm", "chi_pcm", "cross_pcm", "other_tracked_pcm"):
                self.assertAlmostEqual(s[rest], 0.0, delta=1e-9, msg=f"{part} {rest}")

    def test_the_weight_splits_into_the_nuclide_chi_and_their_cross_term(self):
        (row,) = cs.compare(tables())["rows"]
        s = row["splits"]["weight"]
        self.assertAlmostEqual(s["whole_pcm"], row["weight_pcm"], delta=1e-9)
        self.assertGreater(abs(s["chi_pcm"]), 1.0)
        self.assertAlmostEqual(s["others_pcm"], 0.0, delta=1e-9)
        self.assertAlmostEqual(
            s["nuclides"]["U-235"] + s["chi_pcm"] + s["cross_pcm"],
            s["whole_pcm"],
            delta=1e-9,
        )

    def test_a_change_of_the_other_cross_sections_falls_on_the_others(self):
        t = tables()
        fixed = t["fixed"]["states"][1]
        # 自己遮蔽の分を、追わない核種の熱群の吸収の差にする（追う核種の σ は fixed と同じ）
        shifted = constants(
            [fixed["total"][0], fixed["total"][1] + 2e-3],
            fixed["scatter"],
            fixed["nu_fission"],
            fixed["chi"],
        )
        t["shifted"]["states"][1] = shifted
        t["shifted"]["micro"]["states"][1] = t["fixed"]["micro"]["states"][1]
        (row,) = cs.compare(t)["rows"]
        s = row["splits"]["shielding"]
        self.assertAlmostEqual(s["others_pcm"], row["shielding_pcm"], delta=1e-9)
        self.assertAlmostEqual(s["nuclides"]["U-235"], 0.0, delta=1e-9)

    def test_small_tracked_contributions_are_summed_as_the_other_tracked(self):
        t = tables()
        summary = cs.compare(t)
        self.assertEqual(summary["shown_nuclides"], ["U-235"])
        original = cs.NUCLIDE_THRESHOLD_PCM
        try:
            cs.NUCLIDE_THRESHOLD_PCM = 1e9
            summary = cs.compare(t)
        finally:
            cs.NUCLIDE_THRESHOLD_PCM = original
        self.assertEqual(summary["shown_nuclides"], [])
        s = summary["rows"][0]["splits"]["shielding"]
        self.assertAlmostEqual(s["other_tracked_pcm"], s["whole_pcm"], delta=1e-9)


class Checks(unittest.TestCase):
    def test_the_replaced_compositions_must_keep_the_tracked_densities(self):
        t = tables()
        t["spread"]["micro"]["densities"] = [[1e-4], [0.9e-4 * (1 + 1e-9)]]
        with self.assertRaises(SystemExit):
            cs.compare(t)

    def test_the_fresh_points_must_agree(self):
        t = tables()
        fresh = copy.deepcopy(t["fixed"]["states"][0])
        fresh["k_inf"] *= 1 + 1e-9
        t["fixed"]["states"][0] = fresh
        with self.assertRaises(SystemExit):
            cs.compare(t)

    def test_every_table_needs_micro_data_of_the_same_nuclides(self):
        t = tables()
        del t["spread_live"]["micro"]
        with self.assertRaises(SystemExit):
            cs.compare(t)
        t = tables()
        t["fixed"]["micro"]["nuclides"] = ["Pu-239"]
        with self.assertRaises(SystemExit):
            cs.compare(t)


class Report(unittest.TestCase):
    def test_the_report_names_each_candidate(self):
        summary = cs.compare(tables())
        # 5 MWd/kg より前と後の両方の組がある要約にする
        late = copy.deepcopy(summary["rows"][0])
        late["from"], late["to"] = 5.0, 6.0
        summary["rows"].append(late)
        text = cs.report(summary)
        for word in ("自己遮蔽", "空間の分布", "少数群の重み", "交差"):
            self.assertIn(word, text)


if __name__ == "__main__":
    unittest.main()
