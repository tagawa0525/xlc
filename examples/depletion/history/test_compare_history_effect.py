"""compare_history_effect.py（履歴の効きと、ミクロ燃焼で追う核種と追わない核種の分）。

実行: python3 -m unittest discover -s examples/depletion/history -p 'test_compare_history_effect.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_history_effect as c

NUCLIDES = list(c.NUCLIDES)


def history(ks: dict, scale: float) -> dict:
    """体積 1 と 3 の 2 領域。核種 i の数密度は領域 r で scale × (i + 1) × (r + 1)"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": NUCLIDES,
        "regions": [
            {"material": "fuel", "cell": 1, "regions": [0], "volume": 1.0},
            {"material": "fuel", "cell": 2, "regions": [1], "volume": 3.0},
        ],
        "points": [
            {
                "burnup": u,
                "k_inf": k,
                "densities": [
                    [scale * (i + 1) * (r + 1) for i in range(len(NUCLIDES))]
                    for r in range(2)
                ],
            }
            for u, k in ks.items()
        ],
    }


def table(ks: dict) -> dict:
    return {
        "format": "xlc-table",
        "axes": [
            {"name": "burnup", "points": list(ks)},
            {"name": "boron", "points": [1300.0]},
            {"name": "moderator", "points": [565.0]},
            {"name": "fuel", "points": [565.0]},
        ],
        "states": [{"k_inf": k} for k in ks.values()],
    }


def pcm(k_from, k_to):
    return (1 / k_from - 1 / k_to) * 1e5


class Compare(unittest.TestCase):
    def rows(self):
        nominal = history({0.0: 1.2, 10.0: 1.1}, 1.0)
        case = {
            "history": history({0.0: 1.0, 10.0: 0.95}, 1.1),
            "branch": table({0.0: 1.0, 10.0: 0.94}),
            "hybrid": table({0.0: 1.0, 10.0: 0.948}),
            "reverse": table({0.0: 1.2, 10.0: 1.105}),
            "reverse_hybrid": table({0.0: 1.2, 10.0: 1.104}),
        }
        return c.compare(nominal, {"density": case})["density"]

    def test_the_effect_is_the_history_against_the_instantaneous_branch(self):
        row = self.rows()[1]
        self.assertAlmostEqual(row["effect_pcm"], pcm(0.94, 0.95), places=9)
        self.assertAlmostEqual(row["between_histories_pcm"], pcm(1.1, 0.95), places=9)

    def test_untracked_and_tracked_parts_sum_to_the_effect(self):
        row = self.rows()[1]
        self.assertAlmostEqual(row["untracked_pcm"], pcm(0.948, 0.95), places=9)
        self.assertAlmostEqual(
            row["tracked_pcm"] + row["untracked_pcm"], row["effect_pcm"], places=9
        )

    def test_the_reverse_is_the_varied_composition_at_the_nominal_conditions(self):
        row = self.rows()[1]
        self.assertAlmostEqual(row["reverse_pcm"], pcm(1.1, 1.105), places=9)

    def test_the_reverse_splits_into_tracked_and_untracked_parts_alike(self):
        row = self.rows()[1]
        self.assertAlmostEqual(
            row["reverse_untracked_pcm"], pcm(1.104, 1.105), places=9
        )
        self.assertAlmostEqual(
            row["reverse_tracked_pcm"] + row["reverse_untracked_pcm"],
            row["reverse_pcm"],
            places=9,
        )

    def test_fresh_points_have_no_effect(self):
        row = self.rows()[0]
        self.assertEqual(row["effect_pcm"], 0.0)
        self.assertEqual(row["reverse_pcm"], 0.0)

    def test_densities_are_volume_weighted_relative_differences(self):
        # 体積で重みを付けた平均の比は scale の比（1.1）
        for value in self.rows()[1]["densities_percent"].values():
            self.assertAlmostEqual(value, 10.0, places=9)

    def test_tables_solving_burnups_the_histories_lack_are_errors(self):
        nominal = history({0.0: 1.2, 10.0: 1.1}, 1.0)
        case = {
            "history": history({0.0: 1.0, 10.0: 0.95}, 1.1),
            "branch": table({0.0: 1.0, 12.0: 0.94}),
            "hybrid": table({0.0: 1.0, 12.0: 0.948}),
            "reverse": table({0.0: 1.2, 12.0: 1.105}),
            "reverse_hybrid": table({0.0: 1.2, 12.0: 1.104}),
        }
        with self.assertRaises(SystemExit):
            c.compare(nominal, {"density": case})

    @staticmethod
    def fresh_case(**fresh) -> dict:
        """燃焼度 0 の k だけを fresh で変えた条件"""
        ks = {
            "branch": {0.0: 1.0, 10.0: 0.94},
            "hybrid": {0.0: 1.0, 10.0: 0.948},
            "reverse": {0.0: 1.2, 10.0: 1.105},
            "reverse_hybrid": {0.0: 1.2, 10.0: 1.104},
        }
        for key, k in fresh.items():
            ks[key][0.0] = k
        return {
            "history": history({0.0: 1.0, 10.0: 0.95}, 1.1),
            **{key: table(k) for key, k in ks.items()},
        }

    def test_fresh_tables_must_reproduce_their_history(self):
        # 燃焼度 0 の組成は 3 本の履歴で同じなので、別の条件で解いたテーブルは燃焼度 0 の k で分かる
        nominal = history({0.0: 1.2, 10.0: 1.1}, 1.0)
        for key, wrong in (
            ("branch", 1.2),
            ("hybrid", 1.2),
            ("reverse", 1.0),
            ("reverse_hybrid", 1.0),
        ):
            with (
                self.subTest(key),
                self.assertRaisesRegex(SystemExit, f"the {key} table misses"),
            ):
                c.compare(nominal, {"density": self.fresh_case(**{key: wrong})})

    def test_fresh_misses_are_the_table_against_its_history(self):
        nominal = history({0.0: 1.2, 10.0: 1.1}, 1.0)
        case = self.fresh_case(branch=1.0000005, reverse=1.1999994)
        fresh = c.compare(nominal, {"density": case})["density"][0]["fresh_pcm"]
        self.assertAlmostEqual(fresh["branch"], pcm(1.0, 1.0000005), places=9)
        self.assertGreater(fresh["branch"], 0)
        self.assertAlmostEqual(fresh["reverse"], pcm(1.2, 1.1999994), places=9)
        self.assertEqual(fresh["hybrid"], 0.0)

    def test_tables_varying_another_axis_are_errors(self):
        t = table({0.0: 1.0})
        t["axes"][1]["points"] = [0.0, 1300.0]
        with self.assertRaises(SystemExit):
            c.table_k(t)

    def test_tables_solving_other_burnups_are_errors(self):
        nominal = history({0.0: 1.2, 10.0: 1.1}, 1.0)
        case = {
            "history": history({0.0: 1.0, 10.0: 0.95}, 1.1),
            "branch": table({0.0: 1.0, 10.0: 0.94}),
            "hybrid": table({0.0: 1.0}),
            "reverse": table({0.0: 1.2, 10.0: 1.105}),
            "reverse_hybrid": table({0.0: 1.2, 10.0: 1.104}),
        }
        with self.assertRaises(SystemExit):
            c.compare(nominal, {"density": case})


if __name__ == "__main__":
    unittest.main()
