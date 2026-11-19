"""compare_chain_versions.py の要約（2 つの版の燃焼チェーンの、量ごとの違い）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_chain_versions.py'
"""

import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_chain_versions


def nuclide(name, z, a, half_life=None, modes=(), energy=(0.0, 0.0, 0.0)):
    """名前で生成物を書いた核種（`build` が添字に直す）。`modes` は (型, 分岐比, [(生成物, 個数)])"""
    decay = None
    if half_life is not None:
        decay = {
            "half_life": half_life,
            "energy": dict(zip(["light", "electromagnetic", "heavy"], energy)),
            "modes": [
                {
                    "mode": m,
                    "ratio": r,
                    "products": [{"nuclide": p, "count": c} for p, c in ps],
                }
                for m, r, ps in modes
            ],
        }
    return {
        "name": name,
        "z": z,
        "a": a,
        "state": 0,
        "excitation": 0.0,
        "decay": decay,
        "reactions": [],
        "fission": None,
    }


def build(nuclides):
    """生成物の名前を添字にした、bin/chain の出力と同じ形のチェーン"""
    index = {n["name"]: i for i, n in enumerate(nuclides)}
    nuclides = copy.deepcopy(nuclides)
    for n in nuclides:
        if n["decay"] is not None:
            for m in n["decay"]["modes"]:
                for p in m["products"]:
                    p["nuclide"] = index[p["nuclide"]]
        for r in n["reactions"]:
            for c in r["channels"]:
                for p in c["products"]:
                    p["nuclide"] = index[p["nuclide"]]
        if n["fission"] is not None:
            for table in n["fission"]["yields"]:
                for y in table:
                    y["nuclide"] = index[y["nuclide"]]
    return {"format": "xlc-chain", "nuclides": nuclides}


def chain_nuclides():
    """U-235 が核分裂し、I-135 → Xe-135 → Cs-135 と崩壊するチェーン（名前で書く）"""
    u235 = nuclide(
        "U-235",
        92,
        235,
        2.2e16,
        [("alpha", 1.0, [("Th-231", 1), ("He-4", 1)])],
        (4.0e4, 1.6e5, 4.4e6),
    )
    u235["reactions"] = [
        {
            "mt": 102,
            "channels": [
                {"final_state": None, "products": [{"nuclide": "U-236", "count": 1}]}
            ],
        }
    ]
    u235["fission"] = {
        "q": 1.934e8,
        "yield_source": "U-235",
        "energies": [0.0253, 5.0e5],
        "yields": [
            [
                {"nuclide": "I-135", "value": 0.0293},
                {"nuclide": "Xe-135", "value": 0.0024},
            ],
            [
                {"nuclide": "I-135", "value": 0.0288},
                {"nuclide": "Xe-135", "value": 0.0030},
            ],
        ],
    }
    return [
        nuclide("He-4", 2, 4),
        nuclide("Th-231", 90, 231, 9.2e4, [("beta-", 1.0, [("Pa-231", 1)])]),
        nuclide("Pa-231", 91, 231),
        u235,
        nuclide("U-236", 92, 236),
        nuclide(
            "I-135",
            53,
            135,
            23652.0,
            [("beta-", 1.0, [("Xe-135", 1)])],
            (4.8e5, 1.6e6, 0.0),
        ),
        nuclide(
            "Xe-135",
            54,
            135,
            32904.0,
            [("beta-", 1.0, [("Cs-135", 1)])],
            (3.2e5, 2.5e5, 0.0),
        ),
        nuclide("Cs-135", 55, 135),
    ]


def by_name(nuclides, name):
    return next(n for n in nuclides if n["name"] == name)


def counts(summary):
    return {q: s["differ"] for q, s in summary["quantities"].items()}


class SameValues(unittest.TestCase):
    def test_a_chain_has_no_differences_from_itself(self):
        a = build(chain_nuclides())
        summary = compare_chain_versions.summarize(a, a)
        self.assertEqual(set(counts(summary).values()), {0})
        self.assertEqual(summary["nuclides"]["only_a"], [])
        self.assertEqual(summary["nuclides"]["only_b"], [])

    def test_the_order_of_nuclides_modes_and_yields_is_not_a_difference(self):
        a = build(chain_nuclides())
        reordered = list(reversed(chain_nuclides()))
        for table in by_name(reordered, "U-235")["fission"]["yields"]:
            table.reverse()
        u235 = by_name(reordered, "U-235")
        u235["decay"]["modes"][0]["products"].reverse()
        b = build(reordered)
        self.assertEqual(
            set(counts(compare_chain_versions.summarize(a, b)).values()), {0}
        )

    def test_a_branch_written_in_two_parts_is_its_sum(self):
        # 同じ型と娘を 2 度書いた評価がある（Es-246 の ec/beta+）
        nuclides = chain_nuclides()
        by_name(nuclides, "I-135")["decay"]["modes"] = [
            {
                "mode": "beta-",
                "ratio": 0.25,
                "products": [{"nuclide": "Xe-135", "count": 1}],
            },
            {
                "mode": "beta-",
                "ratio": 0.75,
                "products": [{"nuclide": "Xe-135", "count": 1}],
            },
        ]
        summary = compare_chain_versions.summarize(
            build(chain_nuclides()), build(nuclides)
        )
        self.assertEqual(summary["quantities"]["decay_modes"]["differ"], 0)

    def test_a_zero_yield_is_the_same_as_no_yield(self):
        nuclides = chain_nuclides()
        by_name(nuclides, "U-235")["fission"]["yields"][0].append(
            {"nuclide": "Cs-135", "value": 0.0}
        )
        summary = compare_chain_versions.summarize(
            build(chain_nuclides()), build(nuclides)
        )
        self.assertEqual(summary["quantities"]["fission_yields"]["differ"], 0)


class Differences(unittest.TestCase):
    def setUp(self):
        self.a = build(chain_nuclides())
        nuclides = chain_nuclides()
        by_name(nuclides, "Xe-135")["decay"]["half_life"] = 32940.0
        by_name(nuclides, "Xe-135")["decay"]["energy"]["light"] = 3.3e5
        by_name(nuclides, "I-135")["decay"]["modes"] = [
            {
                "mode": "beta-",
                "ratio": 0.9,
                "products": [{"nuclide": "Xe-135", "count": 1}],
            },
            {"mode": "beta-,n", "ratio": 0.1, "products": []},
        ]
        u235 = by_name(nuclides, "U-235")
        u235["fission"]["q"] = 1.935e8
        u235["fission"]["yields"][0][0]["value"] = 0.0300
        u235["reactions"][0]["channels"][0]["products"] = []
        self.b = build(nuclides)
        self.summary = compare_chain_versions.summarize(self.a, self.b)

    def test_each_quantity_counts_its_differing_items(self):
        self.assertEqual(
            counts(self.summary),
            {
                "excitation": 0,
                "half_life": 1,
                "decay_energy": 1,
                # I-135 の beta- の分岐比と、b だけの beta-,n
                "decay_modes": 2,
                "fission_q": 1,
                "yield_sources": 0,
                "yield_energies": 0,
                "fission_yields": 1,
                "reactions": 1,
            },
        )

    def test_differences_list_both_values_and_the_relative_difference(self):
        (item,) = self.summary["quantities"]["half_life"]["largest"]
        self.assertEqual(item["nuclide"], "Xe-135")
        self.assertEqual((item["a"], item["b"]), (32904.0, 32940.0))
        self.assertAlmostEqual(item["relative"], 36.0 / 32940.0)
        (y,) = self.summary["quantities"]["fission_yields"]["largest"]
        self.assertEqual(
            (y["nuclide"], y["energy"], y["product"]), ("U-235", 0.0253, "I-135")
        )

    def test_the_largest_differences_come_first(self):
        modes = self.summary["quantities"]["decay_modes"]["largest"]
        # b だけの分岐（a は 0）の相対差 1 が、0.1 の差より先
        self.assertEqual([m["mode"] for m in modes], ["beta-,n", "beta-"])
        self.assertEqual(modes[0]["a"], 0.0)
        self.assertEqual(modes[0]["relative"], 1.0)

    def test_yields_are_summed_per_parent_and_energy(self):
        table = self.summary["yield_tables"]
        row = next(
            r for r in table if r["nuclide"] == "U-235" and r["energy"] == 0.0253
        )
        self.assertEqual(row["differ"], 1)
        self.assertAlmostEqual(row["sum_abs"], 0.0007)

    def test_the_watched_nuclides_show_both_values_even_when_equal(self):
        watched = self.summary["watched"]
        xe = next(w for w in watched["nuclides"] if w["nuclide"] == "Xe-135")
        self.assertEqual(xe["half_life"], {"a": 32904.0, "b": 32940.0})
        yields = [
            y
            for y in watched["yields"]
            if y["nuclide"] == "U-235"
            and y["product"] == "Xe-135"
            and y["energy"] == 0.0253
        ]
        self.assertEqual([(y["a"], y["b"]) for y in yields], [(0.0024, 0.0024)])

    def test_the_watched_parents_show_their_decay_and_fission_q(self):
        watched = self.summary["watched"]
        u235 = next(w for w in watched["nuclides"] if w["nuclide"] == "U-235")
        self.assertEqual(u235["half_life"], {"a": 2.2e16, "b": 2.2e16})
        self.assertEqual(u235["fission_q"], {"a": 1.934e8, "b": 1.935e8})
        xe = next(w for w in watched["nuclides"] if w["nuclide"] == "Xe-135")
        self.assertIsNone(xe["fission_q"])

    def test_a_watched_nuclide_stable_in_both_has_no_half_life(self):
        nuclides = chain_nuclides() + [nuclide("Nd-143", 60, 143)]
        summary = compare_chain_versions.summarize(build(nuclides), build(nuclides))
        nd = next(w for w in summary["watched"]["nuclides"] if w["nuclide"] == "Nd-143")
        self.assertEqual(nd["half_life"], {"a": None, "b": None})
        self.assertEqual(nd["modes"], [])

    def test_a_stable_and_an_unstable_state_differ_in_half_life(self):
        nuclides = chain_nuclides()
        by_name(nuclides, "Cs-135")["decay"] = {
            "half_life": 7.5e13,
            "energy": {"light": 0.0, "electromagnetic": 0.0, "heavy": 0.0},
            "modes": [],
        }
        summary = compare_chain_versions.summarize(self.a, build(nuclides))
        (item,) = summary["quantities"]["half_life"]["largest"]
        self.assertEqual(
            (item["nuclide"], item["a"], item["b"]), ("Cs-135", None, 7.5e13)
        )
        self.assertEqual(summary["quantities"]["half_life"]["one_sided"], 1)

    def test_a_reaction_in_only_one_chain_is_counted_apart(self):
        # a の U-235 の (n,γ) を b は持たない。経路の違う反応（上の setUp の U-235）とは分けて数える
        nuclides = chain_nuclides()
        by_name(nuclides, "U-235")["reactions"] = []
        summary = compare_chain_versions.summarize(self.a, build(nuclides))
        reactions = summary["quantities"]["reactions"]
        self.assertEqual((reactions["differ"], reactions["one_sided"]), (1, 1))
        self.assertEqual(self.summary["quantities"]["reactions"]["one_sided"], 0)


class NuclideSets(unittest.TestCase):
    def test_nuclides_in_only_one_chain_are_listed(self):
        nuclides = [n for n in chain_nuclides() if n["name"] != "Th-231"]
        by_name(nuclides, "U-235")["decay"]["modes"] = [
            {
                "mode": "alpha",
                "ratio": 1.0,
                "products": [{"nuclide": "He-4", "count": 1}],
            }
        ]
        nuclides.append(nuclide("Cs-136", 55, 136))
        summary = compare_chain_versions.summarize(
            build(chain_nuclides()), build(nuclides)
        )
        self.assertEqual(summary["nuclides"]["only_a"], ["Th-231"])
        self.assertEqual(summary["nuclides"]["only_b"], ["Cs-136"])
        self.assertEqual(summary["nuclides"]["common"], 7)


class Report(unittest.TestCase):
    def test_the_report_gives_the_count_of_each_quantity(self):
        a = build(chain_nuclides())
        text = "\n".join(
            compare_chain_versions.report(
                {"inputs": [], **compare_chain_versions.summarize(a, a)}
            )
        )
        self.assertIn("half_life", text)
        self.assertIn("fission_yields", text)


if __name__ == "__main__":
    unittest.main()
