"""compare_branch_cross_terms.py の要約（分岐のテーブルを公称の状態からの 1 軸ずつの分岐と、選んだ交差の項の部分格子から組み立てた
ときの、直積の状態との差。ADR-0200）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_branch_cross_terms.py'
"""

import itertools
import json
import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_branch_cross_terms as cct

BORON = [900.0, 1300.0, 1700.0]
MODERATOR = [550.0, 565.0, 580.0]
FUEL = [550.0, 565.0, 580.0]
NOMINAL = {"boron": 1300.0, "moderator": 565.0, "fuel": 565.0}


def state(nu_fission2: float, absorption2: float = 0.08) -> dict:
    """2 群の状態。上方散乱なし、k∞ = (νΣf1 + νΣf2 Σs12 / Σr2) / Σr1（Σr1 = 0.03、Σs12 = 0.02、Σr2 = Σa2 + 0.02）"""
    return {
        "k_inf": (0.005 + nu_fission2 * 0.02 / (absorption2 + 0.02)) / 0.03,
        "total": [0.5, 0.92 + absorption2],
        "absorption": [0.01, absorption2],
        "nu_fission": [0.005, nu_fission2],
        "fission": [0.002, nu_fission2 / 2.4],
        "chi": [1.0, 0.0],
        "diffusion": [1.4, 0.4],
        "scatter": [[0.47, 0.02], [0.0, 0.9]],
        "adf": [[1.0, 1.0]] * 6,
        "vacuum_ratio": [None] * 6,
        "normalized_adf": [False] * 6,
        "height": None,
        "pins": None,
    }


def axes(burnup: list[float] | None = None) -> list[dict]:
    out = [
        {"name": "boron", "points": BORON, "b10": 0.199},
        {"name": "moderator", "points": MODERATOR, "density": [0.77, 0.743, 0.71]},
        {"name": "fuel", "points": FUEL},
    ]
    if burnup is not None:
        out.insert(0, {"name": "burnup", "points": burnup})
    return out


def table(g, burnup: list[float] | None = None) -> dict:
    """直積の順（最初の軸がいちばん遅く変わる）で、各状態が g(燃焼度, ホウ素, 減速材, 燃料) の状態"""
    points = burnup if burnup is not None else [None]
    states = [
        g(u, b, m, f)
        for u, b, m, f in itertools.product(points, BORON, MODERATOR, FUEL)
    ]
    return {"format": "xlc-table", "split": 1, "axes": axes(burnup), "states": states}


def additive(u, b, m, f) -> dict:
    return state(0.1 + 1e-5 * (b - 1300) + 1e-4 * (m - 565) - 2e-4 * (f - 565))


def boron_moderator_cross(u, b, m, f) -> dict:
    """νΣf2 に、ホウ素と減速材の積の項を持つ"""
    base = additive(u, b, m, f)
    return state(base["nu_fission"][1] + 1e-6 * (b - 1300) * (m - 565))


MODELS = [
    {"name": "1 軸ずつ", "keep": []},
    {"name": "ホウ素 × 減速材", "keep": [["boron", "moderator"]]},
    {"name": "直積", "keep": [["boron", "moderator", "fuel"]]},
]


def result(summary: dict, model: str, point: int = 0) -> dict:
    (t,) = summary["tables"]
    return t["points"][point]["models"][model]


class Reconstruction(unittest.TestCase):
    def test_one_axis_branches_reproduce_tables_without_cross_terms(self):
        s = cct.summarize([("t", table(additive), None, NOMINAL)], MODELS)
        self.assertLess(result(s, "1 軸ずつ")["max_abs_pcm"], 1e-9)
        self.assertLess(result(s, "1 軸ずつ")["max_relative"]["nu_fission"], 1e-12)

    def test_a_kept_cross_term_removes_its_error(self):
        s = cct.summarize([("t", table(boron_moderator_cross), None, NOMINAL)], MODELS)
        one = result(s, "1 軸ずつ")
        # 落とす積の項は、ホウ素と減速材の角の点で νΣf2 の ±1e-6 × 400 × 15 = ±6e-3。Δρ がいちばん大きいのは、そのうち
        # k∞ の小さい点
        corners = []
        for b, m, f in itertools.product((900.0, 1700.0), (550.0, 580.0), FUEL):
            truth = boron_moderator_cross(None, b, m, f)
            cross = 1e-6 * (b - 1300) * (m - 565)
            dropped = state(truth["nu_fission"][1] - cross)
            corners.append(((1 / truth["k_inf"] - 1 / dropped["k_inf"]) * 1e5, b, m, f))
        want = max(corners, key=lambda c: abs(c[0]))
        self.assertAlmostEqual(one["worst_pcm"], want[0], places=6)
        self.assertEqual(
            [one["worst"][a] for a in ("boron", "moderator", "fuel")], list(want[1:])
        )
        self.assertLess(result(s, "ホウ素 × 減速材")["max_abs_pcm"], 1e-9)
        self.assertLess(result(s, "直積")["max_abs_pcm"], 1e-9)

    def test_states_count_the_nominal_the_axes_and_the_kept_sub_grids(self):
        s = cct.summarize([("t", table(additive), None, NOMINAL)], MODELS)
        # 1 + 3 軸 × 2 点 = 7、ホウ素 × 減速材の 2 × 2 = 4 を足して 11、直積は 27
        self.assertEqual(result(s, "1 軸ずつ")["states"], 7)
        self.assertEqual(result(s, "ホウ素 × 減速材")["states"], 11)
        self.assertEqual(result(s, "直積")["states"], 27)
        self.assertEqual(s["tables"][0]["grid_states"], 27)

    def test_every_burnup_point_is_its_own_nominal(self):
        def g(u, b, m, f):
            base = boron_moderator_cross(u, b, m, f)
            return state(base["nu_fission"][1] + 0.01 * u)

        s = cct.summarize([("t", table(g, [0.0, 10.0]), None, NOMINAL)], MODELS)
        (t,) = s["tables"]
        self.assertEqual([p["burnup"] for p in t["points"]], [0.0, 10.0])
        for p in t["points"]:
            self.assertLess(p["models"]["ホウ素 × 減速材"]["max_abs_pcm"], 1e-9)
            self.assertGreater(p["models"]["1 軸ずつ"]["max_abs_pcm"], 1.0)


class Rods(unittest.TestCase):
    def rodded(self, cross: float):
        """制御棒を入れた状態: 熱群の吸収が 0.02 増え、ホウ素に比例する交差の項 cross × (b − 1300) を持つ"""

        def g(u, b, m, f):
            base = additive(u, b, m, f)
            return state(base["nu_fission"][1], 0.1 + cross * (b - 1300))

        return table(g)

    def models(self):
        return [
            {"name": "1 軸ずつ", "keep": []},
            {"name": "ホウ素 × 制御棒", "keep": [["boron", "rods"]]},
        ]

    def test_rod_branches_at_the_nominal_state_reproduce_a_constant_rod_effect(self):
        s = cct.summarize(
            [("t", table(additive), self.rodded(0.0), NOMINAL)], self.models()
        )
        self.assertLess(result(s, "1 軸ずつ")["max_abs_pcm"], 1e-9)
        # 制御棒の軸を足した 2 × 27 状態。1 軸ずつは 1 + 2 × 3 + 1
        self.assertEqual(s["tables"][0]["grid_states"], 54)
        self.assertEqual(result(s, "1 軸ずつ")["states"], 8)

    def test_a_rod_effect_that_depends_on_boron_needs_the_rod_cross_term(self):
        s = cct.summarize(
            [("t", table(additive), self.rodded(1e-5), NOMINAL)], self.models()
        )
        one = result(s, "1 軸ずつ")
        self.assertGreater(one["max_abs_pcm"], 1.0)
        self.assertEqual(one["worst"]["rods"], "in")
        self.assertLess(result(s, "ホウ素 × 制御棒")["max_abs_pcm"], 1e-9)
        # 制御棒の軸の点は 2 つなので、ホウ素 × 制御棒の部分格子は 2 × 1 = 2 状態を足す
        self.assertEqual(result(s, "ホウ素 × 制御棒")["states"], 10)

    def test_models_with_rods_are_skipped_for_tables_without_rods(self):
        s = cct.summarize([("t", table(additive), None, NOMINAL)], self.models())
        self.assertNotIn("ホウ素 × 制御棒", s["tables"][0]["points"][0]["models"])

    def test_rodded_tables_must_share_the_axes(self):
        rodded = self.rodded(0.0)
        rodded["axes"][0]["points"] = [900.0, 1300.0, 1800.0]
        with self.assertRaises(SystemExit):
            cct.summarize([("t", table(additive), rodded, NOMINAL)], self.models())


THINNING_BURNUP = [0.0, 1.0, 2.0, 4.0, 6.0, 10.0]


class Thinning(unittest.TestCase):
    """分岐を燃焼度の一部の点だけで解き、公称の状態との差を燃焼度について内挿する"""

    def branch_table(self, delta) -> dict:
        """公称の νΣf2 は燃焼度で 2 次に変わり、分岐の状態はそこに delta(燃焼度) × (ホウ素の偏差) を足す"""

        def g(u, b, m, f):
            base = additive(u, b, m, f)["nu_fission"][1] - 1e-3 * u * u
            return state(base + delta(u) * (b - 1300) / 400)

        return table(g, THINNING_BURNUP)

    def test_branch_differences_linear_in_burnup_are_reproduced_from_any_points(self):
        t = self.branch_table(lambda u: 1e-3 * u)
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        (sel,) = s["selections"]
        self.assertEqual(sel["removed"], [1.0, 2.0, 6.0])
        for method in cct.THINNING_METHODS:
            self.assertLess(sel["results"][method]["max_abs_pcm"], 1e-9, method)

    def test_removed_points_take_the_full_table_as_the_truth(self):
        # 差が 6 で折れる。選んだ 4 と 10 の間にある折れを、線形も Hermite（区間の点が 2 つなので線形）も拾えない
        t = self.branch_table(lambda u: 1e-3 * max(u - 6.0, 0.0))
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        (sel,) = s["selections"]
        linear = sel["results"]["linear"]
        self.assertGreater(linear["max_abs_pcm"], 1.0)
        # 間引いた点 1、2 は折れの手前で、差は 0
        self.assertLess(max(abs(r) for r in linear["rho_pcm"][:2]), 1e-9)
        # 誤差を出すのは公称でない状態だけで、そのうち最大の状態を書く
        self.assertNotEqual(linear["worst"]["boron"], 1300.0)
        self.assertEqual(linear["worst"]["burnup"], 6.0)

    def test_selections_must_be_points_from_the_first_to_the_last(self):
        t = self.branch_table(lambda u: 0.0)
        for burnup in ([0.0, 3.0, 10.0], [0.0, 4.0], [1.0, 4.0, 10.0]):
            with self.assertRaises(SystemExit, msg=str(burnup)):
                cct.thin("t", t, NOMINAL, [{"name": "bad", "burnup": burnup}])


class Errors(unittest.TestCase):
    def test_rodded_tables_are_checked_like_the_table(self):
        for change in (
            {"split": 4},
            {"format": "xlc-fewgroup"},
            {"states": table(additive)["states"][:-1]},
        ):
            rodded = {**table(additive), **change}
            with self.assertRaises(SystemExit, msg=str(change)):
                cct.summarize([("t", table(additive), rodded, NOMINAL)], MODELS)

    def test_thinning_must_remove_a_point(self):
        t = table(additive, [0.0, 1.0, 2.0])
        for burnup in ([0.0, 1.0, 2.0], []):
            with self.assertRaises(SystemExit, msg=str(burnup)):
                cct.thin("t", t, NOMINAL, [{"name": "bad", "burnup": burnup}])

    def test_thinning_without_error_names_the_first_removed_state(self):
        t = table(additive, [0.0, 1.0, 2.0])
        s = cct.thin("t", t, NOMINAL, [{"name": "2 点", "burnup": [0.0, 2.0]}])
        hermite = s["selections"][0]["results"]["hermite"]
        self.assertLess(hermite["max_abs_pcm"], 1e-9)
        self.assertEqual(hermite["worst"]["burnup"], 1.0)

    def test_migration_areas_need_two_groups(self):
        three = state(0.1)
        three["total"] = [0.5, 0.6, 1.0]
        with self.assertRaises(SystemExit):
            cct.migration_area(three)

    def test_items_that_are_zero_in_every_group_compare_as_zero(self):
        truth = {**state(0.1), "fission": [0.0, 0.0]}
        self.assertEqual(cct.relative_errors(truth, truth)["fission"], 0.0)
        with self.assertRaises(SystemExit):
            cct.relative_errors({**truth, "fission": [0.0, 1e-3]}, truth)

    def test_the_nominal_state_must_be_a_grid_point(self):
        for nominal in (
            {**NOMINAL, "boron": 1250.0},
            {"boron": 1300.0, "moderator": 565.0},
        ):
            with self.assertRaises(SystemExit, msg=str(nominal)):
                cct.summarize([("t", table(additive), None, nominal)], MODELS)

    def test_models_must_name_axes_of_the_table(self):
        with self.assertRaises(SystemExit):
            cct.summarize(
                [("t", table(additive), None, NOMINAL)],
                [{"name": "bad", "keep": [["boron", "density"]]}],
            )


class Report(unittest.TestCase):
    def test_the_report_needs_the_inputs(self):
        s = cct.summarize([("t", table(additive), None, NOMINAL)], MODELS)
        with self.assertRaises(SystemExit):
            cct.report(s)
        lines = cct.report({**s, "inputs": [{"path": "a.msgpack", "md5": "0" * 32}]})
        self.assertIn("a.msgpack 00000000", lines[0])


class Main(unittest.TestCase):
    def test_a_file_read_twice_is_listed_once_in_the_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "t.msgpack").write_bytes(msgpack.packb(table(additive)))
            nominal = "nominal = { boron = 1300, moderator = 565, fuel = 565 }"
            (root / "c.toml").write_text(
                f'[[table]]\nname = "a"\nfile = "t.msgpack"\n{nominal}\n'
                f'[[table]]\nname = "b"\nfile = "t.msgpack"\n{nominal}\n'
                '[[model]]\nname = "1 軸ずつ"\nkeep = []\n'
            )
            out = root / "s.json"
            argv = sys.argv
            sys.argv = ["x", "--json", str(out), str(root / "c.toml")]
            try:
                cct.main()
            finally:
                sys.argv = argv
            paths = [i["path"] for i in json.loads(out.read_text())["inputs"]]
            self.assertEqual(len(paths), len(set(paths)))
            self.assertEqual(len(paths), 2)


if __name__ == "__main__":
    unittest.main()
