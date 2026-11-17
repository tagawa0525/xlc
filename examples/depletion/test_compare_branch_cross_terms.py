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
    sizes = [len(a["points"]) for a in axes(burnup)]
    return {
        "format": "xlc-table",
        "split": 1,
        "axes": axes(burnup),
        "points": [list(p) for p in itertools.product(*[range(n) for n in sizes])],
        "states": states,
    }


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

    def test_errors_are_split_by_the_axes_off_their_nominal_values(self):
        # ホウ素だけを変えた状態の差は燃焼度に線形で、ホウ素と減速材の両方を変えた状態だけが 6 で折れる差を持つ
        def g(u, b, m, f):
            base = additive(u, b, m, f)["nu_fission"][1] - 1e-3 * u * u
            base += 1e-3 * u * (b - 1300) / 400
            if b != 1300 and m != 565:
                base += 1e-3 * max(u - 6.0, 0.0)
            return state(base)

        t = table(g, THINNING_BURNUP)
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        for method in cct.THINNING_METHODS:
            result = s["selections"][0]["results"][method]
            kinds = result["kinds"]
            self.assertEqual(
                list(kinds),
                [
                    "boron",
                    "moderator",
                    "fuel",
                    "boron+moderator",
                    "boron+fuel",
                    "moderator+fuel",
                    "boron+moderator+fuel",
                ],
                method,
            )
            for kind, r in kinds.items():
                if "boron" in kind and "moderator" in kind:
                    self.assertGreater(r["max_abs_pcm"], 1.0, (method, kind))
                    self.assertEqual(r["worst"]["burnup"], 6.0, (method, kind))
                else:
                    self.assertLess(r["max_abs_pcm"], 1e-9, (method, kind))
            self.assertEqual(
                result["max_abs_pcm"],
                max(r["max_abs_pcm"] for r in kinds.values()),
                method,
            )
        lines = cct.report({"inputs": [], "tables": [], "thinning": [s]})
        self.assertTrue(any("boron+moderator:" in line for line in lines))

    def test_reports_of_thinning_without_kinds_ask_to_rerun(self):
        t = self.branch_table(lambda u: 1e-3 * u)
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        for r in s["selections"][0]["results"].values():
            del r["kinds"]
        with self.assertRaisesRegex(SystemExit, "rerun"):
            cct.report({"inputs": [], "tables": [], "thinning": [s]})

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


PAIRS = [["boron", "moderator"], ["boron", "fuel"], ["moderator", "fuel"]]


def lattice_table(burnup: list[float] | None = None) -> dict:
    """lattice が書く形の直積のテーブル（点の番号、公称の値、内挿の変数を持つ。ADR-0200）"""
    t = table(boron_moderator_cross, burnup)
    for a in t["axes"]:
        if a["name"] in NOMINAL:
            a["nominal"] = NOMINAL[a["name"]]
            a["transform"] = "sqrt" if a["name"] == "fuel" else "linear"
    t["points"] = [
        list(p)
        for p in itertools.product(*[range(len(a["points"])) for a in t["axes"]])
    ]
    return t


class Products(unittest.TestCase):
    def test_compared_tables_must_hold_the_product_points(self):
        # 解いた点だけのテーブル（ADR-0200）は、直積の状態の番号で読めない
        t = table(additive)
        t["points"][0], t["points"][1] = t["points"][1], t["points"][0]
        with self.assertRaises(SystemExit):
            cct.table_grid("a", t, NOMINAL)


class Subgrids(unittest.TestCase):
    def test_pair_subgrids_keep_the_points_with_a_nominal_axis(self):
        product = lattice_table()
        out = cct.subgrid_table(product, PAIRS)
        # 1 + 6 + 12 点。3 つの軸がどれも公称でない点だけを除く
        self.assertEqual(len(out["points"]), 19)
        self.assertTrue(all(1 in p for p in out["points"]))
        for p, s in zip(out["points"], out["states"], strict=True):
            self.assertEqual(s, product["states"][product["points"].index(p)])
        self.assertEqual(out["axes"], product["axes"])

    def test_one_axis_branches_keep_the_nominal_state_and_the_axes(self):
        out = cct.subgrid_table(lattice_table(), [])
        self.assertEqual(len(out["points"]), 7)

    def test_burnup_and_rod_points_keep_every_value(self):
        burnt = lattice_table([0.0, 10.0])
        out = cct.subgrid_table(burnt, PAIRS)
        self.assertEqual(len(out["points"]), 2 * 19)
        rodded = lattice_table()
        rodded["axes"].append({"name": "rod", "values": ["out", "aic"]})
        rodded["points"] = [p + [r] for p in rodded["points"] for r in range(2)]
        rodded["states"] = [s for s in rodded["states"] for _ in range(2)]
        out = cct.subgrid_table(rodded, PAIRS)
        self.assertEqual(len(out["points"]), 2 * 19)

    def test_tables_must_be_products_with_nominal_values(self):
        bare = lattice_table()
        del bare["axes"][0]["nominal"]
        with self.assertRaises(SystemExit):
            cct.subgrid_table(bare, PAIRS)
        short = lattice_table()
        short["points"] = short["points"][:-1]
        short["states"] = short["states"][:-1]
        with self.assertRaises(SystemExit):
            cct.subgrid_table(short, PAIRS)
        with self.assertRaises(SystemExit):
            cct.subgrid_table(lattice_table(), [["boron", "rod"]])

    def test_the_command_writes_the_subgrid_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "product.msgpack"
            source.write_bytes(msgpack.packb(lattice_table()))
            output = root / "pairs.msgpack"
            argv = sys.argv
            sys.argv = [
                "compare_branch_cross_terms.py",
                "--subgrids",
                str(source),
                str(output),
                "boron,moderator",
                "boron,fuel",
                "moderator,fuel",
            ]
            try:
                cct.main()
            finally:
                sys.argv = argv
            self.assertEqual(len(msgpack.unpackb(output.read_bytes())["points"]), 19)


def micro_nuclide(
    total: float, below: float = 0.5, derivative: float | None = None
) -> dict:
    """lattice が書く形の 2 群の核種の微視的断面積（ADR-0216）。`below` は境界より下の核分裂の割合、`derivative` は
    数密度の微分の全断面積の 1 群目（ADR-0221）"""
    out = {
        "total": [total, 2 * total],
        "nu_fission": [0.1, 0.4],
        "fission": [0.04, 0.16],
        "scatter": [[total / 2, 0.01], [0.0, total]],
        "reactions": [
            {"mt": 102, "channels": [[total / 10, total / 5]]},
            {"mt": 17, "channels": [[0.0, 0.0]]},
        ],
        "thermal_fission": below,
    }
    if derivative is not None:
        out["derivative"] = {
            "total": [derivative, 1.0],
            "nu_fission": [0.0, 0.0],
            "fission": [0.0, 0.0],
            "scatter": [[0.0, 0.0], [0.0, 0.0]],
            "rate": {
                "fission": [0.0, 0.0],
                "reactions": [{"mt": 102, "channels": [[derivative, 0.0]]}],
                "thermal_fission": 0.0,
            },
        }
    return out


def with_micro(t: dict, nuclides) -> dict:
    """`nuclides(燃焼度, 状態の点の番号)` を状態ごとの微視的断面積にしたテーブル"""
    burnup = t["axes"][0]["points"]
    t["micro"] = {
        "nuclides": ["U-235", "Pu-240"],
        "densities": [[1e-4, 1e-6] for _ in burnup],
        "states": [nuclides(burnup[p[0]], tuple(p[1:])) for p in t["points"]],
    }
    return t


# 状態の点の番号の、公称の状態（NOMINAL の点）
AT_NOMINAL = (1, 1, 1)


class ThinningSubgrids(unittest.TestCase):
    """ADR-0200 の部分格子のテーブル（2 軸の組すべての点を、すべての燃焼点で解いたもの）で燃焼点を間引く"""

    def burnt(self) -> dict:
        def g(u, b, m, f):
            base = additive(u, b, m, f)["nu_fission"][1] - 1e-3 * u * u
            return state(base + 1e-3 * max(u - 3.0, 0.0) * (b - 1300 + m - 565) / 400)

        t = table(g, THINNING_BURNUP)
        for a in t["axes"]:
            if a["name"] in NOMINAL:
                a["nominal"] = NOMINAL[a["name"]]
        return t

    def test_subgrid_tables_are_thinned_like_the_same_states_of_the_product(self):
        product = self.burnt()
        pairs = cct.subgrid_table(product, PAIRS)
        self.assertLess(len(pairs["points"]), len(product["points"]))
        selection = [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}]
        full = cct.thin("t", product, NOMINAL, selection)["selections"][0]
        sub = cct.thin("t", pairs, NOMINAL, selection)["selections"][0]
        for method in cct.THINNING_METHODS:
            kinds = sub["results"][method]["kinds"]
            self.assertNotIn("boron+moderator+fuel", kinds, method)
            for kind, r in kinds.items():
                self.assertEqual(
                    r, full["results"][method]["kinds"][kind], (method, kind)
                )
            self.assertGreater(kinds["boron+moderator"]["max_abs_pcm"], 1.0, method)

    def test_tables_must_solve_every_state_at_every_burnup(self):
        t = cct.subgrid_table(self.burnt(), PAIRS)
        last = len(t["points"]) - 1
        del t["points"][last]
        del t["states"][last]
        with self.assertRaisesRegex(SystemExit, "every burnup"):
            cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])

    def test_tables_must_hold_the_nominal_state(self):
        t = self.burnt()
        keep = [n for n, p in enumerate(t["points"]) if tuple(p[1:]) != AT_NOMINAL]
        t["points"] = [t["points"][n] for n in keep]
        t["states"] = [t["states"][n] for n in keep]
        with self.assertRaisesRegex(SystemExit, "nominal"):
            cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])

    def test_nominal_values_of_the_table_must_be_the_given_ones(self):
        t = self.burnt()
        t["axes"][1]["nominal"] = 900.0
        with self.assertRaisesRegex(SystemExit, "nominal"):
            cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])


class ThinningMicro(unittest.TestCase):
    """微視的断面積と数密度の微分を、少数群定数と同じ内挿（公称の状態との差を燃焼度について内挿する）で比べる"""

    def thinned(self, nuclides) -> dict:
        t = with_micro(Thinning().branch_table(lambda u: 1e-3 * u), nuclides)
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        return s["selections"][0]["results"]

    def test_differences_linear_in_burnup_give_no_micro_error(self):
        def nuclides(u, at):
            off = 0.0 if at == AT_NOMINAL else 1.0
            return [
                micro_nuclide(10 + u * u + off * u),
                micro_nuclide(5 + off * u, derivative=2 + u * off),
            ]

        for method, r in self.thinned(nuclides).items():
            for item, e in r["micro"].items():
                self.assertLess(e["max_relative"], 1e-12, (method, item))
            for item, e in r["derivative"].items():
                self.assertLess(e["max_relative"], 1e-12, (method, item))

    def test_micro_errors_name_the_item_the_nuclide_and_the_state(self):
        # Pu-240 の全断面積の差だけが 6 で折れる
        def nuclides(u, at):
            off = 0.0 if at == AT_NOMINAL else 1.0
            return [micro_nuclide(10.0), micro_nuclide(5 + off * max(u - 6.0, 0.0))]

        linear = self.thinned(nuclides)["linear"]
        total = linear["micro"]["total"]
        self.assertGreater(total["max_relative"], 1e-3)
        self.assertEqual(total["nuclide"], "Pu-240")
        self.assertEqual(total["worst"]["burnup"], 6.0)
        self.assertEqual(linear["micro"]["reactions"]["nuclide"], "Pu-240 mt102")
        self.assertLess(linear["micro"]["nu_fission"]["max_relative"], 1e-12)
        # 状態の種類ごとにも書き、全体はその最大
        kinds = linear["kinds"]
        self.assertEqual(
            total["max_relative"],
            max(k["micro"]["total"]["max_relative"] for k in kinds.values()),
        )

    def test_fractions_below_the_cutoff_interpolate_as_cross_sections(self):
        # 境界より下の核分裂の断面積（割合 × 群の和の σf）の差が燃焼度に線形なら、割合が線形でなくても誤差は 0（core の
        # `interpolate_micro` と同じく、断面積で内挿して割合に戻す）
        def nuclides(u, at):
            off = 0.0 if at == AT_NOMINAL else 1.0
            n = micro_nuclide(10.0, below=0.5)
            scale = 1 + 0.1 * u * off
            n["fission"] = [0.04 * scale, 0.16 * scale]
            n["thermal_fission"] = (0.1 + 0.01 * u * off) / (0.2 * scale)
            return [n, micro_nuclide(5.0)]

        for method, r in self.thinned(nuclides).items():
            self.assertLess(
                r["micro"]["thermal_fission"]["max_relative"], 1e-12, method
            )

    def test_density_derivatives_are_compared_to_the_largest_value_of_their_item(self):
        # 微分の全断面積は [d, 1]。d が 0 を横切っても、要素ごとの比でなく項目の最大の大きさとの比で書く
        def nuclides(u, at):
            off = 0.0 if at == AT_NOMINAL else 1.0
            d = -0.5 + 0.1 * u + off * 1e-3 * max(u - 6.0, 0.0)
            return [micro_nuclide(10.0), micro_nuclide(5.0, derivative=d)]

        linear = self.thinned(nuclides)["linear"]
        total = linear["derivative"]["total"]
        # 6 の差 0 から、4 と 10 の線形の内挿で 2/6 × 4e-3 を外す。項目の最大の大きさは 1
        self.assertAlmostEqual(total["max_relative"], 4e-3 / 3, places=12)
        self.assertEqual(total["nuclide"], "Pu-240")
        self.assertEqual(
            linear["derivative"]["rate_reactions"]["nuclide"], "Pu-240 mt102"
        )

    def test_micro_values_off_zero_where_the_table_is_zero_are_errors(self):
        def nuclides(u, at):
            n = micro_nuclide(10.0)
            if at != AT_NOMINAL and u != 6.0:
                n["reactions"][1]["channels"] = [[1e-3, 0.0]]
            return [n, micro_nuclide(5.0)]

        with self.assertRaisesRegex(SystemExit, "zero"):
            self.thinned(nuclides)

    def test_tables_without_micro_data_compare_only_the_few_group_constants(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        s = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}])
        for r in s["selections"][0]["results"].values():
            self.assertNotIn("micro", r)
            self.assertNotIn("derivative", r)


# 燃焼点ごとの全体 100 秒、自己遮蔽を解く状態 3 秒、ホウ素を変えて自己遮蔽を再利用する状態 1 秒
# 3 つの軸の直積の状態の種類（`state_kind`）
KINDS = [
    "boron",
    "moderator",
    "fuel",
    "boron+moderator",
    "boron+fuel",
    "moderator+fuel",
    "boron+moderator+fuel",
]
TIMING = {"total": 100.0, "solve": 3.0, "reuse": 1.0, "reuse-axis": "boron"}


class ThinningTime(unittest.TestCase):
    """燃焼点を間引いたテーブルを作る時間の見込み（公称の状態はすべての燃焼点で解く）"""

    def summary(self, plans=None) -> dict:
        t = Thinning().branch_table(lambda u: 1e-3 * max(u - 3.0, 0.0))
        return cct.thin(
            "t",
            t,
            NOMINAL,
            [
                {"name": "3 点", "burnup": [0.0, 4.0, 10.0]},
                {"name": "4 点", "burnup": [0.0, 2.0, 4.0, 10.0]},
            ],
            TIMING,
            plans or [],
        )

    def test_hours_scale_the_off_nominal_states_by_the_chosen_points(self):
        s = self.summary()
        # 27 状態（3 × 3 × 3）のうち、ホウ素を公称から変えた 18 状態が再利用、残りの 9 状態（公称を含む）が自己遮蔽を解く。
        # 燃焼点ごとの状態によらない秒は 100 − (9 × 3 + 18 × 1) = 55
        self.assertAlmostEqual(s["other_seconds"], 55.0)
        self.assertAlmostEqual(s["hours"], 6 * 100.0 / 3600)
        three, four = s["selections"]
        self.assertAlmostEqual(
            three["hours"], (6 * (55 + 3) + 3 * (8 * 3 + 18 * 1)) / 3600
        )
        self.assertAlmostEqual(
            four["hours"], (6 * (55 + 3) + 4 * (8 * 3 + 18 * 1)) / 3600
        )
        kinds = three["results"]["hermite"]["kinds"]
        self.assertEqual(kinds["boron"]["states"], 2)
        self.assertAlmostEqual(kinds["boron"]["seconds_per_point"], 2 * 1.0)
        self.assertAlmostEqual(kinds["moderator+fuel"]["seconds_per_point"], 4 * 3.0)

    def test_plans_take_each_kind_from_its_selection(self):
        kinds = KINDS
        choice = dict.fromkeys(kinds, "4 点") | {
            "boron": "3 点",
            "fuel": "3 点",
            "boron+fuel": "3 点",
        }
        s = self.summary([{"name": "混ぜる", "kinds": choice}])
        (plan,) = s["plans"]
        three, four = s["selections"]
        for method in cct.THINNING_METHODS:
            expected = max(
                (three if choice[k] == "3 点" else four)["results"][method]["kinds"][k][
                    "max_abs_pcm"
                ]
                for k in kinds
            )
            self.assertEqual(plan["results"][method]["max_abs_pcm"], expected, method)
        seconds = sum(
            len((three if choice[k] == "3 点" else four)["burnup"])
            * three["results"]["hermite"]["kinds"][k]["seconds_per_point"]
            for k in kinds
        )
        self.assertAlmostEqual(plan["hours"], (6 * (55 + 3) + seconds) / 3600)

    def test_plans_must_choose_a_selection_for_every_kind(self):
        with self.assertRaisesRegex(SystemExit, "kind"):
            self.summary([{"name": "足りない", "kinds": {"boron": "3 点"}}])

    def test_reusing_states_need_the_points_of_the_state_they_reuse(self):
        # ホウ素と減速材を変えた状態は、減速材だけを変えた状態の自己遮蔽を再利用するので、その燃焼点の中で解く
        kinds = KINDS
        choice = dict.fromkeys(kinds, "3 点") | {"boron+moderator": "4 点"}
        with self.assertRaisesRegex(SystemExit, "reuse"):
            self.summary([{"name": "広い", "kinds": choice}])
        # 公称の状態はすべての燃焼点で解くので、ホウ素だけを変えた状態はどの燃焼点でもよい
        self.summary(
            [
                {
                    "name": "ホウ素",
                    "kinds": dict.fromkeys(kinds, "3 点") | {"boron": "4 点"},
                }
            ]
        )


class MainThinningOnly(unittest.TestCase):
    def test_inputs_may_hold_only_the_thinning(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            t = Thinning().branch_table(lambda u: 1e-3 * u)
            (root / "t.msgpack").write_bytes(msgpack.packb(t))
            (root / "c.toml").write_text(
                '[[thinning]]\nname = "a"\nfile = "t.msgpack"\n'
                "nominal = { boron = 1300, moderator = 565, fuel = 565 }\n"
                'selections = [{ name = "3 点", burnup = [0, 4, 10] }]\n'
                "[thinning.timing]\ntotal = 100.0\nsolve = 3.0\nreuse = 1.0\n"
                'reuse-axis = "boron"\n'
            )
            out = root / "s.json"
            argv = sys.argv
            sys.argv = ["x", "--json", str(out), str(root / "c.toml")]
            try:
                cct.main()
            finally:
                sys.argv = argv
            summary = json.loads(out.read_text())
            self.assertEqual(summary["tables"], [])
            self.assertIn("hours", summary["thinning"][0]["selections"][0])


if __name__ == "__main__":
    unittest.main()
