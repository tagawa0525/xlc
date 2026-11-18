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
# 別のテーブルの状態の数で見込む時間（`states` はどの種類も書き、`total` でなく `other`）
ESTIMATE = {
    "other": 55.0,
    "solve": 3.0,
    "reuse": 1.0,
    "reuse-axis": "boron",
    "states": dict.fromkeys(KINDS, 0) | {"boron": 5, "moderator+fuel": 36},
}


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

    def test_given_states_and_other_seconds_estimate_another_table(self):
        # 少ない状態で測った間引きの誤差に、本番の形の状態の数（`states`）と、状態によらない秒（`other`）で時間を見込む
        timing = ESTIMATE
        t = Thinning().branch_table(lambda u: 1e-3 * max(u - 3.0, 0.0))
        s = cct.thin(
            "t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}], timing
        )
        self.assertAlmostEqual(s["other_seconds"], 55.0)
        self.assertAlmostEqual(s["hours"], 6 * (55 + 3 + 5 * 1 + 36 * 3) / 3600)
        (three,) = s["selections"]
        self.assertAlmostEqual(
            three["hours"], (6 * (55 + 3) + 3 * (5 * 1 + 36 * 3)) / 3600
        )
        kinds = three["results"]["hermite"]["kinds"]
        self.assertAlmostEqual(kinds["boron"]["seconds_per_point"], 5.0)
        self.assertAlmostEqual(kinds["moderator+fuel"]["seconds_per_point"], 108.0)
        self.assertAlmostEqual(kinds["boron+moderator+fuel"]["seconds_per_point"], 0.0)
        # 誤差を測った状態の数はテーブルのまま
        self.assertEqual(kinds["boron"]["states"], 2)

    def test_timing_needs_either_the_total_or_the_other_seconds(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        selection = [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}]
        for timing in (
            {"solve": 3.0, "reuse": 1.0, "reuse-axis": "boron"},
            TIMING | {"other": 55.0},
        ):
            with self.assertRaisesRegex(SystemExit, "total"):
                cct.thin("t", t, NOMINAL, selection, timing)

    def test_given_states_must_be_counts(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        selection = [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}]
        for count in (5.9, -3, True):
            timing = ESTIMATE | {"states": ESTIMATE["states"] | {"boron": count}}
            with self.assertRaisesRegex(SystemExit, "count", msg=count):
                cct.thin("t", t, NOMINAL, selection, timing)

    def test_the_kinds_carry_the_state_counts_of_the_time_estimate(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        selection = [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}]
        (s,) = cct.thin("t", t, NOMINAL, selection, ESTIMATE)["selections"]
        kinds = s["results"]["hermite"]["kinds"]
        self.assertEqual(kinds["boron"]["timed_states"], 5)
        self.assertEqual(kinds["moderator"]["timed_states"], 0)
        (s,) = cct.thin("t", t, NOMINAL, selection, TIMING)["selections"]
        kinds = s["results"]["hermite"]["kinds"]
        self.assertEqual(kinds["boron"]["timed_states"], kinds["boron"]["states"])

    def test_given_states_must_be_kinds_of_the_table(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        timing = ESTIMATE | {"states": ESTIMATE["states"] | {"boron+rods": 3}}
        with self.assertRaisesRegex(SystemExit, "kind"):
            cct.thin(
                "t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}], timing
            )

    def test_given_states_must_list_every_kind_of_the_table(self):
        # 書き漏らした種類を 0 にすると、その種類の時間が黙って落ちる
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        states = dict(ESTIMATE["states"])
        del states["moderator+fuel"]
        with self.assertRaisesRegex(SystemExit, "kind"):
            cct.thin(
                "t",
                t,
                NOMINAL,
                [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}],
                ESTIMATE | {"states": states},
            )

    def test_given_states_need_the_other_seconds_not_the_total(self):
        # `total` は誤差を測ったテーブルの全体の秒で、別のテーブルの状態の秒を引いても状態によらない秒にならない
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        timing = TIMING | {"states": ESTIMATE["states"]}
        with self.assertRaisesRegex(SystemExit, "other"):
            cct.thin(
                "t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}], timing
            )

    def test_unknown_timing_keys_are_errors(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        timing = TIMING | {"state": ESTIMATE["states"]}
        with self.assertRaisesRegex(SystemExit, "key"):
            cct.thin(
                "t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}], timing
            )

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
            self.assertEqual(summary["nominal"], [])
            self.assertIn("hours", summary["thinning"][0]["selections"][0])


class NominalBurnupPoints(unittest.TestCase):
    """公称の状態を燃焼度の一部の点だけで解いたとして、その点から燃焼度について内挿した公称の状態を、全点のテーブルの公称の状態と
    比べる（公称でない状態は使わない）"""

    def kinked(self, nominal_kink: float, off_kink: float = 0.0) -> dict:
        """公称の νΣf2 は燃焼度で線形に変わり、6 で nominal_kink × (u − 6) の折れを持つ。公称でない状態は off_kink の折れを足す"""

        def g(u, b, m, f):
            base = 0.1 - 1e-3 * u + nominal_kink * max(u - 6.0, 0.0)
            if (b, m, f) != (1300.0, 565.0, 565.0):
                base += 1e-5 * (b - 1300) + off_kink * max(u - 6.0, 0.0)
            return state(base)

        return table(g, THINNING_BURNUP)

    def compared(self, t: dict, burnup=(0.0, 4.0, 10.0)) -> dict:
        s = cct.nominal_points(
            "t", t, NOMINAL, [{"name": "3 点", "burnup": list(burnup)}]
        )
        (sel,) = s["selections"]
        return sel

    def test_nominal_states_linear_in_burnup_are_reproduced_from_any_points(self):
        sel = self.compared(self.kinked(0.0, off_kink=1e-3))
        self.assertEqual(sel["removed"], [1.0, 2.0, 6.0])
        for method in cct.THINNING_METHODS:
            r = sel["results"][method]
            self.assertLess(r["max_abs_pcm"], 1e-9, method)
            self.assertEqual(len(r["rho_pcm"]), 3, method)

    def test_removed_points_compare_the_nominal_state_of_the_full_table(self):
        # 公称の状態が 6 で折れる。選んだ 4 と 10 の間の折れを、線形も Hermite（区間の点が 2 つなので線形）も拾えない
        sel = self.compared(self.kinked(1e-3))
        linear = sel["results"]["linear"]
        self.assertGreater(linear["max_abs_pcm"], 1.0)
        self.assertEqual(linear["worst_burnup"], 6.0)
        self.assertLess(max(abs(r) for r in linear["rho_pcm"][:2]), 1e-9)
        # 公称でない状態の折れは、公称の状態の比較に効かない
        other = self.compared(self.kinked(1e-3, off_kink=5e-3))
        self.assertEqual(other["results"], sel["results"])

    def test_relative_errors_name_the_few_group_items(self):
        linear = self.compared(self.kinked(1e-3))["results"]["linear"]
        relative = linear["max_relative"]
        self.assertEqual(sorted(relative), sorted(cct.RELATIVE))
        self.assertGreater(relative["nu_fission"], 1e-3)
        self.assertLess(relative["absorption"], 1e-12)

    def test_micro_errors_name_the_item_the_nuclide_and_the_burnup(self):
        # 公称の状態の Pu-240 の全断面積だけが 6 で折れる（公称でない状態の折れは効かない）
        def nuclides(u, at):
            kink = max(u - 6.0, 0.0) * (1.0 if at == AT_NOMINAL else 3.0)
            return [micro_nuclide(10 + u), micro_nuclide(5 + kink)]

        t = with_micro(self.kinked(0.0), nuclides)
        linear = self.compared(t)["results"]["linear"]
        total = linear["micro"]["total"]
        self.assertGreater(total["max_relative"], 1e-3)
        self.assertEqual(total["nuclide"], "Pu-240")
        self.assertEqual(total["worst"]["burnup"], 6.0)
        self.assertLess(linear["micro"]["nu_fission"]["max_relative"], 1e-12)

    def test_selections_must_be_points_from_the_first_to_the_last_leaving_some_out(
        self,
    ):
        t = self.kinked(0.0)
        for burnup in ([1.0, 4.0, 10.0], [0.0, 4.0], [0.0, 3.0, 10.0], THINNING_BURNUP):
            with self.assertRaisesRegex(SystemExit, "increasing", msg=burnup):
                self.compared(t, burnup)

    def test_inputs_may_hold_the_nominal_points_and_report_them(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "t.msgpack").write_bytes(msgpack.packb(self.kinked(1e-3)))
            (root / "c.toml").write_text(
                '[[nominal]]\nname = "a"\nfile = "t.msgpack"\n'
                "nominal = { boron = 1300, moderator = 565, fuel = 565 }\n"
                'selections = [{ name = "3 点", burnup = [0, 4, 10] }]\n'
            )
            out = root / "s.json"
            argv = sys.argv
            sys.argv = ["x", "--json", str(out), str(root / "c.toml")]
            try:
                cct.main()
            finally:
                sys.argv = argv
            summary = json.loads(out.read_text())
            (entry,) = summary["nominal"]
            self.assertEqual(entry["name"], "a")
            self.assertEqual(entry["selections"][0]["removed"], [1.0, 2.0, 6.0])
            lines = cct.report(summary)
            self.assertTrue(
                any("公称の状態" in line and "3 点" in line for line in lines)
            )


class ThinningNominalBurnups(unittest.TestCase):
    """間引きの真値のテーブルより少ない燃焼点で公称の状態を解く、本番の形のテーブルの時間の見込み"""

    def thinned(self, nominal_burnup, burnup=(0.0, 4.0, 10.0)) -> dict:
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        timing = ESTIMATE | {"nominal-burnup": nominal_burnup}
        return cct.thin(
            "t", t, NOMINAL, [{"name": "3 点", "burnup": list(burnup)}], timing
        )

    def test_nominal_points_count_the_nominal_and_other_seconds(self):
        s = self.thinned([0, 2, 4, 10])
        self.assertEqual(s["nominal_burnup"], [0.0, 2.0, 4.0, 10.0])
        self.assertAlmostEqual(s["hours"], 4 * (55 + 3 + 5 * 1 + 36 * 3) / 3600)
        (three,) = s["selections"]
        self.assertAlmostEqual(
            three["hours"], (4 * (55 + 3) + 3 * (5 * 1 + 36 * 3)) / 3600
        )

    def test_tables_without_nominal_points_solve_the_nominal_state_at_every_point(self):
        t = Thinning().branch_table(lambda u: 1e-3 * u)
        s = cct.thin(
            "t", t, NOMINAL, [{"name": "3 点", "burnup": [0.0, 4.0, 10.0]}], ESTIMATE
        )
        self.assertEqual(s["nominal_burnup"], THINNING_BURNUP)

    def test_selections_must_be_nominal_points(self):
        # 公称でない状態は、公称の状態を解く燃焼点でしか解けない（公称の状態との差を持つため）
        with self.assertRaisesRegex(SystemExit, "nominal-burnup"):
            self.thinned([0, 2, 4, 10], burnup=(0.0, 6.0, 10.0))

    def test_selections_off_the_table_are_named_before_the_nominal_points(self):
        # 表に無い点を選んだ誤りを、公称の状態の点の外として出さない
        with self.assertRaisesRegex(SystemExit, "increasing points of the table"):
            self.thinned([0, 2, 4, 10], burnup=(0.0, 3.0, 10.0))

    def test_nominal_points_must_be_points_of_the_table_from_the_first_to_the_last(
        self,
    ):
        for points in ([0, 3, 10], [1, 4, 10], [0, 4], [10, 4, 0], [0, 4, 4, 10], 4):
            with self.assertRaisesRegex(SystemExit, "nominal-burnup", msg=points):
                self.thinned(points)


class CombinedBurnupPoints(unittest.TestCase):
    """公称の状態を一部の燃焼点で、公称でない状態をその中の一部の燃焼点で解いたとして、core の組み立て（公称の状態を公称の点から
    内挿し、公称でない状態の公称の状態との差を公称でない点から内挿して足す。ADR-0200）の値を、全点のテーブルの状態と比べる"""

    def kinked(self, nominal_kink: float, off_kink: float = 0.0, at: float = 6.0):
        """公称の νΣf2 は燃焼度で線形で、`at` で nominal_kink の折れを持つ。公称でない状態の差はホウ素の偏差に比例し、燃焼度で
        線形に変わり、6 で off_kink の折れを持つ"""

        def g(u, b, m, f):
            base = 0.1 - 1e-3 * u + nominal_kink * max(u - at, 0.0)
            if (b, m, f) != (1300.0, 565.0, 565.0):
                base += 1e-5 * (b - 1300) + 1e-6 * u + off_kink * max(u - 6.0, 0.0)
            return state(base)

        return table(g, THINNING_BURNUP)

    def combined(self, t, nominal, off, timing=None) -> dict:
        selections = [
            {"name": "公称", "burnup": list(nominal)},
            {"name": "公称でない", "burnup": list(off)},
        ]
        pairs = [{"nominal": "公称", "off-nominal": "公称でない"}]
        return cct.combine("t", t, NOMINAL, selections, pairs, timing)

    def test_states_linear_in_burnup_are_reproduced_from_any_points(self):
        s = self.combined(self.kinked(0.0), [0.0, 4.0, 10.0], [0.0, 4.0, 10.0])
        (pair,) = s["pairs"]
        self.assertEqual(pair["removed"], [1.0, 2.0, 6.0])
        for method in cct.THINNING_METHODS:
            self.assertLess(pair["results"][method]["max_abs_pcm"], 1e-9, method)

    def test_off_nominal_states_carry_the_error_of_the_nominal_interpolation(self):
        # 公称の状態だけが 6 で折れ、公称でない状態の差は線形。公称の点にも公称でない点にも無い 6 で、公称でない状態も公称の
        # 状態の内挿の誤差を持つ（間引きの比較は公称の状態を全点の値にするので、この誤差を持たない）
        t = self.kinked(1e-3)
        s = self.combined(t, [0.0, 4.0, 10.0], [0.0, 4.0, 10.0])
        thinned = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0, 4, 10]}])
        nominal = cct.nominal_points(
            "t", t, NOMINAL, [{"name": "3 点", "burnup": [0, 4, 10]}]
        )
        for method in cct.THINNING_METHODS:
            r = s["pairs"][0]["results"][method]
            self.assertLess(
                thinned["selections"][0]["results"][method]["max_abs_pcm"], 1e-9
            )
            for kind, k in r["kinds"].items():
                self.assertGreater(k["max_abs_pcm"], 1.0, (method, kind))
                self.assertEqual(k["worst"]["burnup"], 6.0, (method, kind))
            # 公称の状態の誤差は、公称の状態だけの比較と同じ
            expected = nominal["selections"][0]["results"][method]
            self.assertAlmostEqual(
                r["nominal"]["max_abs_pcm"], expected["max_abs_pcm"], places=9
            )
            self.assertEqual(r["nominal"]["worst_burnup"], 6.0)
            self.assertEqual(
                r["max_abs_pcm"],
                max(
                    [r["nominal"]["max_abs_pcm"]]
                    + [k["max_abs_pcm"] for k in r["kinds"].values()]
                ),
            )
            self.assertEqual(sorted(r["max_relative"]), sorted(cct.RELATIVE))
            self.assertGreater(r["max_relative"]["nu_fission"], 1e-3)

    def test_nominal_points_take_the_nominal_state_of_the_table(self):
        # 公称の状態は 1 で折れる。公称の点の 2 と 6 では公称の状態はテーブルの値で、誤差は間引きの比較と同じ。どちらの点にも
        # 無い 1 だけが、公称の状態の内挿の誤差（0 から 2 の線形）を足す
        t = self.kinked(1e-3, off_kink=1e-3, at=1.0)
        s = self.combined(t, [0.0, 2.0, 4.0, 6.0, 10.0], [0.0, 4.0, 10.0])
        thinned = cct.thin("t", t, NOMINAL, [{"name": "3 点", "burnup": [0, 4, 10]}])
        for method in cct.THINNING_METHODS:
            r = s["pairs"][0]["results"][method]
            expected = thinned["selections"][0]["results"][method]["rho_pcm"]
            self.assertEqual(len(r["rho_pcm"]), 3)
            for got, want in zip(r["rho_pcm"][1:], expected[1:], strict=True):
                self.assertAlmostEqual(got, want, places=9, msg=method)
            self.assertGreater(abs(r["rho_pcm"][0] - expected[0]), 1.0, method)
            self.assertLess(r["nominal"]["max_abs_pcm"] - abs(r["rho_pcm"][0]), 1e-9)

    def test_off_nominal_points_must_be_nominal_points(self):
        # core は公称でない状態の公称の状態との差を、その状態を解いた燃焼点で作る。公称の状態はテーブルの燃焼度の軸の点の
        # すべてで解く（ADR-0200）ので、公称でない状態の点は公称の点の中から選ぶ
        with self.assertRaisesRegex(SystemExit, "nominal"):
            self.combined(self.kinked(0.0), [0.0, 4.0, 10.0], [0.0, 6.0, 10.0])

    def test_pairs_name_selections_of_the_entry(self):
        selections = [{"name": "3 点", "burnup": [0, 4, 10]}]
        for pair in (
            {"nominal": "3 点", "off-nominal": "4 点"},
            {"nominal": "4 点", "off-nominal": "3 点"},
        ):
            with self.assertRaisesRegex(SystemExit, "4 点", msg=pair):
                cct.combine("t", self.kinked(0.0), NOMINAL, selections, [pair])

    def test_hours_count_the_nominal_points_and_the_off_nominal_points(self):
        s = self.combined(
            self.kinked(0.0), [0.0, 2.0, 4.0, 10.0], [0.0, 4.0, 10.0], ESTIMATE
        )
        (pair,) = s["pairs"]
        self.assertAlmostEqual(
            pair["hours"], (4 * (55 + 3) + 3 * (5 * 1 + 36 * 3)) / 3600
        )
        # 公称の状態の燃焼点は組が持つので、時間の見込みの `nominal-burnup` は書かない
        with self.assertRaisesRegex(SystemExit, "nominal-burnup"):
            self.combined(
                self.kinked(0.0),
                [0.0, 2.0, 4.0, 10.0],
                [0.0, 4.0, 10.0],
                ESTIMATE | {"nominal-burnup": [0, 2, 4, 10]},
            )

    def test_micro_errors_carry_the_nominal_interpolation(self):
        # 公称の状態の Pu-240 の全断面積だけが 6 で折れ、公称でない状態の差は線形
        def nuclides(u, at):
            off = 0.0 if at == AT_NOMINAL else u
            return [micro_nuclide(10 + u), micro_nuclide(5 + max(u - 6.0, 0.0) + off)]

        t = with_micro(self.kinked(0.0), nuclides)
        s = self.combined(t, [0.0, 4.0, 10.0], [0.0, 4.0, 10.0])
        linear = s["pairs"][0]["results"]["linear"]
        total = linear["micro"]["total"]
        self.assertGreater(total["max_relative"], 1e-3)
        self.assertEqual(total["nuclide"], "Pu-240")
        self.assertEqual(total["worst"]["burnup"], 6.0)
        self.assertIn("boron", total["worst"])
        self.assertLess(linear["micro"]["nu_fission"]["max_relative"], 1e-12)
        for kind in linear["kinds"].values():
            self.assertGreater(kind["micro"]["total"]["max_relative"], 1e-3)

    def test_inputs_may_hold_the_combined_points_and_report_them(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "t.msgpack").write_bytes(msgpack.packb(self.kinked(1e-3)))
            (root / "c.toml").write_text(
                '[[combined]]\nname = "a"\nfile = "t.msgpack"\n'
                "nominal = { boron = 1300, moderator = 565, fuel = 565 }\n"
                'selections = [{ name = "4 点", burnup = [0, 2, 4, 10] }, '
                '{ name = "3 点", burnup = [0, 4, 10] }]\n'
                'pairs = [{ nominal = "4 点", off-nominal = "3 点" }]\n'
                "[combined.timing]\nother = 55.0\nsolve = 3.0\nreuse = 1.0\n"
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
            (entry,) = summary["combined"]
            self.assertEqual(entry["name"], "a")
            (pair,) = entry["pairs"]
            self.assertEqual(pair["name"], "4 点 / 3 点")
            self.assertEqual(pair["nominal_burnup"], [0.0, 2.0, 4.0, 10.0])
            self.assertEqual(pair["off_nominal_burnup"], [0.0, 4.0, 10.0])
            self.assertIn("hours", pair)
            lines = cct.report(summary)
            self.assertTrue(any("4 点 / 3 点" in line for line in lines))


if __name__ == "__main__":
    unittest.main()
