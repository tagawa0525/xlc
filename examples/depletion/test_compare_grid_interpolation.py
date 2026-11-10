"""compare_grid_interpolation.py の要約（直積の分岐のテーブルの格子点を、その軸の両隣の格子点から core と同じ式で内挿したときの差と、
本番の間隔の誤差の見込み。ADR-0239）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_grid_interpolation.py'
"""

import itertools
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_grid_interpolation as cgi
from compare_burnup_interpolation import k_inf, unvector, vector

# 減速材と燃料は間隔が等しくない（減速材の 570 K は 560 K と 575 K から重み 2/3、燃料の 600 K は √T で内挿する）
BORON = [0.0, 650.0, 1300.0]
MODERATOR = [560.0, 570.0, 575.0]
FUEL = [560.0, 600.0, 750.0]
BURNUP = [0.0, 10.0]


def state(nu_fission2: float) -> dict:
    """2 群の状態（上方散乱なし）"""
    return {
        "k_inf": (0.005 + nu_fission2 * 0.02 / 0.1) / 0.03,
        "total": [0.5, 1.0],
        "absorption": [0.01, 0.08],
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


def axes() -> list[dict]:
    return [
        {"name": "burnup", "points": BURNUP},
        {
            "name": "boron",
            "points": BORON,
            "b10": 0.199,
            "nominal": 650.0,
            "transform": "linear",
        },
        {
            "name": "moderator",
            "points": MODERATOR,
            "density": [0.75, 0.73, 0.72],
            "nominal": 570.0,
            "transform": "linear",
        },
        {"name": "fuel", "points": FUEL, "nominal": 600.0, "transform": "sqrt"},
    ]


def table(g) -> dict:
    """直積の順（最初の軸がいちばん遅く変わる）で、各状態の νΣf2 が g(燃焼度, ホウ素, 減速材, 燃料) の状態"""
    states = [
        state(g(u, b, m, f))
        for u, b, m, f in itertools.product(BURNUP, BORON, MODERATOR, FUEL)
    ]
    sizes = [len(a["points"]) for a in axes()]
    return {
        "format": "xlc-table",
        "split": 1,
        "axes": axes(),
        "points": [list(p) for p in itertools.product(*[range(n) for n in sizes])],
        "states": states,
    }


def linear_in_transforms(u, b, m, f) -> float:
    """どの軸についても内挿の変数（ppm、温度、√T）の 1 次"""
    return 0.1 + 1e-6 * u + 1e-5 * b + 1e-4 * m - 2e-3 * math.sqrt(f)


def moderator_square(u, b, m, f) -> float:
    """減速材の温度の 2 次の項を持つ（燃焼度とともに大きくなる）"""
    return linear_in_transforms(u, b, m, f) + (1 + u / 10) * 1e-6 * (m - 560) ** 2


def center(summary: dict, axis: str, point: float) -> dict:
    (c,) = [c for c in summary["centers"] if c["axis"] == axis and c["point"] == point]
    return c


def interval(summary: dict, axis: str, low: float) -> dict:
    (i,) = [
        i for i in summary["intervals"] if i["axis"] == axis and i["interval"][0] == low
    ]
    return i


class Measured(unittest.TestCase):
    def test_constants_linear_in_the_axis_transforms_have_no_error(self):
        # 温度で線形に内挿すると √T の 1 次の項に誤差が出るので、燃料は √T で内挿していること、減速材は両隣の重み（2/3）で内挿して
        # いることを確かめる
        s = cgi.summarize(table(linear_in_transforms))
        self.assertEqual(
            [(c["axis"], c["point"]) for c in s["centers"]],
            [("boron", 650.0), ("moderator", 570.0), ("fuel", 600.0)],
        )
        for c in s["centers"]:
            for p in c["points"]:
                self.assertLess(p["max_abs_pcm"], 1e-9)

    def test_the_error_is_the_state_minus_its_neighbours_interpolated_by_core(self):
        s = cgi.summarize(table(moderator_square))
        c = center(s, "moderator", 570.0)
        self.assertEqual(c["neighbours"], [560.0, 575.0])
        self.assertAlmostEqual(c["weight"], 2 / 3, places=15)
        for p, u in zip(c["points"], BURNUP, strict=True):
            self.assertEqual(p["burnup"], u)
            # 状態ごとに、両隣の少数群定数を重み 2/3 で足した状態の k∞ と、格子点の状態の k∞ の Δρ
            worst = None
            for b, f in itertools.product(BORON, FUEL):
                truth = state(moderator_square(u, b, 570.0, f))
                rebuilt = unvector(
                    vector(state(moderator_square(u, b, 560.0, f))) / 3
                    + 2 * vector(state(moderator_square(u, b, 575.0, f))) / 3,
                    truth,
                )
                rho = (1 / k_inf(truth) - 1 / k_inf(rebuilt)) * 1e5
                if worst is None or abs(rho) > abs(worst[0]):
                    worst = (rho, b, f)
            self.assertAlmostEqual(p["worst_pcm"], worst[0], places=8)
            self.assertAlmostEqual(p["max_abs_pcm"], abs(worst[0]), places=8)
            self.assertEqual(
                p["worst"], {"boron": worst[1], "moderator": 570.0, "fuel": worst[2]}
            )
        # 減速材の 2 次の項は、ホウ素と燃料の格子点の内挿に誤差を出さない
        for p in (
            center(s, "fuel", 600.0)["points"] + center(s, "boron", 650.0)["points"]
        ):
            self.assertLess(p["max_abs_pcm"], 1e-9)

    def test_axes_without_an_inner_point_have_no_centers(self):
        t = table(linear_in_transforms)
        keep = [i for i, p in enumerate(t["points"]) if p[1] != 2]
        t["axes"][1]["points"] = BORON[:2]
        t["points"] = [t["points"][i] for i in keep]
        t["states"] = [t["states"][i] for i in keep]
        s = cgi.summarize(t)
        self.assertNotIn("boron", [c["axis"] for c in s["centers"]])
        self.assertNotIn("boron", [i["axis"] for i in s["intervals"]])


class Projected(unittest.TestCase):
    def test_each_interval_scales_the_measured_error_of_its_end_points(self):
        # 2 次の項 c (x − a)² の線形の内挿の誤差は、区間 [a, b] の中の x で c (x − a)(b − x) なので、両隣の間隔 h1、h2 の格子点の
        # 誤差 c h1 h2 から、幅 h の区間の中点の誤差 c h² / 4 を h² / (4 h1 h2) 倍で見込む
        s = cgi.summarize(table(moderator_square))
        measured = center(s, "moderator", 570.0)
        low = interval(s, "moderator", 560.0)
        high = interval(s, "moderator", 570.0)
        self.assertEqual(low["interval"], [560.0, 570.0])
        self.assertEqual(high["interval"], [570.0, 575.0])
        (from_low,) = low["from"]
        (from_high,) = high["from"]
        self.assertEqual(from_low["center"], 570.0)
        self.assertAlmostEqual(from_low["factor"], 10**2 / (4 * 10 * 5), places=15)
        self.assertAlmostEqual(from_high["factor"], 5**2 / (4 * 10 * 5), places=15)
        for n, p in enumerate(measured["points"]):
            self.assertAlmostEqual(
                from_low["max_abs_pcm"][n],
                p["max_abs_pcm"] * from_low["factor"],
                places=12,
            )
            self.assertEqual(low["max_abs_pcm"][n], from_low["max_abs_pcm"][n])

    def test_the_fuel_scales_in_the_square_root_of_the_temperature(self):
        s = cgi.summarize(table(moderator_square))
        (f,) = interval(s, "fuel", 600.0)["from"]
        h1, h2 = math.sqrt(600) - math.sqrt(560), math.sqrt(750) - math.sqrt(600)
        self.assertAlmostEqual(f["factor"], h2 / (4 * h1), places=12)

    def test_inner_intervals_take_the_larger_of_the_two_end_points(self):
        t = table(moderator_square)
        # 4 点の軸にして、区間 [570, 575] の両端を内側の格子点にする
        t["axes"][2]["points"] = [560.0, 570.0, 575.0, 590.0]
        t["axes"][2]["density"] = [0.75, 0.73, 0.72, 0.69]
        sizes = [len(a["points"]) for a in t["axes"]]
        t["points"] = [list(p) for p in itertools.product(*[range(n) for n in sizes])]
        t["states"] = [
            state(moderator_square(u, b, m, f))
            for u, b, m, f in itertools.product(
                BURNUP, BORON, t["axes"][2]["points"], FUEL
            )
        ]
        s = cgi.summarize(t)
        inner = interval(s, "moderator", 570.0)
        self.assertEqual([f["center"] for f in inner["from"]], [570.0, 575.0])
        for n in range(len(BURNUP)):
            self.assertEqual(
                inner["max_abs_pcm"][n], max(f["max_abs_pcm"][n] for f in inner["from"])
            )


class Inputs(unittest.TestCase):
    def test_axes_need_a_known_transform(self):
        t = table(linear_in_transforms)
        t["axes"][3]["transform"] = "log"
        with self.assertRaises(SystemExit):
            cgi.summarize(t)
        del t["axes"][3]["transform"]
        with self.assertRaises(SystemExit):
            cgi.summarize(t)

    def test_a_rod_axis_is_refused_as_table_grid_refuses_it(self):
        t = table(linear_in_transforms)
        t["axes"].append({"name": "rod", "values": ["out", "in"]})
        with self.assertRaisesRegex(SystemExit, "second table"):
            cgi.summarize(t)

    def test_the_table_must_hold_the_product_of_its_axes(self):
        t = table(linear_in_transforms)
        t["points"], t["states"] = t["points"][:-1], t["states"][:-1]
        with self.assertRaises(SystemExit):
            cgi.summarize(t)

    def test_the_command_writes_the_summary_with_its_inputs(self):
        # 入力の TOML のテーブルは TOML のディレクトリからのパスで、要約の `inputs` は TOML とテーブルのパスと md5
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "table.msgpack").write_bytes(msgpack.packb(wide_table(fuel_cubic)))
            config = root / "grid.toml"
            config.write_text(
                'table = "table.msgpack"\n'
                '[[scheme]]\naxis = "fuel"\nvariable = "sqrt"\nmethod = "lagrange"\n'
            )
            output = root / "grid.json"
            argv = sys.argv
            sys.argv = [
                "compare_grid_interpolation.py",
                "--json",
                str(output),
                str(config),
            ]
            try:
                cgi.main()
            finally:
                sys.argv = argv
            summary = json.loads(output.read_text())
            self.assertEqual(summary["format"], "xlc-compare-grid-interpolation")
            self.assertEqual(
                [i["path"] for i in summary["inputs"]],
                [str(config), str(root / "table.msgpack")],
            )
            self.assertEqual(
                [(x["axis"], x["variable"], x["method"]) for x in summary["schemes"]],
                [("fuel", "sqrt", "lagrange")],
            )
            self.assertTrue(cgi.report(summary))
            del summary["inputs"]
            with self.assertRaises(SystemExit):
                cgi.report(summary)


# 内側の格子点を除いた格子に 4 点の方式（`lagrange`）が使えるよう、どの軸も 6 点。減速材の密度は温度の線形でない
WIDE_BORON = [0.0, 325.0, 650.0, 975.0, 1300.0, 1625.0]
WIDE_MODERATOR = [560.0, 570.0, 580.0, 590.0, 600.0, 610.0]
WIDE_DENSITY = [1.2 - 1.2e-6 * t * t for t in WIDE_MODERATOR]
WIDE_FUEL = [560.0, 600.0, 750.0, 900.0, 1050.0, 1200.0]


def wide_table(g) -> dict:
    """6 点の軸の直積で、各状態の νΣf2 が g(燃焼度, ホウ素, 減速材の番号, 燃料) の状態"""
    t = {
        "format": "xlc-table",
        "split": 1,
        "axes": [
            {"name": "burnup", "points": BURNUP},
            {
                "name": "boron",
                "points": WIDE_BORON,
                "b10": 0.199,
                "nominal": 650.0,
                "transform": "linear",
            },
            {
                "name": "moderator",
                "points": WIDE_MODERATOR,
                "density": WIDE_DENSITY,
                "nominal": 580.0,
                "transform": "linear",
            },
            {
                "name": "fuel",
                "points": WIDE_FUEL,
                "nominal": 900.0,
                "transform": "sqrt",
            },
        ],
    }
    sizes = [len(a["points"]) for a in t["axes"]]
    t["points"] = [list(p) for p in itertools.product(*[range(n) for n in sizes])]
    t["states"] = [
        state(g(u, b, m, f))
        for u, b, m, f in itertools.product(
            BURNUP, WIDE_BORON, range(len(WIDE_MODERATOR)), WIDE_FUEL
        )
    ]
    return t


def boron_square(u, b, m, f) -> float:
    return 0.1 + 1e-8 * (1 + u) * (b - 300.0) ** 2


def fuel_cubic(u, b, m, f) -> float:
    r = math.sqrt(f) - 30.0
    return 0.1 + 1e-5 * (1 + u) * (r**3 + r * r)


def linear_in_density(u, b, m, f) -> float:
    return 0.1 + 0.05 * WIDE_DENSITY[m]


def linear_in_log(u, b, m, f) -> float:
    return 0.1 + 0.01 * math.log(f)


def scheme(summary: dict, axis: str, variable: str, method: str) -> dict:
    (s,) = [
        s
        for s in summary["schemes"]
        if (s["axis"], s["variable"], s["method"]) == (axis, variable, method)
    ]
    return s


def scheme_max(summary: dict, axis: str, variable: str, method: str) -> float:
    return max(
        p["max_abs_pcm"]
        for c in scheme(summary, axis, variable, method)["centers"]
        for p in c["points"]
    )


def schemes(*specs) -> list[dict]:
    return [{"axis": a, "variable": v, "method": m} for a, v, m in specs]


class Schemes(unittest.TestCase):
    def test_each_scheme_has_the_error_at_each_inner_point_and_burnup(self):
        s = cgi.summarize(
            wide_table(boron_square), schemes(("boron", "value", "hermite"))
        )
        (x,) = s["schemes"]
        self.assertEqual([c["point"] for c in x["centers"]], WIDE_BORON[1:-1])
        for c in x["centers"]:
            self.assertEqual([p["burnup"] for p in c["points"]], BURNUP)

    def test_the_linear_scheme_in_the_axis_transform_is_the_core_interpolation(self):
        specs = schemes(("moderator", "value", "linear"), ("fuel", "sqrt", "linear"))
        s = cgi.summarize(wide_table(fuel_cubic), specs)
        for spec in specs:
            centers = [c for c in s["centers"] if c["axis"] == spec["axis"]]
            x = scheme(s, spec["axis"], spec["variable"], spec["method"])
            self.assertEqual(
                [[p["worst_pcm"] for p in c["points"]] for c in x["centers"]],
                [[p["worst_pcm"] for p in c["points"]] for c in centers],
            )

    def test_hermite_reproduces_a_quadratic_in_its_variable(self):
        # 傾きを 3 点の放物線で取る 3 次 Hermite（core の燃焼度の軸。ADR-0192）は 2 次の関数を厳密に再現する
        s = cgi.summarize(
            wide_table(boron_square),
            schemes(("boron", "value", "hermite"), ("boron", "value", "linear")),
        )
        self.assertLess(scheme_max(s, "boron", "value", "hermite"), 1e-9)
        self.assertGreater(scheme_max(s, "boron", "value", "linear"), 1e-3)

    def test_lagrange_reproduces_a_cubic_in_its_variable(self):
        s = cgi.summarize(
            wide_table(fuel_cubic),
            schemes(("fuel", "sqrt", "lagrange"), ("fuel", "sqrt", "hermite")),
        )
        self.assertLess(scheme_max(s, "fuel", "sqrt", "lagrange"), 1e-9)
        self.assertGreater(scheme_max(s, "fuel", "sqrt", "hermite"), 1e-3)

    def test_spline_is_the_natural_cubic_spline_of_the_remaining_points(self):
        from scipy.interpolate import CubicSpline

        s = cgi.summarize(wide_table(fuel_cubic), schemes(("fuel", "sqrt", "spline")))
        g = [math.sqrt(f) for f in WIDE_FUEL]
        for c in scheme(s, "fuel", "sqrt", "spline")["centers"]:
            i = WIDE_FUEL.index(c["point"])
            for p, u in zip(c["points"], BURNUP, strict=True):
                worst = None
                for b, m in itertools.product(WIDE_BORON, range(len(WIDE_MODERATOR))):
                    rest = [k for k in range(len(g)) if k != i]
                    spline = CubicSpline(
                        [g[k] for k in rest],
                        [fuel_cubic(u, b, m, WIDE_FUEL[k]) for k in rest],
                        bc_type="natural",
                    )
                    truth = state(fuel_cubic(u, b, m, WIDE_FUEL[i]))
                    rebuilt = state(float(spline(g[i])))
                    rho = (1 / k_inf(truth) - 1 / k_inf(rebuilt)) * 1e5
                    if worst is None or abs(rho) > abs(worst):
                        worst = rho
                self.assertAlmostEqual(p["worst_pcm"], worst, places=8)

    def test_the_density_variable_interpolates_in_the_axis_densities(self):
        s = cgi.summarize(
            wide_table(linear_in_density),
            schemes(
                ("moderator", "density", "linear"), ("moderator", "value", "linear")
            ),
        )
        self.assertLess(scheme_max(s, "moderator", "density", "linear"), 1e-9)
        self.assertGreater(scheme_max(s, "moderator", "value", "linear"), 1e-3)

    def test_the_log_variable_interpolates_in_the_logarithm(self):
        s = cgi.summarize(
            wide_table(linear_in_log),
            schemes(("fuel", "log", "linear"), ("fuel", "sqrt", "linear")),
        )
        self.assertLess(scheme_max(s, "fuel", "log", "linear"), 1e-9)
        self.assertGreater(scheme_max(s, "fuel", "sqrt", "linear"), 1e-3)


class SchemeInputs(unittest.TestCase):
    def test_unknown_axes_variables_and_methods_are_errors(self):
        t = wide_table(boron_square)
        for spec in (
            ("rod", "value", "linear"),
            ("boron", "cube", "linear"),
            ("boron", "value", "akima"),
            ("fuel", "density", "linear"),
        ):
            with self.assertRaises(SystemExit):
                cgi.summarize(t, schemes(spec))

    def test_a_method_needs_enough_points_without_the_inner_point(self):
        # 4 点の軸では、内側の格子点を除くと 3 点で、4 点の lagrange は使えない（黙って次数を下げない）
        t = table(linear_in_transforms)
        t["axes"][2]["points"] = [560.0, 570.0, 575.0, 590.0]
        t["axes"][2]["density"] = [0.75, 0.73, 0.72, 0.69]
        sizes = [len(a["points"]) for a in t["axes"]]
        t["points"] = [list(p) for p in itertools.product(*[range(n) for n in sizes])]
        t["states"] = [
            state(linear_in_transforms(u, b, m, f))
            for u, b, m, f in itertools.product(
                BURNUP, BORON, t["axes"][2]["points"], FUEL
            )
        ]
        cgi.summarize(t, schemes(("moderator", "value", "hermite")))
        with self.assertRaises(SystemExit):
            cgi.summarize(t, schemes(("moderator", "value", "lagrange")))


if __name__ == "__main__":
    unittest.main()
