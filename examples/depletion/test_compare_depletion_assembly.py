"""compare_depletion_assembly.py の要約（集合体の燃焼の、燃焼点ごとの k の差、物質の種類ごとの核種の数密度の相対差、ピン出力の差。
ADR-0178、ADR-0181）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_depletion_assembly.py'
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_depletion_assembly as compare

NAMES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149", "B-10", "Gd-155", "Gd-157"]


# (種類, 体積, xlc の 2 点の数密度, OpenMC の 2 点の数密度)。xlc の名前は `U-235`、OpenMC の名前は `U235`
REGIONS = [
    (
        "fuel",
        1.0,
        {"U-235": [2.0, 1.0], "Pu-239": [0.0, 0.1]},
        {"U235": [2.0, 1.0], "Pu239": [0.0, 0.1]},
    ),
    (
        "fuel",
        3.0,
        {"U-235": [2.0, 1.0], "Pu-239": [0.0, 0.1]},
        {"U235": [2.0, 1.0], "Pu239": [0.0, 0.1]},
    ),
    ("ifba", 0.5, {"B-10": [2.0, 1.0]}, {"B10": [2.0, 0.9]}),
    ("gd0", 0.25, {"Gd-155": [2.0, 1.0]}, {"Gd155": [2.0, 1.0]}),
    ("gd1", 0.75, {"Gd-155": [2.0, 2.0]}, {"Gd155": [2.0, 1.0]}),
]


# xlc のセルの番号は水のギャップの 1 周を含む（3×3 の配置は 5 列）: 角 (0, 0) は 6、辺 (0, 1) は 7
CELLS = [6, 7, 6, 7, 6]
ORBITS = [[1, 1], [1, 0], [1, 1], [1, 0], [1, 1]]
PINS = [[1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0]]


def pins(edge: float, corner: float) -> list[list[float]]:
    """3×3 の配置のピン出力（辺の中央 `edge`、角 `corner`、中央は燃料のないセル）。"""
    return [[corner, edge, corner], [edge, 0.0, edge], [corner, edge, corner]]


# 最後の点のピン出力: xlc は辺 1.03、角 0.97。OpenMC は組の平均が辺 1.01、角 0.99（辺の組の中に揺らぎ）
XLC_PINS = [PINS, pins(1.03, 0.97)]
OPENMC_LAST = [[0.99, 1.02, 0.99], [1.00, 0.0, 1.01], [0.99, 1.01, 0.99]]


def std(edge: float, corner: float) -> list[list[float]]:
    """3×3 の配置の統計誤差（辺の中央 `edge`、角 `corner`）。"""
    return [[corner, edge, corner], [edge, 0.0, edge], [corner, edge, corner]]


def histories():
    """2 点（燃焼度 0、1）の、燃料 2 領域（体積 1 と 3）、IFBA 1 領域、Gd のリング 2 つの履歴。"""
    xlc = {
        "format": "xlc-depletion-history",
        "nuclides": NAMES,
        "transportable": ["U-235", "B-10"],
        "regions": [
            {"material": kind, "cell": CELLS[r], "regions": [r], "volume": volume}
            for r, (kind, volume, _, _) in enumerate(REGIONS)
        ],
        "seconds": {
            "shielding": 9.0,
            "moc": 0.5,
            "rates": 0.01,
            "matrix": 0.02,
            "cram": 0.03,
        },
        "solves": 3,
        "points": [
            {
                "burnup": burnup,
                "k_inf": k,
                "pin_powers": XLC_PINS[p],
                "densities": [
                    [values.get(n, [0.0, 0.0])[p] for n in NAMES]
                    for _, _, values, _ in REGIONS
                ],
            }
            for p, (burnup, k) in enumerate([(0.0, 1.17), (1.0, 1.12)])
        ],
    }
    omc = {
        "format": "openmc-depletion-assembly",
        "problem": "2O",
        "burnup": [0.0, 1.0],
        "k": [1.171, 1.1215],
        "k-std": [0.0002, 0.0002],
        "transportable": ["B10", "U235"],
        "transport-seconds": [100.0, 110.0, 120.0],
        "total-seconds": 331.0,
        "cram-seconds": [0.5, 0.25],
        # 3×3 の簡単な配置: 中央は燃料のないセル、辺の中央 4 つと角 4 つがそれぞれ 1 つの組
        "pins": [PINS, OPENMC_LAST],
        "pins-std": [std(0.002, 0.004), std(0.004, 0.008)],
        "regions": [
            {
                "material": kind,
                "orbit": ORBITS[r],
                "count": 1,
                "volume": volume,
                "densities": {
                    **{
                        compare.openmc_name(n): [0.0, 0.0]
                        for n in compare.nuclides_of(kind)
                    },
                    **values,
                },
            }
            for r, (kind, volume, _, values) in enumerate(REGIONS)
        ],
    }
    return xlc, omc


def with_pyrex_in_the_center(pair, center: float):
    """3×3 の配置の中央のセル（xlc の番号 12、組 [0, 0]）に、燃焼する Pyrex の領域を足す。OpenMC の最後の点の中央のピン出力は `center`。"""
    xlc, omc = pair
    xlc["regions"].append(
        {"material": "pyrex", "cell": 12, "regions": [len(REGIONS)], "volume": 0.5}
    )
    for point in xlc["points"]:
        point["densities"].append([2.0 if n == "B-10" else 0.0 for n in NAMES])
    omc["regions"].append(
        {
            "material": "pyrex",
            "orbit": [0, 0],
            "count": 1,
            "volume": 0.5,
            "densities": {"B10": [2.0, 1.0]},
        }
    )
    omc["pins"][1] = [row[:] for row in omc["pins"][1]]
    omc["pins"][1][1][1] = center
    omc["pins-std"][1][1][1] = center * 1e-3
    return xlc, omc


class Summary(unittest.TestCase):
    def test_k_differences_and_their_change_from_the_first_point(self):
        s = compare.summarize(*histories())
        first, last = s["k"]
        self.assertAlmostEqual(first["diff_pcm"], -100.0)
        self.assertAlmostEqual(first["change_pcm"], 0.0)
        self.assertAlmostEqual(last["diff_pcm"], -150.0)
        self.assertAlmostEqual(last["change_pcm"], -50.0)
        self.assertAlmostEqual(last["openmc_std_pcm"], 20.0)

    def test_densities_are_volume_weighted_means_over_the_regions_of_a_kind(self):
        s = compare.summarize(*histories())
        # 燃料: xlc と OpenMC は同じ。Pu-239 は 0 の点で OpenMC が 0 なので比べない
        self.assertAlmostEqual(s["kinds"]["fuel"]["U-235"][1], 0.0)
        self.assertIsNone(s["kinds"]["fuel"]["Pu-239"][0])
        self.assertAlmostEqual(s["kinds"]["fuel"]["Pu-239"][1], 0.0)

    def test_poison_kinds_are_compared_by_their_own_nuclides(self):
        s = compare.summarize(*histories())
        # IFBA の B-10: xlc 1.0、OpenMC 0.9（最後の点）
        self.assertAlmostEqual(s["kinds"]["ifba"]["B-10"][1], 1.0 / 0.9 - 1.0)
        self.assertEqual(set(s["kinds"]["ifba"]), {"B-10"})
        self.assertEqual(set(s["kinds"]["gd0"]), {"Gd-155", "Gd-157"})
        self.assertEqual(
            set(s["kinds"]["fuel"]),
            {"U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"},
        )

    def test_the_ring_of_gadolinium_has_its_own_difference(self):
        s = compare.summarize(*histories())
        self.assertAlmostEqual(s["kinds"]["gd0"]["Gd-155"][1], 0.0)
        # gd1 は xlc が OpenMC の 2 倍
        self.assertAlmostEqual(s["kinds"]["gd1"]["Gd-155"][1], 1.0)

    def test_the_last_point_has_the_spread_over_the_regions_of_each_kind(self):
        s = compare.summarize(*histories())
        fuel = s["spread"]["fuel"]["U-235"]
        self.assertAlmostEqual(fuel["max_abs"], 0.0)
        self.assertEqual(fuel["regions"], 2)
        self.assertAlmostEqual(s["spread"]["gd1"]["Gd-155"]["max_abs"], 1.0)
        self.assertAlmostEqual(s["spread"]["gd1"]["Gd-155"]["rms"], 1.0)

    def test_the_pin_power_error_is_that_of_the_mean_over_each_symmetry_orbit(self):
        """xlc は対称な格子を 1/8 で解くので、OpenMC のピン出力は対称操作の組で平均して比べる。組の平均の誤差は
        sqrt(Σ σ²) / 個数、平均が正の組の二乗平均平方根。"""
        s = compare.summarize(*histories())
        # 辺の組: sqrt(4 × 0.002²) / 4 = 0.002 / 2、角の組: 0.004 / 2
        self.assertAlmostEqual(
            s["pin_std_rms"][0], math.sqrt((0.001**2 + 0.002**2) / 2)
        )
        self.assertAlmostEqual(
            s["pin_std_rms"][1],
            # 最後の点の組の平均は辺 1.01、角 0.99
            math.sqrt(((0.002 / 1.01) ** 2 + (0.004 / 0.99) ** 2) / 2),
        )

    def test_pin_powers_are_compared_as_means_over_each_symmetry_orbit(self):
        """xlc の 1/8 の解と、OpenMC の組の平均を比べる。組の相対差の、核分裂のある組の二乗平均平方根と、絶対値が最大の組（位置は
        組のうち上の行から最初のセル）とその組の統計誤差。"""
        s = compare.summarize(*histories())
        first, last = s["pins"]
        self.assertEqual(first["orbits"], 2)
        self.assertAlmostEqual(first["rms"], 0.0)
        edge, corner = 1.03 / 1.01 - 1.0, 0.97 / 0.99 - 1.0
        self.assertAlmostEqual(last["rms"], math.sqrt((edge**2 + corner**2) / 2))
        self.assertAlmostEqual(last["max"], corner)
        self.assertEqual(last["at"], [0, 0])
        # 角の組の統計誤差: sqrt(4 × 0.008²) / 4 を組の平均 0.99 で割る
        self.assertAlmostEqual(last["max_std"], 0.004 / 0.99)

    def test_pins_with_fission_in_other_cells_stop_the_comparison(self):
        xlc, omc = histories()
        xlc["points"][1]["pin_powers"] = [
            [0.97, 1.03, 0.97],
            [1.03, 0.1, 1.03],
            [0.97, 1.03, 0.97],
        ]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_fission_of_a_burnt_absorber_in_openmc_is_left_out_of_the_pin_comparison(
        self,
    ):
        """燃焼する吸収体（Pyrex、WABA）だけのセルは核分裂のないセルで、比べる組に入れない。OpenMC の燃焼はチェーンの全核種を
        CRAM で解くので、吸収体にも丸めの誤差の大きさのアクチニドができ、核分裂の tally が 1e-24〜5e-16 ほどの値を持つ。その最大は要約に残す。"""
        xlc, omc = with_pyrex_in_the_center(histories(), 7e-24)
        s = compare.summarize(xlc, omc)
        expected = compare.summarize(*histories())
        self.assertEqual(s["pins"][1]["orbits"], 2)
        self.assertAlmostEqual(s["pins"][1]["rms"], expected["pins"][1]["rms"])
        self.assertAlmostEqual(s["pin_std_rms"][1], expected["pin_std_rms"][1])
        self.assertEqual(s["pins"][0]["absorber_max"], 0.0)
        self.assertEqual(s["pins"][1]["absorber_max"], 7e-24)

    def test_pin_powers_are_compared_relative_to_the_mean_over_the_fuel_pins(self):
        """OpenMC の抽出は核分裂の値が正のセルの平均で割るので、吸収体のセルに丸めの誤差の核分裂があると、燃料のピンの平均でなく
        吸収体のセルも数えた平均で割られる（3×3 では 9/8 倍）。比べる前に、燃料のピンの平均で割り直す。"""
        xlc, omc = with_pyrex_in_the_center(histories(), 7e-24)
        omc["pins"][1] = [[v * 9 / 8 for v in row] for row in omc["pins"][1]]
        omc["pins-std"][1] = [[v * 9 / 8 for v in row] for row in omc["pins-std"][1]]
        s = compare.summarize(xlc, omc)
        expected = compare.summarize(*histories())
        self.assertAlmostEqual(s["pins"][1]["rms"], expected["pins"][1]["rms"])
        self.assertAlmostEqual(s["pins"][1]["max"], expected["pins"][1]["max"])
        self.assertAlmostEqual(s["pins"][1]["max_std"], expected["pins"][1]["max_std"])
        self.assertAlmostEqual(s["pin_std_rms"][1], expected["pin_std_rms"][1])
        # 燃料のピンの平均（9/8）で割り直した値
        self.assertAlmostEqual(s["pins"][1]["absorber_max"] / 7e-24, 1.0)

    def test_openmc_fission_in_a_cell_without_burnable_material_stops_the_comparison(
        self,
    ):
        """燃焼する物質のない（燃焼しない水や構造材だけの）セルに核分裂があれば、配置の取り違えなので止める。"""
        xlc, omc = histories()
        omc["pins"][1] = [[0.99, 1.02, 0.99], [1.00, 7e-24, 1.01], [0.99, 1.01, 0.99]]
        with self.assertRaisesRegex(SystemExit, "fission"):
            compare.summarize(xlc, omc)

    def test_fuel_sized_fission_in_a_burnt_absorber_stops_the_comparison(self):
        """丸めの誤差の桁でない吸収体の核分裂は、配置の取り違えとして止める。"""
        xlc, omc = with_pyrex_in_the_center(histories(), 0.5)
        with self.assertRaisesRegex(SystemExit, "round-off"):
            compare.summarize(xlc, omc)

    def test_negative_round_off_fission_of_a_burnt_absorber_counts_by_its_size(self):
        xlc, omc = with_pyrex_in_the_center(histories(), -3e-15)
        s = compare.summarize(xlc, omc)
        self.assertAlmostEqual(s["pins"][1]["absorber_max"] / 3e-15, 1.0)

    def test_burnable_materials_without_fuel_stop_the_comparison(self):
        with self.assertRaisesRegex(SystemExit, "no fuel"):
            compare.cell_kinds([{"material": "pyrex", "orbit": [0, 0]}])

    def test_xlc_fission_in_a_burnt_absorber_stops_the_comparison(self):
        xlc, omc = with_pyrex_in_the_center(histories(), 0.0)
        xlc["points"][1]["pin_powers"] = [
            [0.97, 1.03, 0.97],
            [1.03, 0.1, 1.03],
            [0.97, 1.03, 0.97],
        ]
        with self.assertRaisesRegex(SystemExit, "fission"):
            compare.summarize(xlc, omc)

    def test_xlc_pins_that_differ_within_an_orbit_stop_the_comparison(self):
        """xlc は対称操作の組の中で同じ値を持つはずで、違えば（対称でない格子や入力の誤り）平均で隠さず止める。"""
        xlc, omc = histories()
        xlc["points"][1]["pin_powers"] = [
            [0.97, 1.03, 0.97],
            [1.03, 0.0, 1.031],
            [0.97, 1.03, 0.97],
        ]
        with self.assertRaisesRegex(SystemExit, "orbit"):
            compare.summarize(xlc, omc)

    def test_a_history_without_pin_powers_stops_the_comparison(self):
        xlc, omc = histories()
        del xlc["points"][0]["pin_powers"]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_pin_maps_of_other_shapes_stop_the_comparison(self):
        xlc, omc = histories()
        xlc["points"][0]["pin_powers"] = [[1.0, 1.0], [1.0, 1.0]]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_pin_errors_of_another_shape_stop_the_comparison(self):
        xlc, omc = histories()
        omc["pins-std"][1] = [[0.01, 0.01], [0.01, 0.01]]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_the_report_of_a_summary_without_pins_asks_to_rerun_the_stage(self):
        summary = {"inputs": [], **compare.summarize(*histories())}
        del summary["pins"]
        with self.assertRaisesRegex(SystemExit, "rerun"):
            compare.report(summary)

    def test_the_report_of_fewer_pin_rows_than_points_stops(self):
        summary = {"inputs": [], **compare.summarize(*histories())}
        summary["pins"] = summary["pins"][:1]
        with self.assertRaisesRegex(SystemExit, "pin"):
            compare.report(summary)

    def test_the_report_has_the_pin_power_differences(self):
        summary = {"inputs": [], **compare.summarize(*histories())}
        lines = compare.report(summary)
        header = lines.index(
            "ピン出力の差（対称操作の組、%）: 燃焼度 | 組 | rms | 最大 ± OpenMC の統計誤差 | 位置（行, 列。1 から） | OpenMC の統計誤差 rms"
        )
        self.assertEqual(
            lines[header + 2],
            f"1 | 2 | {100 * math.sqrt(((1.03 / 1.01 - 1) ** 2 + (0.97 / 0.99 - 1) ** 2) / 2):.3f} | "
            f"{100 * (0.97 / 0.99 - 1):+.3f} ± {100 * 0.004 / 0.99:.3f} | (1, 1) | "
            f"{100 * math.sqrt(((0.002 / 1.01) ** 2 + (0.004 / 0.99) ** 2) / 2):.3f}",
        )

    def test_the_report_has_the_largest_fission_of_burnt_absorbers(self):
        xlc, omc = with_pyrex_in_the_center(histories(), 7e-24)
        summary = {"inputs": [], **compare.summarize(xlc, omc)}
        self.assertIn(
            "燃焼する吸収体だけのセルの OpenMC の核分裂（ピン出力。比べない）の最大: 7.0e-24",
            compare.report(summary),
        )

    def test_times_of_both_codes(self):
        s = compare.summarize(*histories())
        self.assertEqual(s["seconds"]["xlc"]["solves"], 3)
        o = s["seconds"]["openmc"]
        self.assertEqual(
            (o["transport"], o["fastest"], o["slowest"]), (330.0, 100.0, 120.0)
        )

    def test_different_transported_nuclides_stop_the_comparison(self):
        xlc, omc = histories()
        omc["transportable"] = ["U235"]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_different_burnup_points_stop_the_comparison(self):
        xlc, omc = histories()
        omc["burnup"] = [0.0, 2.0]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_a_region_in_another_orbit_stops_the_comparison(self):
        """体積と種類が同じでも、組が違う領域（軸上と対角線上は同じセル数）を取り違えない。"""
        xlc, omc = histories()
        omc["regions"][0]["orbit"] = [1, 0]
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_regions_must_agree_in_kind_and_volume(self):
        xlc, omc = histories()
        omc["regions"][2]["material"] = "pyrex"
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)
        xlc, omc = histories()
        omc["regions"][1]["volume"] = 3.1
        with self.assertRaises(SystemExit):
            compare.summarize(xlc, omc)

    def test_the_names_of_xlc_nuclides_become_openmc_names(self):
        self.assertEqual(compare.openmc_name("Am-242m"), "Am242_m1")
        self.assertEqual(compare.openmc_name("Gd-155"), "Gd155")
        self.assertTrue(math.isclose(compare.relative(1.1, 1.0), 0.1))


if __name__ == "__main__":
    unittest.main()
