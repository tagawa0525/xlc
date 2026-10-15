"""openmc_depletion_assembly.py の、集合体の問題の定義とモデルの組み立て（ADR-0178）。

OpenMC のある環境で実行する（モデルを組み立てるだけで、輸送は解かない）:
nix develop .#references -c python3 -m unittest discover -s examples/depletion -p 'test_openmc_depletion_assembly.py'
"""

import itertools
import math
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import openmc_depletion_assembly as asm

PELLET_AREA = math.pi * 0.4096**2


def run_of(problem: str, **extra) -> dict:
    """入力を読んだ形（読み込みを通さず、モデルの組み立てに要るものだけ）。"""
    return {"problem": problem, **extra}


def regions_of(problem: str, **extra) -> list[asm.Region]:
    return asm.regions(run_of(problem, **extra))


def pins_of(regions: list[asm.Region], kind: str) -> int:
    """`kind` の領域が表すセルの数（組の大きさの和）。"""
    return sum(r.count for r in regions if r.kind == kind)


class Orbits(unittest.TestCase):
    def test_octant_positions_cover_the_17_by_17_lattice(self):
        sizes = [asm.orbit_size((a, b)) for a in range(9) for b in range(a + 1)]
        self.assertEqual(sum(sizes), 289)

    def test_a_position_belongs_to_the_orbit_of_its_distances_from_the_centre(self):
        self.assertEqual(asm.orbit(8, 8), (0, 0))
        self.assertEqual(asm.orbit(0, 0), (8, 8))
        self.assertEqual(asm.orbit(11, 5), (3, 3))
        self.assertEqual(asm.orbit(5, 14), (6, 3))
        self.assertEqual(asm.orbit(14, 5), (6, 3))

    def test_orbit_sizes_are_1_4_or_8(self):
        self.assertEqual(asm.orbit_size((0, 0)), 1)
        self.assertEqual(asm.orbit_size((3, 0)), 4)
        self.assertEqual(asm.orbit_size((3, 3)), 4)
        self.assertEqual(asm.orbit_size((3, 1)), 8)


class Cases(unittest.TestCase):
    def kinds(self, problem: str, **extra) -> dict[str, int]:
        found: dict[str, int] = {}
        for i in range(17):
            for j in range(17):
                kind = asm.cell_kind(problem, i, j)
                found[kind] = found.get(kind, 0) + 1
        return found

    def test_2a_has_264_fuel_pins_24_guide_tubes_and_an_instrument_tube(self):
        self.assertEqual(self.kinds("2A"), {"fuel": 264, "gt": 24, "it": 1})

    def test_the_history_cases_change_only_what_they_name(self):
        self.assertEqual(self.kinds("2A-density"), self.kinds("2A"))
        self.assertEqual(self.kinds("2A-rods"), {"fuel": 264, "aic": 24, "it": 1})
        self.assertEqual(asm.CASES["2A-density"].temperature, 565.0)
        self.assertEqual(asm.CASES["2A-rods"].temperature, 565.0)
        self.assertEqual(asm.CASES["2A-density"].moderator, asm.MODERATOR)
        self.assertEqual(asm.CASES["2A-rods"].moderator, asm.MODERATOR_DENSE)

    def test_2f_puts_pyrex_in_all_24_guide_tubes_at_600_k(self):
        self.assertEqual(self.kinds("2F"), {"fuel": 264, "pyrex": 24, "it": 1})
        self.assertEqual(asm.CASES["2F"].temperature, 600.0)

    def test_2m_has_128_ifba_rods(self):
        self.assertEqual(
            self.kinds("2M"), {"fuel": 136, "ifba": 128, "gt": 24, "it": 1}
        )

    def test_2n_has_104_ifba_rods_and_20_waba_in_the_guide_tubes(self):
        self.assertEqual(
            self.kinds("2N"),
            {"fuel": 160, "ifba": 104, "waba": 20, "gt": 4, "it": 1},
        )

    def test_2o_and_2p_have_12_and_24_gadolinia_rods(self):
        self.assertEqual(self.kinds("2O")["gd"], 12)
        self.assertEqual(self.kinds("2P")["gd"], 24)
        self.assertEqual(self.kinds("2O")["fuel"], 252)
        self.assertEqual(self.kinds("2P")["fuel"], 240)

    def test_every_layout_is_symmetric_under_the_eight_operations(self):
        for problem in asm.CASES:
            for i in range(17):
                for j in range(17):
                    a, b = asm.orbit(i, j)
                    # 組の代表（中央から (a, b)）と同じ種類
                    self.assertEqual(
                        asm.cell_kind(problem, i, j),
                        asm.cell_kind(problem, 8 + a, 8 + b),
                        f"{problem} {(i, j)}",
                    )

    def test_an_unknown_problem_stops(self):
        with self.assertRaises(SystemExit):
            asm.cell_kind("2L", 0, 0)


class Burnup(unittest.TestCase):
    def test_the_two_grids_have_40_and_70_steps_up_to_60_mwd_per_kg(self):
        self.assertEqual(len(asm.BURNUP_40), 40)
        self.assertEqual(len(asm.BURNUP_70), 70)
        self.assertEqual(asm.BURNUP_40[-1], 60.0)
        self.assertEqual(asm.BURNUP_70[-1], 60.0)

    def test_the_grids_increase_from_a_positive_first_point(self):
        for grid in (asm.BURNUP_40, asm.BURNUP_70):
            self.assertGreater(grid[0], 0)
            self.assertTrue(all(b > a for a, b in itertools.pairwise(grid)))

    def test_the_70_step_grid_is_fine_where_the_poisons_burn(self):
        """IFBA の B-10 と Gd が燃える 16 MWd/kg までは、40 点の格子の 2 倍近い点を持つ。"""
        under = lambda grid: sum(1 for b in grid if b <= 16)
        self.assertGreaterEqual(under(asm.BURNUP_70), 1.7 * under(asm.BURNUP_40))


class Regions(unittest.TestCase):
    def test_2a_has_39_fuel_regions_of_the_octant_positions(self):
        regions = regions_of("2A")
        self.assertEqual(len(regions), 39)
        self.assertEqual({r.kind for r in regions}, {"fuel"})
        self.assertEqual(pins_of(regions, "fuel"), 264)

    def test_the_volume_is_the_pellet_area_times_the_pins_of_the_orbit(self):
        for r in regions_of("2A"):
            self.assertAlmostEqual(r.volume, r.count * PELLET_AREA, places=12)

    def test_regions_are_in_the_order_of_first_appearance_scanning_rows_from_the_top(
        self,
    ):
        regions = regions_of("2A")
        # 左上の角のピン（行 0、列 0）の組が最初。次は同じ行の隣
        self.assertEqual(regions[0].orbit, (8, 8))
        self.assertEqual(regions[1].orbit, (8, 7))
        self.assertEqual(regions[2].orbit, (8, 6))

    def test_ifba_has_its_own_region_of_the_coating_annulus(self):
        regions = regions_of("2M")
        coating = math.pi * (0.4106**2 - 0.4096**2)
        self.assertEqual(pins_of(regions, "ifba"), 128)
        for r in regions:
            if r.kind == "ifba":
                self.assertAlmostEqual(r.volume, r.count * coating, places=12)
        self.assertEqual(pins_of(regions, "fuel"), 264)

    def test_pyrex_and_waba_are_regions_of_their_annuli(self):
        pyrex = regions_of("2F")
        self.assertEqual(pins_of(pyrex, "pyrex"), 24)
        for r in pyrex:
            if r.kind == "pyrex":
                self.assertAlmostEqual(
                    r.volume, r.count * math.pi * (0.427**2 - 0.241**2), places=12
                )
        waba = regions_of("2N")
        self.assertEqual(pins_of(waba, "waba"), 20)
        self.assertEqual(pins_of(waba, "ifba"), 104)

    def test_the_rods_of_the_history_case_do_not_burn(self):
        self.assertEqual({r.kind for r in regions_of("2A-rods")}, {"fuel"})

    def test_gadolinia_pins_burn_in_equal_area_rings_named_gd0_to_gd_n(self):
        regions = regions_of("2O", **{"gd-rings": 8})
        kinds = [f"gd{k}" for k in range(8)]
        for kind in kinds:
            self.assertEqual(pins_of(regions, kind), 12)
        self.assertEqual({r.kind for r in regions}, {"fuel", *kinds})
        by_orbit: dict[tuple, set[float]] = {}
        for r in regions:
            if r.kind in kinds:
                by_orbit.setdefault(r.orbit, set()).add(round(r.volume / r.count, 12))
        for orbit, areas in by_orbit.items():
            self.assertEqual(areas, {round(PELLET_AREA / 8, 12)}, orbit)
        # 組のリングは物質の順（内側から）に並ぶ
        first = [
            r.kind for r in regions if r.kind.startswith("gd") and r.orbit == (4, 2)
        ]
        self.assertEqual(first, kinds)

    def test_ring_radii_split_the_pellet_into_equal_areas(self):
        radii = asm.gd_radii(4)
        self.assertEqual(radii[-1], 0.4096)
        areas = [math.pi * (b * b - a * a) for a, b in zip([0.0, *radii], radii)]
        for area in areas:
            self.assertAlmostEqual(area, PELLET_AREA / 4, places=12)


class Reading(unittest.TestCase):
    def write(self, text: str) -> Path:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "openmc-2A.toml"
        path.write_text(text)
        return path

    BASE = """
problem = "{problem}"
cross-sections = "openmc-library/cross_sections.xml"
chain = "../openmc-chain.xml"
power = 38.41
burnup = [0.1, 0.5]
particles = 1000
batches = 20
inactive = 5
seed = 1
threads = 2
"""

    def test_reads_an_input_and_makes_the_paths_absolute(self):
        run = asm.read(self.write(self.BASE.format(problem="2A")))
        self.assertEqual(run["problem"], "2A")
        self.assertTrue(Path(run["cross-sections"]).is_absolute())
        self.assertTrue(Path(run["chain"]).is_absolute())

    def test_a_gadolinia_problem_needs_its_ring_count_and_others_refuse_it(self):
        with self.assertRaises(SystemExit):
            asm.read(self.write(self.BASE.format(problem="2O")))
        run = asm.read(self.write(self.BASE.format(problem="2O") + "gd-rings = 8\n"))
        self.assertEqual(run["gd-rings"], 8)
        with self.assertRaises(SystemExit):
            asm.read(self.write(self.BASE.format(problem="2A") + "gd-rings = 8\n"))

    def test_unknown_keys_and_problems_stop(self):
        with self.assertRaises(SystemExit):
            asm.read(self.write(self.BASE.format(problem="2A") + "typo = 1\n"))
        with self.assertRaises(SystemExit):
            asm.read(self.write(self.BASE.format(problem="1C")))

    def test_start_from_needs_a_non_negative_start_point(self):
        base = self.BASE.format(problem="2A") + 'start-from = "x.msgpack"\n'
        with self.assertRaises(SystemExit):
            asm.read(self.write(base))
        with self.assertRaises(SystemExit):
            asm.read(self.write(base + "start-point = -1\n"))

    def test_the_burnup_points_must_increase_from_a_positive_value(self):
        for burnup in ("[0.0, 1.0]", "[1.0, 0.5]", "[]"):
            text = self.BASE.format(problem="2A").replace("[0.1, 0.5]", burnup)
            with self.assertRaises(SystemExit):
                asm.read(self.write(text))


class Model(unittest.TestCase):
    def test_only_the_burnable_materials_are_depletable_and_all_have_volumes(self):
        model, regions = asm.build(
            run_of("2N", particles=1000, batches=20, inactive=5, seed=1)
        )
        depletable = [m for m in model.materials if m.depletable]
        self.assertEqual(len(depletable), len(regions))
        for m in depletable:
            self.assertGreater(m.volume, 0)
        names = {m.name for m in model.materials if not m.depletable}
        self.assertEqual(names, {"moderator", "zircaloy-4"})

    def test_the_rods_of_2a_rods_are_a_fixed_material(self):
        model, _ = asm.build(
            run_of("2A-rods", particles=1000, batches=20, inactive=5, seed=1)
        )
        fixed = {m.name for m in model.materials if not m.depletable}
        self.assertIn("aic", fixed)
        self.assertIn("ss304", fixed)

    def test_the_volume_of_a_region_is_the_sum_over_its_cells(self):
        model, _ = asm.build(
            run_of("2A", particles=1000, batches=20, inactive=5, seed=1)
        )
        self.assertAlmostEqual(
            sum(m.volume for m in model.materials if m.depletable),
            264 * PELLET_AREA,
            places=9,
        )

    def test_a_mesh_tally_scores_the_fission_of_every_pin(self):
        model, _ = asm.build(
            run_of("2A", particles=1000, batches=20, inactive=5, seed=1)
        )
        (tally,) = model.tallies
        self.assertEqual(tally.name, "fission")
        self.assertEqual(tally.scores, ["fission"])

    def test_the_moderator_density_of_the_history_case_is_the_spec_0_661(self):
        model, _ = asm.build(
            run_of("2A-density", particles=1000, batches=20, inactive=5, seed=1)
        )
        water = next(m for m in model.materials if m.name == "moderator")
        h1 = water.get_nuclide_atom_densities()["H1"]
        self.assertAlmostEqual(h1, 4.41459e-02)

    def test_a_565_k_problem_is_solved_with_interpolated_temperatures(self):
        model, _ = asm.build(
            run_of("2A", particles=1000, batches=20, inactive=5, seed=1)
        )
        self.assertEqual(model.settings.temperature["method"], "interpolation")
        model, _ = asm.build(
            run_of("2F", particles=1000, batches=20, inactive=5, seed=1)
        )
        self.assertEqual(model.settings.temperature["method"], "nearest")


class Compositions(unittest.TestCase):
    def test_ifba_has_the_zirconium_diboride_of_the_spec(self):
        # 表 P2-4: B-10 2.16410E-02、B-11 1.96824E-02、Zr 2.06617E-02（表 12 の合計）
        self.assertAlmostEqual(asm.IFBA["B10"], 2.16410e-2)
        self.assertAlmostEqual(asm.IFBA["B11"], 1.96824e-2)
        zr = sum(v for k, v in asm.IFBA.items() if k.startswith("Zr"))
        self.assertAlmostEqual(zr, 2.06617e-2, delta=1e-6)

    def test_gadolinia_fuel_has_1_8_percent_uranium_and_the_seven_gadolinium_isotopes(
        self,
    ):
        self.assertAlmostEqual(asm.GD_FUEL["U235"], 3.90500e-4)
        gd = {k for k in asm.GD_FUEL if k.startswith("Gd")}
        self.assertEqual(gd, {f"Gd{a}" for a in (152, 154, 155, 156, 157, 158, 160)})

    def test_waba_splits_natural_carbon_into_c12_and_c13(self):
        self.assertNotIn("C", asm.WABA)
        self.assertAlmostEqual(asm.WABA["C12"] + asm.WABA["C13"], 3.77001e-3, places=12)
        self.assertAlmostEqual(asm.WABA["Al27"], 3.90223e-2)

    def test_waba_radii_are_those_of_table_13(self):
        self.assertEqual(asm.WABA_RADII, (0.286, 0.339, 0.353, 0.404, 0.418, 0.484))


class IntervalPins(unittest.TestCase):
    """区間に分けた実行は、区間ごとに前の区間までのピン出力に続けて pins.json を書く（再開した VM に前の区間はない）。"""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.work = Path(tmp.name)
        saved = asm.bos_points
        self.addCleanup(setattr, asm, "bos_points", saved)

    def point(self, value: float) -> dict:
        return {
            "mean": [[value]],
            "std": [[0.0]],
            "k": value,
            "k-std": 0.0,
            "runtime": {},
        }

    def interval(self, k: int) -> Path:
        step = self.work / f"step-{k}"
        step.mkdir()
        return step

    def test_the_first_interval_writes_its_own_points(self):
        asm.bos_points = lambda directory: [self.point(0.0)]
        step = self.interval(1)
        asm.collect_interval(step, 1, None)
        self.assertEqual([p["k"] for p in asm.read_pins(step)], [0.0])

    def test_an_interval_continues_the_previous_points(self):
        asm.bos_points = lambda directory: [self.point(0.0)]
        first = self.interval(1)
        asm.collect_interval(first, 1, None)
        asm.bos_points = lambda directory: [self.point(1.0)]
        second = self.interval(2)
        asm.collect_interval(second, 2, first)
        self.assertEqual([p["k"] for p in asm.read_pins(second)], [0.0, 1.0])

    def test_the_last_interval_adds_the_final_point_too(self):
        asm.bos_points = lambda directory: [self.point(0.0)]
        first = self.interval(1)
        asm.collect_interval(first, 1, None)
        asm.bos_points = lambda directory: [self.point(1.0), self.point(2.0)]
        second = self.interval(2)
        asm.collect_interval(second, 2, first)
        self.assertEqual([p["k"] for p in asm.read_pins(second)], [0.0, 1.0, 2.0])

    def test_a_rerun_interval_replaces_the_points_it_made_before(self):
        """前の区間の pins.json が後の点まで持っていても（再開で区間を再実行する）、区間 k の分で置き換える。"""
        first = self.interval(1)
        asm.write_pins(first, [self.point(0.0), self.point(9.0)])
        asm.bos_points = lambda directory: [self.point(1.0)]
        second = self.interval(2)
        asm.collect_interval(second, 2, first)
        self.assertEqual([p["k"] for p in asm.read_pins(second)], [0.0, 1.0])


if __name__ == "__main__":
    unittest.main()
