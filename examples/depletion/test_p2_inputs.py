"""p2_inputs.py が書く、集合体の燃焼の入力（OpenMC と xlc。ADR-0178）の性質。

OpenMC のある環境で実行する: nix develop .#references -c python3 -m unittest discover -s examples/depletion -p 'test_p2_inputs.py'
"""

import re
import sys
import unittest
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import openmc_depletion_assembly as asm
import p2_inputs

P2 = HERE / "p2"


def load(path: Path) -> dict:
    return tomllib.loads(path.read_text())


def lattice(job: str) -> dict:
    return load(P2 / f"lattice-{job}.toml")


def openmc_input(job: str) -> dict:
    return load(P2 / f"openmc-{job}.toml")


class WrittenInputs(unittest.TestCase):
    def test_the_inputs_on_disk_are_what_the_tables_make(self):
        for path, text in p2_inputs.files().items():
            self.assertTrue(path.exists(), path)
            self.assertEqual(path.read_text(), text, f"{path} is stale")

    def test_there_are_six_cases_and_two_history_runs(self):
        self.assertEqual(
            set(p2_inputs.JOBS),
            {"2A", "2F", "2M", "2N", "2O", "2P", "2A-density", "2A-rods"},
        )

    def test_burnup_points_and_power_are_the_same_in_both_codes(self):
        for job in p2_inputs.JOBS:
            omc, xlc = openmc_input(job), lattice(job)
            self.assertEqual(omc["burnup"], xlc["depletion"]["burnup"], job)
            self.assertEqual(omc["power"], xlc["depletion"]["power"], job)

    def test_the_step_counts_follow_the_plan_table(self):
        counts = {job: len(openmc_input(job)["burnup"]) for job in p2_inputs.JOBS}
        self.assertEqual(
            counts,
            {
                "2A": 40,
                "2F": 40,
                "2N": 40,
                "2A-density": 40,
                "2A-rods": 40,
                "2M": 70,
                "2O": 70,
                "2P": 70,
            },
        )

    def test_every_openmc_input_is_one_the_tool_reads(self):
        for job in p2_inputs.JOBS:
            run = asm.read(P2 / f"openmc-{job}.toml")
            self.assertEqual(run["problem"], p2_inputs.JOBS[job][0])
            self.assertEqual(run["particles"], p2_inputs.PARTICLES)

    def test_job_names_match_the_aws_definition(self):
        text = (P2 / "aws.toml").read_text()
        pattern = re.search(r"^jobs = '([^']*)'", text, re.MULTILINE).group(1)
        for job in p2_inputs.JOBS:
            self.assertRegex(f"openmc-{job}", pattern)


class TimeFactorInputs(unittest.TestCase):
    """時間の倍率 s を測る短い実行の入力（2A、粒子 2 万・世代 40、手元の Slurm で `-c 8`）。"""

    def test_the_two_inputs_differ_only_in_the_starting_composition(self):
        fresh = load(P2 / "s-2A-fresh.toml")
        burned = load(P2 / "s-2A-burned.toml")
        self.assertNotIn("start-from", fresh)
        self.assertEqual(burned["start-point"], 6)
        rest = {
            k: v for k, v in burned.items() if k not in ("start-from", "start-point")
        }
        self.assertEqual(rest, fresh)

    def test_they_run_one_short_step_of_2a_with_20000_particles_and_40_batches(self):
        run = asm.read(P2 / "s-2A-burned.toml")
        self.assertEqual(run["problem"], "2A")
        self.assertEqual((run["particles"], run["batches"]), (20000, 40))
        self.assertEqual(run["burnup"], [0.1])
        self.assertEqual(run["threads"], 8)

    def test_the_composition_is_the_5_mwd_per_kg_point_of_the_1c_burn(self):
        run = asm.read(P2 / "s-2A-burned.toml")
        self.assertEqual(run["start-from"].name, "openmc-1C-rings1.msgpack")
        burnup = [0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
        self.assertEqual(burnup[load(P2 / "s-2A-burned.toml")["start-point"]], 5.0)


class Lattice(unittest.TestCase):
    def test_rows_are_17_by_17_cells_that_are_defined(self):
        for job in p2_inputs.JOBS:
            xlc = lattice(job)
            rows = [row.split() for row in xlc["rows"]]
            self.assertEqual([len(r) for r in rows], [17] * 17, job)
            self.assertLessEqual({c for r in rows for c in r}, set(xlc["cell"]), job)

    def test_the_cells_of_the_rows_are_the_kinds_of_the_openmc_layout(self):
        for job, (problem, _, _) in p2_inputs.JOBS.items():
            rows = [row.split() for row in lattice(job)["rows"]]
            for i in range(17):
                for j in range(17):
                    self.assertEqual(
                        rows[i][j],
                        p2_inputs.CELL_NAMES[asm.cell_kind(problem, i, j)],
                        f"{job} {(i, j)}",
                    )

    def test_zone_materials_are_defined_and_radii_increase(self):
        for job in p2_inputs.JOBS:
            xlc = lattice(job)
            for name, cell in xlc["cell"].items():
                radii = [z["radius"] for z in cell["zones"] if "radius" in z]
                self.assertEqual(radii, sorted(radii), f"{job} {name}")
                for z in cell["zones"]:
                    self.assertTrue(
                        z["material"] == "void" or z["material"] in xlc["material"],
                        f"{job} {name} {z['material']}",
                    )

    def test_the_depleted_materials_are_the_regions_of_the_openmc_run(self):
        for job, (problem, _, _) in p2_inputs.JOBS.items():
            run = {
                "problem": problem,
                **({"gd-rings": 8} if asm.CASES[problem].gd else {}),
            }
            kinds = list(dict.fromkeys(r.kind for r in asm.regions(run)))
            self.assertEqual(lattice(job)["depletion"]["materials"], kinds, job)

    def test_gadolinium_rings_are_equal_area_rings_of_the_pellet(self):
        for job in ("2O", "2P"):
            cell = lattice(job)["cell"]["gd"]
            radii = [
                z["radius"] for z in cell["zones"] if z["material"].startswith("gd")
            ]
            self.assertEqual(len(radii), p2_inputs.GD_RINGS)
            for got, want in zip(radii, asm.gd_radii(p2_inputs.GD_RINGS)):
                self.assertAlmostEqual(got, want, places=9)

    def test_every_input_uses_the_default_track_table(self):
        """Gd のリングごとの物質も飛程の表の 1 つの軸を共有して解ける（ADR-0182）ので、どの入力も既定の表を使う。"""
        for job in p2_inputs.JOBS:
            self.assertNotIn("track-table", lattice(job), job)

    def test_every_material_has_the_temperature_of_the_case(self):
        for job, (problem, _, _) in p2_inputs.JOBS.items():
            for name, m in lattice(job)["material"].items():
                self.assertEqual(
                    m["temperature"], asm.CASES[problem].temperature, f"{job} {name}"
                )

    def test_the_history_cases_change_one_condition_of_2a(self):
        base = lattice("2A")
        self.assertEqual(lattice("2A-density")["rows"], base["rows"])
        self.assertNotEqual(
            lattice("2A-density")["material"]["water"], base["material"]["water"]
        )
        rods = lattice("2A-rods")
        self.assertEqual(rods["material"]["water"], base["material"]["water"])
        self.assertEqual(rods["depletion"]["materials"], ["fuel"])
        self.assertIn("aic", rods["material"])


class Libraries(unittest.TestCase):
    """入力の核種は、ライブラリ（p2/xs.toml と p2/openmc-library.toml）が入力の温度で持つ。"""

    @classmethod
    def setUpClass(cls):
        cls.xs = load(P2 / "xs.toml")["nuclide"]
        cls.omc = load(P2 / "openmc-library.toml")["neutron"]

    def covered(self, table: dict, name: str, temperature: float) -> bool:
        if name not in table:
            return False
        temperatures = table[name]["temperatures"]
        if name in ("H-1", "H1"):
            # S(α,β) の温度の内挿（ADR-0080）。範囲の中にあればよい
            return min(temperatures) <= temperature <= max(temperatures)
        return temperature in temperatures

    def test_the_xlc_library_has_every_nuclide_of_the_lattice_inputs_at_their_temperature(
        self,
    ):
        for job in p2_inputs.JOBS:
            for name, m in lattice(job)["material"].items():
                nuclides = m.get("nuclides", {})
                for nuclide in nuclides:
                    self.assertTrue(
                        self.covered(self.xs, nuclide, m["temperature"]),
                        f"{job} {name} {nuclide}",
                    )

    def test_the_openmc_library_has_every_nuclide_of_the_models_at_their_temperature(
        self,
    ):
        for job, (problem, _, _) in p2_inputs.JOBS.items():
            run = {
                "problem": problem,
                "particles": 1,
                "batches": 2,
                "inactive": 1,
                "seed": 1,
            }
            if asm.CASES[problem].gd:
                run["gd-rings"] = 2
            model, _ = asm.build(run)
            for m in model.materials:
                for nuclide in m.get_nuclide_atom_densities():
                    self.assertTrue(
                        self.covered(self.omc, nuclide, m.temperature),
                        f"{job} {m.name} {nuclide}",
                    )

    def test_every_nuclide_has_two_temperatures_for_the_interpolation_mode(self):
        for name, item in self.omc.items():
            self.assertGreaterEqual(len(item["temperatures"]), 2, name)


if __name__ == "__main__":
    unittest.main()
