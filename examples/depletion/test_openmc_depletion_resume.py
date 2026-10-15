"""openmc_depletion_resume.py の区間（ステップ）の置き場の規則（ADR-0174）。

OpenMC を使わない部分だけを確かめる。区間に分けた実行が分けない実行と同じ結果を出すことは、OpenMC のある環境で
小さい入力を流して確かめる（ADR-0174 の「確認」）。

実行方法: python3 -m unittest discover -s examples/depletion -p 'test_openmc_depletion_resume.py'
"""

import tempfile
import unittest
from pathlib import Path

import openmc_depletion_resume as resume


def make_step(work: Path, k: int, *, complete: bool, total: int = 6) -> Path:
    step = work / f"step-{k}"
    step.mkdir(parents=True)
    (step / "depletion_results.h5").write_bytes(b"results")
    if complete:
        resume.write_complete(step, k, total)
    return step


class Steps(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work = Path(self.tmp.name) / "run"
        self.work.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def test_only_steps_with_the_complete_marker_count(self):
        make_step(self.work, 1, complete=True)
        make_step(self.work, 2, complete=True)
        make_step(self.work, 3, complete=False)
        self.assertEqual(resume.complete_steps(self.work), [1, 2])

    def test_steps_are_ordered_by_number_not_by_name(self):
        for k in (2, 10, 1):
            make_step(self.work, k, complete=True)
        self.assertEqual(resume.complete_steps(self.work), [1, 2, 10])

    def test_no_steps_before_the_first_run(self):
        self.assertEqual(resume.complete_steps(self.work), [])
        self.assertEqual(resume.complete_steps(self.work / "missing"), [])

    def test_the_complete_marker_records_the_step_and_the_total(self):
        step = make_step(self.work, 3, complete=True, total=6)
        self.assertEqual(resume.read_complete(step), (3, 6))

    def test_a_fresh_run_starts_at_step_1(self):
        self.assertEqual(resume.first_step(self.work, restart=False), 1)

    def test_a_fresh_run_refuses_a_directory_with_steps(self):
        make_step(self.work, 1, complete=True)
        with self.assertRaises(SystemExit):
            resume.first_step(self.work, restart=False)

    def test_a_restart_continues_after_the_last_complete_step(self):
        make_step(self.work, 1, complete=True)
        make_step(self.work, 2, complete=True)
        self.assertEqual(resume.first_step(self.work, restart=True), 3)

    def test_a_restart_drops_the_incomplete_step(self):
        make_step(self.work, 1, complete=True)
        incomplete = make_step(self.work, 2, complete=False)
        self.assertEqual(resume.first_step(self.work, restart=True), 2)
        self.assertFalse(incomplete.exists())
        self.assertTrue((self.work / "step-1" / "complete").exists())

    def test_a_restart_never_touches_a_complete_step(self):
        step = make_step(self.work, 1, complete=True)
        before = (step / "depletion_results.h5").read_bytes()
        make_step(self.work, 2, complete=False)
        resume.first_step(self.work, restart=True)
        self.assertEqual((step / "depletion_results.h5").read_bytes(), before)

    def test_a_restart_without_a_complete_step_is_an_error(self):
        make_step(self.work, 1, complete=False)
        with self.assertRaises(SystemExit):
            resume.first_step(self.work, restart=True)

    def test_a_restart_needs_only_the_last_complete_step(self):
        """VM は再開で最後の区間だけを S3 から戻すので、前の区間はなくてよい。"""
        make_step(self.work, 4, complete=True)
        self.assertEqual(resume.first_step(self.work, restart=True), 5)

    def test_a_step_directory_without_a_number_is_not_a_step(self):
        (self.work / "step-x").mkdir()
        (self.work / "step-1.bak").mkdir()
        self.assertEqual(resume.step_numbers(self.work), [])


class Backends(unittest.TestCase):
    """区間に分けた実行は、入力の問題で、ピンセル（1B〜1D）と集合体（2 で始まる）のツールを選ぶ（ADR-0178）。"""

    def write(self, problem: str) -> Path:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "input.toml"
        path.write_text(f'problem = "{problem}"\n')
        return path

    def test_a_pin_cell_problem_uses_the_pin_cell_tools(self):
        backend = resume.backend_of(self.write("1C"))
        self.assertEqual(
            backend.result_files, ("depletion_results.h5", "transport-seconds.json")
        )
        self.assertIsNone(backend.collect)

    def test_an_assembly_problem_also_carries_the_pin_powers(self):
        backend = resume.backend_of(self.write("2A"))
        self.assertEqual(
            backend.result_files,
            ("depletion_results.h5", "transport-seconds.json", "pins.json"),
        )
        self.assertIsNotNone(backend.collect)

    def test_a_missing_problem_is_an_error(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "input.toml"
        path.write_text("particles = 1\n")
        with self.assertRaises(SystemExit):
            resume.backend_of(path)


if __name__ == "__main__":
    unittest.main()
