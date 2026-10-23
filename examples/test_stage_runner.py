"""stage_runner.py の検査（スレッドの割り当てと順序、投入、DVC の出力の読み方、記録、deps の指紋、実行）。

実行: python3 -m unittest discover -s examples -p 'test_stage_runner.py'
"""

import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stage_runner as sr

YAML = "examples/p/dvc.yaml"


class PlanTest(unittest.TestCase):
    def test_orders_longest_first_and_gives_threads_in_proportion(self):
        record = {"g@a": 10.0, "g@b": 50.0, "g@c": 20.0, "g@d": 20.0}
        # C = 100、B = 4: 単位 25 秒。b は ceil(50 / 25) = 2、c と d は ceil(0.8) = 1、a は ceil(0.4) = 1
        self.assertEqual(
            sr.plan(["g@a", "g@b", "g@c", "g@d"], record, 4),
            [("g@b", 2), ("g@c", 1), ("g@d", 1), ("g@a", 1)],
        )

    def test_clamps_threads_to_the_budget(self):
        self.assertEqual(sr.plan(["g@a"], {"g@a": 5.0}, 16), [("g@a", 16)])

    def test_unrecorded_stage_takes_the_longest_recorded_time_of_the_group(self):
        record = {"g@a": 30.0, "g@b": 10.0, "h@x": 1000.0}
        # 新しい n は g の最大の 30 秒とみなす（h の記録は別の組）。C = 30 + 10 + 30 = 70、B = 7: 単位 10 秒
        self.assertEqual(
            sr.plan(["g@n", "g@b", "g@a"], record, 7),
            [("g@a", 3), ("g@n", 3), ("g@b", 1)],
        )

    def test_without_records_the_budget_is_shared_equally(self):
        self.assertEqual(
            sr.plan(["g@b", "g@a", "g@c"], {}, 6),
            [("g@a", 2), ("g@b", 2), ("g@c", 2)],
        )


class LaunchableTest(unittest.TestCase):
    def test_takes_stages_in_order_that_fit_and_skips_the_rest(self):
        queue = [("a", 4), ("b", 3), ("c", 1), ("d", 2)]
        # 空き 5: a（4）を入れて残り 1、b（3）と d（2）は飛ばし、c（1）を入れる
        self.assertEqual(sr.launchable(queue, 5), [("a", 4), ("c", 1)])

    def test_nothing_fits(self):
        self.assertEqual(sr.launchable([("a", 4)], 3), [])


class DryRunTest(unittest.TestCase):
    def test_maps_each_stage_to_its_command(self):
        text = (
            "Stage 'build' didn't change, skipping\n"
            f"Running stage '{YAML}:g@a':\n"
            "> ../../bin/xs nuclide xs.toml xs-range.toml a xs-nuclides/a.msgpack\n"
            "\n"
            f"Running stage '{YAML}:g@b':\n"
            "> ../../bin/xs nuclide xs.toml xs-range.toml b xs-nuclides/b.msgpack\n"
            "Use `dvc push` to send your updates to remote storage.\n"
        )
        self.assertEqual(
            sr.parse_dry_run(text),
            {
                f"{YAML}:g@a": "../../bin/xs nuclide xs.toml xs-range.toml a xs-nuclides/a.msgpack",
                f"{YAML}:g@b": "../../bin/xs nuclide xs.toml xs-range.toml b xs-nuclides/b.msgpack",
            },
        )

    def test_rejects_a_stage_with_several_commands(self):
        text = f"Running stage '{YAML}:g@a':\n> echo 1\n> echo 2\n"
        with self.assertRaises(sr.RunnerError):
            sr.parse_dry_run(text)


class WithThreadsTest(unittest.TestCase):
    def test_puts_the_option_right_after_the_executable(self):
        self.assertEqual(
            sr.with_threads("../../bin/xs nuclide xs.toml  'a b'", 3),
            "../../bin/xs --threads 3 nuclide xs.toml  'a b'",
        )

    def test_rejects_other_executables(self):
        with self.assertRaises(sr.RunnerError):
            sr.with_threads("python3 run.py", 2)


STATUS = {
    f"{YAML}:g@a": [{"changed deps": {"bin/xs": "modified"}}],
    f"{YAML}:g@b": [
        {"changed deps": {"examples/p/xs.toml": {"nuclide.b": "modified"}}},
        {"changed outs": {"examples/p/xs-nuclides/b.msgpack": "deleted"}},
    ],
    f"{YAML}:range": [{"changed deps": {"bin/xs": "modified"}}],
    "dvc.yaml:build": [{"changed deps": {"crates": "modified"}}],
}


class StatusTest(unittest.TestCase):
    def test_separates_the_group_from_upstream(self):
        group, others = sr.split_status(STATUS, YAML, "g")
        self.assertEqual(sorted(group), [f"{YAML}:g@a", f"{YAML}:g@b"])
        self.assertEqual(sorted(others), ["dvc.yaml:build", f"{YAML}:range"])

    def test_collects_changed_dependency_paths(self):
        self.assertEqual(
            sr.changed_dep_paths(STATUS),
            {"bin/xs", "examples/p/xs.toml", "crates"},
        )


LOCK = """schema: '2.0'
stages:
  g@a:
    cmd: ../../bin/xs nuclide xs.toml xs-range.toml a
      xs-nuclides/a.msgpack
    deps:
    - path: ../../bin/xs
      hash: md5
      md5: 6b502338fb57b4a096f3e960e2fb8970
      size: 2846736
      isexec: true
    - path: ../../references/endf/thermal
      hash: md5
      md5: dbf65c78ca90c9c9266163dc97d2cbe4.dir
      size: 225043459
      nfiles: 91
    params:
      xs.toml:
        nuclide.a:
          temperatures:
          - 565
        sigma0:
        - 10.0
    outs:
    - path: xs-nuclides/a.msgpack
      hash: md5
      md5: de5d942f0b6dd214a26dcdc08b18340f
      size: 3761182
  g@b:
    cmd: echo b
"""


class LockTest(unittest.TestCase):
    def test_reads_dependency_and_params_paths_from_the_root(self):
        self.assertEqual(
            sr.lock_dep_paths(LOCK, "g@a", Path("examples/p")),
            {"bin/xs", "references/endf/thermal", "examples/p/xs.toml"},
        )

    def test_stage_without_dependencies(self):
        self.assertEqual(sr.lock_dep_paths(LOCK, "g@b", Path("examples/p")), set())

    def test_stage_not_in_the_lock(self):
        self.assertIsNone(sr.lock_dep_paths(LOCK, "g@c", Path("examples/p")))


class FingerprintTest(unittest.TestCase):
    def test_changes_when_a_file_is_rewritten_or_added(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "dir").mkdir()
            (d / "dir" / "f").write_text("1")
            (d / "file").write_text("1")
            paths = ["dir", "file", "missing"]
            before = sr.fingerprints(paths, d)
            self.assertEqual(sr.fingerprints(paths, d), before)
            self.assertIsNone(before["missing"])

            time.sleep(0.01)
            (d / "file").write_text("2")
            self.assertEqual(sr.changed_paths(before, d), ["file"])

            (d / "dir" / "g").write_text("1")
            self.assertEqual(sr.changed_paths(before, d), ["dir", "file"])


class PhysicalCoresTest(unittest.TestCase):
    def test_counts_cores_not_hardware_threads_within_the_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # 2 パッケージ × 2 コア × 2 スレッド
            for cpu in range(8):
                topo = root / f"cpu{cpu}" / "topology"
                topo.mkdir(parents=True)
                (topo / "physical_package_id").write_text(f"{cpu // 4}\n")
                (topo / "core_id").write_text(f"{(cpu // 2) % 2}\n")
            (root / "cpufreq").mkdir()
            self.assertEqual(sr.physical_cores(root, allowed=8), 4)
            self.assertEqual(sr.physical_cores(root, allowed=3), 3)


class RecordTest(unittest.TestCase):
    def test_missing_record_is_empty_and_writes_sorted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stage-cpu-seconds.json"
            self.assertEqual(sr.read_record(path), {})
            sr.write_record(path, {"g@b": 2.0, "g@a": 1.25})
            self.assertEqual(sr.read_record(path), {"g@a": 1.25, "g@b": 2.0})
            self.assertEqual(path.read_text(), '{\n  "g@a": 1.25,\n  "g@b": 2.0\n}\n')


# 引数（--threads <数> <名前> <終了コード> [<秒>]）を受けて、<秒> 待ってから <名前>.out にスレッドの数を書き、終了コードで
# 終わる偽の xs
FAKE_XS = """#!/bin/sh
[ "$1" = --threads ] || exit 99
sleep "${5:-0}"
echo "$2" > "$3.out"
exit "$4"
"""


class ExecuteTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        xs = self.dir / "xs"
        xs.write_text(FAKE_XS)
        xs.chmod(0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def test_runs_every_stage_with_its_threads_and_reports_each_finish(self):
        finished = []
        cmds = {name: f"./xs {name} 0" for name in "abc"}
        failures = sr.execute(
            [("a", 2), ("b", 1), ("c", 1)],
            cmds,
            self.dir,
            3,
            lambda done: finished.extend(f.name for f in done),
        )
        self.assertEqual(failures, [])
        self.assertEqual(sorted(finished), ["a", "b", "c"])
        for name, threads in [("a", 2), ("b", 1), ("c", 1)]:
            self.assertEqual((self.dir / f"{name}.out").read_text(), f"{threads}\n")

    def test_waits_for_free_threads_before_launching(self):
        # 予算 2: a（2）が終わるまで b を投入しない。b の投入は a の終わりの後
        finished = []
        cmds = {"a": "./xs a 0 0.3", "b": "./xs b 0"}
        sr.execute(
            [("a", 2), ("b", 1)],
            cmds,
            self.dir,
            2,
            lambda done: finished.extend(f.name for f in done),
        )
        self.assertEqual(finished, ["a", "b"])

    def test_failure_stops_new_launches_but_lets_running_stages_finish(self):
        finished = []
        # 予算 2: a（すぐ失敗）と b（0.5 秒）を同時に実行する。c は a の失敗の後に空きができるが、投入しない
        cmds = {"a": "./xs a 3", "b": "./xs b 0 0.5", "c": "./xs c 0"}
        failures = sr.execute(
            [("a", 1), ("b", 1), ("c", 1)],
            cmds,
            self.dir,
            2,
            lambda done: finished.extend(f.name for f in done),
        )
        self.assertEqual(failures, ["a"])
        self.assertEqual(finished, ["b"])
        self.assertFalse((self.dir / "c.out").exists())

    def test_error_in_the_finish_handler_stops_running_stages(self):
        cmds = {"a": "./xs a 0", "b": "./xs b 0 30"}

        def finish(done):
            raise sr.RunnerError("deps changed")

        start = time.monotonic()
        with self.assertRaises(sr.RunnerError):
            sr.execute([("a", 1), ("b", 1)], cmds, self.dir, 2, finish)
        self.assertLess(time.monotonic() - start, 10)
        self.assertFalse((self.dir / "b.out").exists())

    def test_sigterm_stops_running_stages(self):
        # 偽の xs は自分の pid を書いてから 30 秒待つ。runner を SIGTERM で止めたとき、自分のセッションに置いた xs も止まること
        (self.dir / "xs").write_text("#!/bin/sh\necho $$ > pid\nsleep 30\n")
        script = (
            "import sys; from pathlib import Path; import stage_runner as sr\n"
            "sr.stop_on_sigterm()\n"
            "sr.execute([('a', 1)], {'a': './xs a'}, Path(sys.argv[1]), 1, lambda done: None)\n"
        )
        runner = subprocess.Popen(
            [sys.executable, "-c", script, str(self.dir)],
            cwd=Path(__file__).resolve().parent,
        )
        pid_file = self.dir / "pid"
        deadline = time.monotonic() + 10
        while not (pid_file.exists() and pid_file.read_text().strip()):
            self.assertLess(time.monotonic(), deadline)
            time.sleep(0.05)
        xs_pid = int(pid_file.read_text())
        runner.send_signal(signal.SIGTERM)
        self.assertNotEqual(runner.wait(timeout=10), 0)
        with self.assertRaises(ProcessLookupError):
            os.kill(xs_pid, 0)


if __name__ == "__main__":
    unittest.main()
