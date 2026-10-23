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

    def test_stages_recorded_as_zero_share_the_budget_equally(self):
        self.assertEqual(
            sr.plan(["g@a", "g@b"], {"g@a": 0.0, "g@b": 0.0}, 4),
            [("g@a", 2), ("g@b", 2)],
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
        self.assertEqual(sr.launchable(queue, 5, budget=8), [("a", 4), ("c", 1)])

    def test_waits_while_others_run_when_nothing_fits(self):
        # 走っているものが終われば収まるので待つ（少ないスレッドで始めると、そのステージが最後まで遅い）
        self.assertEqual(sr.launchable([("a", 4), ("b", 5)], 3, budget=8), [])

    def test_stage_larger_than_the_budget_takes_the_whole_budget_when_idle(self):
        # 収まる b を先に入れ、a は何も走っていなくなったとき予算の全部で入る
        self.assertEqual(sr.launchable([("a", 10), ("b", 1)], 8, budget=8), [("b", 1)])
        self.assertEqual(sr.launchable([("a", 10)], 8, budget=8), [("a", 8)])


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
    "build": [{"changed deps": {"crates": "modified"}}],
}


class StatusTest(unittest.TestCase):
    def test_separates_the_group_from_upstream(self):
        group, others = sr.split_status(STATUS, YAML, "g")
        self.assertEqual(sorted(group), [f"{YAML}:g@a", f"{YAML}:g@b"])
        self.assertEqual(sorted(others), ["build", f"{YAML}:range"])

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

            # 同じ大きさで書き直す。mtime の刻みが粗いファイルシステムでも違うように、mtime を 1 秒進める
            mtime = (d / "file").stat().st_mtime_ns
            (d / "file").write_text("2")
            os.utime(d / "file", ns=(mtime + 10**9, mtime + 10**9))
            self.assertEqual(sr.changed_paths(before, d), ["file"])

            (d / "dir" / "g").write_text("1")
            self.assertEqual(sr.changed_paths(before, d), ["dir", "file"])

    def test_unreadable_directory_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "dir" / "sub").mkdir(parents=True)
            (d / "dir" / "sub").chmod(0)
            try:
                with self.assertRaises(PermissionError):
                    sr.fingerprints(["dir"], d)
            finally:
                (d / "dir" / "sub").chmod(0o755)

    def test_sees_changes_under_a_symlinked_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "target").mkdir()
            (d / "target" / "f").write_text("1")
            (d / "dir").mkdir()
            (d / "dir" / "link").symlink_to(d / "target")
            before = sr.fingerprints(["dir"], d)
            (d / "target" / "g").write_text("1")
            self.assertEqual(sr.changed_paths(before, d), ["dir"])


class PhysicalCoresTest(unittest.TestCase):
    def test_counts_cores_not_hardware_threads_within_the_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # 2 パッケージ × 2 コア × 2 スレッド。core_id はダイごとに 0 から振られうるので、すべて 0 にしておく
            for cpu in range(8):
                topo = root / f"cpu{cpu}" / "topology"
                topo.mkdir(parents=True)
                (topo / "physical_package_id").write_text("0\n")
                (topo / "core_id").write_text("0\n")
                first = cpu - cpu % 2
                (topo / "thread_siblings_list").write_text(f"{first}-{first + 1}\n")
            (root / "cpufreq").mkdir()
            self.assertEqual(sr.physical_cores(root, allowed=set(range(8))), 4)
            # CPU 0 と 1 は同じコア、2 は別のコア
            self.assertEqual(sr.physical_cores(root, allowed={0, 1, 2}), 2)
            self.assertEqual(sr.physical_cores(root, allowed={0, 2, 4, 6}), 4)


class RecordTest(unittest.TestCase):
    def test_rejects_times_that_are_not_non_negative_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stage-cpu-seconds.json"
            for bad in ["-1.5", "true", '"3"', "null", "NaN", "Infinity"]:
                path.write_text(f'{{"g@a": 2.0, "g@b": {bad}}}')
                with self.assertRaises(sr.RunnerError, msg=bad):
                    sr.read_record(path)

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

    def test_removes_outputs_just_before_launching_each_stage(self):
        # a は失敗し、b は投入されない。b の出力は消さない
        (self.dir / "a-old").write_text("old")
        (self.dir / "b-old").write_text("old")
        failures = sr.execute(
            [("a", 1), ("b", 1)],
            {"a": "./xs a 3", "b": "./xs b 0"},
            self.dir,
            1,
            lambda done: None,
            outputs={"a": [str(self.dir / "a-old")], "b": [str(self.dir / "b-old")]},
        )
        self.assertEqual(failures, ["a"])
        self.assertFalse((self.dir / "a-old").exists())
        self.assertTrue((self.dir / "b-old").exists())

    def test_stage_ignoring_sigterm_is_killed_after_the_grace(self):
        (self.dir / "stubborn").write_text("#!/bin/sh\ntrap '' TERM\nsleep 30\n")
        (self.dir / "stubborn").chmod(0o755)
        (self.dir / "xs").write_text(
            '#!/bin/sh\n[ "$3" = b ] && exec ./stubborn\nexit 0\n'
        )

        def finish(done):
            raise sr.RunnerError("deps changed")

        start = time.monotonic()
        with self.assertRaises(sr.RunnerError):
            sr.execute(
                [("a", 1), ("b", 1)],
                {"a": "./xs a", "b": "./xs b"},
                self.dir,
                2,
                finish,
                kill_grace=0.5,
            )
        self.assertLess(time.monotonic() - start, 10)

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

    def test_stage_asking_more_than_the_budget_runs_with_the_budget(self):
        failures = sr.execute(
            [("a", 3)], {"a": "./xs a 0"}, self.dir, 2, lambda done: None
        )
        self.assertEqual(failures, [])
        self.assertEqual((self.dir / "a.out").read_text(), "2\n")

    def test_stages_the_finish_handler_rejects_are_failures(self):
        # finish が返したステージ（commit できなかったもの）は失敗として数え、その後は投入しない。b は a の commit の前に
        # 空いたスレッドへ投入済み
        cmds = {"a": "./xs a 0", "b": "./xs b 0", "c": "./xs c 0"}
        failures = sr.execute(
            [("a", 1), ("b", 1), ("c", 1)],
            cmds,
            self.dir,
            1,
            lambda done: [f.name for f in done],
        )
        self.assertEqual(failures, ["a", "b"])
        self.assertFalse((self.dir / "c.out").exists())


# 偽の dvc: 自分の pid を書き、SIGTERM を受けたら got-term を書いて終わる。引数 fail なら 3 で終わる
FAKE_DVC = """#!/bin/sh
[ "$1" = fail ] && exit 3
trap 'echo term > got-term; exit 0' TERM
echo $$ > pid
sleep 30 &
wait
"""


class DvcTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        dvc = self.dir / "dvc"
        dvc.write_text(FAKE_DVC)
        dvc.chmod(0o755)
        self.env = {**os.environ, "PATH": f"{self.dir}:{os.environ['PATH']}"}

    def tearDown(self):
        self.tmp.cleanup()

    def test_failure_is_a_runner_error(self):
        old = os.environ["PATH"]
        os.environ["PATH"] = self.env["PATH"]
        try:
            with self.assertRaises(sr.RunnerError):
                sr.dvc("fail")
        finally:
            os.environ["PATH"] = old

    def test_sigterm_reaches_dvc_as_sigterm(self):
        # runner が SIGTERM で止まるとき、待っていた dvc を SIGKILL せず SIGTERM で止め、終わるのを待つこと
        script = "import stage_runner as sr\nsr.stop_on_sigterm()\nsr.dvc('wait')\n"
        env = {**self.env, "PYTHONPATH": str(Path(__file__).resolve().parent)}
        # 偽の dvc は cwd に pid と got-term を書く
        runner = subprocess.Popen([sys.executable, "-c", script], cwd=self.dir, env=env)
        pid_file = self.dir / "pid"
        deadline = time.monotonic() + 10
        while not (pid_file.exists() and pid_file.read_text().strip()):
            self.assertLess(time.monotonic(), deadline)
            time.sleep(0.05)
        runner.send_signal(signal.SIGTERM)
        self.assertNotEqual(runner.wait(timeout=10), 0)
        self.assertTrue((self.dir / "got-term").exists())


class UpstreamTest(unittest.TestCase):
    def test_runs_upstream_until_only_the_group_is_changed(self):
        # build を実行すると xs-range が changed になる（dvc status は上流が出力を書き換えた後の下流を示さない）
        statuses = [
            {"build": [], f"{YAML}:g@a": []},
            {f"{YAML}:range": [], f"{YAML}:g@a": []},
            {f"{YAML}:g@a": [{"changed deps": {"bin/xs": "modified"}}]},
        ]
        repros = []
        status = sr.settle_upstream(
            lambda: statuses.pop(0), repros.append, YAML, "g", lambda m: None
        )
        self.assertEqual(repros, [["build"], [f"{YAML}:range"]])
        self.assertEqual(list(status), [f"{YAML}:g@a"])

    def test_upstream_in_another_order_is_still_stuck(self):
        statuses = [{"a": [], "b": []}, {"b": [], "a": []}]
        with self.assertRaises(sr.RunnerError):
            sr.settle_upstream(
                lambda: statuses.pop(0), lambda keys: None, YAML, "g", lambda m: None
            )

    def test_upstream_that_comes_back_is_an_error(self):
        # always_changed の上流が出力を変え、下流が 1 回おきに出る（[a] → [a, b] → [a]）
        statuses = [{"a": []}, {"a": [], "b": []}, {"a": []}, {"a": [], "b": []}]
        repros = []
        with self.assertRaises(sr.RunnerError):
            sr.settle_upstream(
                lambda: statuses.pop(0), repros.append, YAML, "g", lambda m: None
            )
        self.assertEqual(len(repros), 2)

    def test_upstream_that_stays_changed_is_an_error(self):
        status = {"always": ["always changed"], f"{YAML}:g@a": []}
        with self.assertRaises(sr.RunnerError):
            sr.settle_upstream(
                lambda: status, lambda keys: None, YAML, "g", lambda m: None
            )


class StagesToRunTest(unittest.TestCase):
    def test_runs_stages_with_changed_deps_or_command(self):
        group = {
            f"{YAML}:g@a": [{"changed deps": {"bin/xs": "modified"}}],
            f"{YAML}:g@b": ["changed command"],
        }
        self.assertEqual(sr.stages_to_run(group), [f"{YAML}:g@a", f"{YAML}:g@b"])

    def test_stage_with_only_changed_outputs_is_an_error(self):
        # deps は lock のままで出力だけが違う（戻していない、手で変えた）。runner は作り直さず止まる
        group = {
            f"{YAML}:g@a": [{"changed deps": {"bin/xs": "modified"}}],
            f"{YAML}:g@b": [{"changed outs": {"examples/p/b.msgpack": "deleted"}}],
        }
        with self.assertRaises(sr.RunnerError):
            sr.stages_to_run(group)


class CommitTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "p").mkdir()
        (self.root / "p" / "in.txt").write_text("1")
        self.lock = self.root / "p" / "dvc.lock"
        self.lock.write_text("schema: '2.0'\nstages:\n")
        self.snapshot = sr.fingerprints(["p/in.txt"], self.root)
        self.run_cache = self.root / "runs"
        (self.run_cache / "aa").mkdir(parents=True)
        (self.run_cache / "aa" / "old").write_text("old")

    def tearDown(self):
        self.tmp.cleanup()

    def fake_commit(self, dep="in.txt", fail=(), touch=False, cmd="./xs {}"):
        # dvc commit の代わり: lock に cmd と deps を書き、run-cache に 1 件足す
        def commit(names):
            if any(n in fail for n in names):
                raise sr.RunnerError("commit failed")
            with self.lock.open("a") as f:
                for n in names:
                    out_md5 = f"{n:f>32}"
                    f.write(f"  {n}:\n    cmd: {cmd.format(n)}\n")
                    f.write(f"    deps:\n    - path: {dep}\n      hash: md5\n")
                    f.write(f"      md5: {'0' * 32}\n")
                    f.write(f"    outs:\n    - path: {n}.out\n      hash: md5\n")
                    f.write(f"      md5: {out_md5}\n")
                    (self.run_cache / "aa" / f"new-{n}").write_text(f"md5: {out_md5}\n")
            # 同じ間にほかの worktree が足した run-cache の項目
            (self.run_cache / "aa" / "other").write_text(f"md5: {'e' * 32}\n")
            if touch:
                mtime = (self.root / "p" / "in.txt").stat().st_mtime_ns + 10**9
                os.utime(self.root / "p" / "in.txt", ns=(mtime, mtime))

        return commit

    def commit(self, names, commit):
        return sr.commit_stages(
            names,
            commit,
            self.lock,
            Path("p"),
            self.snapshot,
            self.root,
            cmds={n: f"./xs {n}" for n in names},
            run_cache=self.run_cache,
        )

    def test_stage_that_cannot_be_committed_is_a_failure_and_the_rest_are_committed(
        self,
    ):
        failed = self.commit(["a", "x", "b"], self.fake_commit(fail={"x"}))
        self.assertEqual(failed, ["x"])
        text = self.lock.read_text()
        self.assertIn("  a:", text)
        self.assertIn("  b:", text)
        self.assertNotIn("  x:", text)

    def test_deps_changed_while_committing_restores_the_lock_and_run_cache(self):
        before = self.lock.read_text()
        with self.assertRaises(sr.RunnerError):
            self.commit(["a"], self.fake_commit(touch=True))
        self.assertEqual(self.lock.read_text(), before)
        # dvc commit が足した run-cache の項目（commit した出力の md5 を持つもの）も消す。前からあるものと、ほかの
        # worktree が同じ間に足したものは残す
        self.assertEqual(
            sorted(p.name for p in self.run_cache.rglob("*") if p.is_file()),
            ["old", "other"],
        )

    def test_command_changed_since_the_run_restores_the_lock(self):
        # 実行した cmd と、commit が dvc.lock に書いた cmd が違う（実行の間に dvc.yaml の vars が変わった）
        before = self.lock.read_text()
        with self.assertRaises(sr.RunnerError):
            self.commit(["a"], self.fake_commit(cmd="./xs {} --other"))
        self.assertEqual(self.lock.read_text(), before)

    def test_unwatched_deps_in_the_lock_restore_the_lock(self):
        before = self.lock.read_text()
        with self.assertRaises(sr.RunnerError):
            self.commit(["a"], self.fake_commit(dep="other.txt"))
        self.assertEqual(self.lock.read_text(), before)

    def test_commit_cut_by_a_held_sigterm_is_not_retried_stage_by_stage(self):
        # SIGTERM を保留している間に batch の commit が落ちたら、1 つずつ commit し直さずに止まる（止める途中で dvc を
        # 増やさない）
        calls = []

        def commit(names):
            calls.append(names)
            raise sr.RunnerError("commit failed")

        previous = signal.signal(signal.SIGTERM, lambda *args: None)
        signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
        try:
            os.kill(os.getpid(), signal.SIGTERM)
            with self.assertRaises(sr.RunnerError):
                self.commit(["a", "b"], commit)
        finally:
            signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})
            signal.signal(signal.SIGTERM, previous)
        self.assertEqual(calls, [["a", "b"]])

    def test_single_stage_is_not_committed_twice(self):
        calls = []

        def commit(names):
            calls.append(names)
            raise sr.RunnerError("commit failed")

        self.assertEqual(self.commit(["x"], commit), ["x"])
        self.assertEqual(calls, [["x"]])

    def test_deps_changed_before_committing_commits_nothing(self):
        mtime = (self.root / "p" / "in.txt").stat().st_mtime_ns + 10**9
        os.utime(self.root / "p" / "in.txt", ns=(mtime, mtime))
        calls = []
        with self.assertRaises(sr.RunnerError):
            self.commit(["a"], calls.append)
        self.assertEqual(calls, [])


class LockCommandTest(unittest.TestCase):
    def test_joins_a_folded_command(self):
        self.assertEqual(
            sr.lock_cmd(LOCK, "g@a"),
            "../../bin/xs nuclide xs.toml xs-range.toml a xs-nuclides/a.msgpack",
        )
        self.assertEqual(sr.lock_cmd(LOCK, "g@b"), "echo b")


class OutputsTest(unittest.TestCase):
    def test_lists_locked_and_reported_outputs_from_the_root(self):
        status = [
            {"changed deps": {"bin/xs": "modified"}},
            {"changed outs": {"examples/p/extra.msgpack": "deleted"}},
        ]
        self.assertEqual(
            sr.stage_outputs(LOCK, "g@a", Path("examples/p"), status),
            {"examples/p/xs-nuclides/a.msgpack", "examples/p/extra.msgpack"},
        )

    def test_removes_files_and_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "f").write_text("1")
            (root / "d").mkdir()
            (root / "d" / "x").write_text("1")
            sr.remove_outputs(["f", "d", "missing"], root)
            self.assertEqual(list(root.iterdir()), [])


class PruneTest(unittest.TestCase):
    def test_reads_stage_names_without_the_dvc_yaml_prefix(self):
        # dvc stage list --name-only の出力の形（DVC 3.67.1）
        text = f"{YAML}:range\n{YAML}:g@a\n{YAML}:g@b\n"
        self.assertEqual(sr.stage_names(text, YAML), {"range", "g@a", "g@b"})

    def test_drops_records_of_stages_no_longer_in_the_group(self):
        record = {"g@a": 1.0, "g@gone": 2278.0, "h@x": 3.0}
        self.assertEqual(
            sr.prune_record(record, "g", {"g@a", "g@b"}), {"g@a": 1.0, "h@x": 3.0}
        )


class ArgumentsTest(unittest.TestCase):
    def test_dvc_yaml_outside_the_repository_is_an_error_not_a_traceback(self):
        self.assertEqual(sr.main(["/nonexistent/dvc.yaml", "g", "--threads", "1"]), 1)

    def test_threads_must_be_positive(self):
        with self.assertRaises(SystemExit):
            sr.main(["examples/p/dvc.yaml", "g", "--threads", "0"])


if __name__ == "__main__":
    unittest.main()
