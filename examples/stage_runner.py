"""foreach のステージの組のうち changed のものを、前回の CPU 時間の長い順に、スレッドの数を決めて同時に実行し、終わったものを
`dvc commit` で dvc.lock に取り込む（ADR-0193）。

使い方（default の devShell で。dvc を references の devShell で実行しない。ADR-0130）:
    python3 examples/stage_runner.py examples/<問題>/dvc.yaml xs-nuclide [--threads <予算>]

- `dvc status --with-deps` で組の changed を拾う。組の外の上流（build、xs-range）が changed なら、先に `dvc repro` で実行する
- cmd は `dvc repro --dry -f -s` が出す解決済みのもの。実行ファイル（xs）の直後に `--threads <数>` を足して、dvc.yaml の
  ディレクトリで実行する。dvc.yaml の cmd と dvc.lock の cmd には `--threads` が入らない（出力はスレッドの数によらない。ADR-0189）
- 前回の CPU 時間は dvc.yaml の隣の stage-cpu-seconds.json（git で追跡する）。commit のたびに書き直す
- 始めに deps の指紋を取り、commit の前後で比べる。変わっていれば止まる
- 標準ライブラリだけ（default の devShell の Python に YAML の読み手はない。dvc.lock は record_provenance.parse_stage で読む）
"""

import argparse
import json
import math
import os
import re
import signal
import stat
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import record_provenance

RECORD_NAME = "stage-cpu-seconds.json"


class RunnerError(Exception):
    pass


def physical_cores(
    sysfs: Path = Path("/sys/devices/system/cpu"), allowed: int | None = None
) -> int:
    """物理コアの数を、このプロセスに許された CPU の数で抑えたもの（xs の既定と同じ。ADR-0067、ADR-0189）"""
    cores = {
        (
            (topo / "physical_package_id").read_text().strip(),
            (topo / "core_id").read_text().strip(),
        )
        for topo in sysfs.glob("cpu[0-9]*/topology")
    }
    if not cores:
        raise RunnerError(f"no CPU topology under {sysfs}")
    if allowed is None:
        allowed = len(os.sched_getaffinity(0))
    return min(len(cores), allowed)


def group_of(stage: str) -> str:
    return stage.split("@", 1)[0]


def plan(
    stages: list[str], record: dict[str, float], budget: int
) -> list[tuple[str, int]]:
    """(ステージ, スレッドの数) を CPU 時間の長い順に。スレッドの数は ceil(c / (C / budget)) を 1〜budget に抑えたもの"""
    times = {}
    for stage in stages:
        if stage in record:
            times[stage] = record[stage]
            continue
        group = [
            v for k, v in record.items() if "@" in k and group_of(k) == group_of(stage)
        ]
        # 記録の無いステージは組の最長とみなす。組に記録が無ければすべて同じ値（予算を等分する）
        times[stage] = max(group) if group else 1.0
    total = sum(times.values())
    if total == 0:
        # どれも記録の刻み（0.1 秒）より短い。同じ長さなので予算を等分する
        times = dict.fromkeys(stages, 1.0)
        total = float(len(stages))
    order = sorted(stages, key=lambda s: (-times[s], s))
    # round: 割り切れる比が浮動小数点の誤差で ceil を 1 つ上げないように
    return [
        (s, min(budget, max(1, math.ceil(round(times[s] * budget / total, 9)))))
        for s in order
    ]


def launchable(queue: list[tuple[str, int]], free: int) -> list[tuple[str, int]]:
    """queue の順に、空きに収まるものを取る（収まらないものは飛ばして後ろを見る）"""
    taken = []
    for name, threads in queue:
        if threads <= free:
            taken.append((name, threads))
            free -= threads
    return taken


def parse_dry_run(text: str) -> dict[str, str]:
    """`dvc repro --dry` の出力から、ステージ → cmd"""
    cmds: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if m := re.fullmatch(r"Running stage '(.+)':", line):
            current = m.group(1)
            cmds[current] = []
        elif line.startswith("> "):
            if current is None:
                raise RunnerError(
                    f"command without a stage in dvc repro --dry: {line!r}"
                )
            cmds[current].append(line[2:])
    for stage, lines in cmds.items():
        if len(lines) != 1:
            raise RunnerError(f"{stage} has {len(lines)} commands; the runner runs one")
    return {stage: lines[0] for stage, lines in cmds.items()}


def with_threads(cmd: str, threads: int) -> str:
    m = re.match(r"(\S+)(.*)\Z", cmd, re.DOTALL)
    if not m or Path(m.group(1)).name != "xs":
        raise RunnerError(
            f"the command does not start with xs, which takes --threads: {cmd!r}"
        )
    return f"{m.group(1)} --threads {threads}{m.group(2)}"


def split_status(
    status: dict, dvc_yaml: str | Path, group: str
) -> tuple[dict, list[str]]:
    """`dvc status --json` を、組のステージと組の外（上流）に分ける"""
    prefix = f"{dvc_yaml}:{group}@"
    mine = {k: v for k, v in status.items() if k.startswith(prefix)}
    return mine, [k for k in status if not k.startswith(prefix)]


def changed_dep_paths(status: dict) -> set[str]:
    """`dvc status --json` が changed と示した deps と params のファイル（リポジトリの直下から）"""
    paths = set()
    for entries in status.values():
        for entry in entries:
            if isinstance(entry, dict):
                paths.update(entry.get("changed deps", {}))
    return paths


def lock_dep_paths(lock_text: str, stage: str, directory: Path) -> set[str] | None:
    """dvc.lock の stage の deps と params のファイル（リポジトリの直下から）。stage が無ければ None"""
    lines = lock_text.split("\n")
    if f"  {stage}:" not in lines:
        return None
    deps = record_provenance.parse_stage(lock_text, stage)["deps"]
    paths = {path for path, _ in deps}
    start = lines.index(f"  {stage}:") + 1
    in_params = False
    for line in lines[start:]:
        if re.match(r"  \S", line):
            break
        if line == "    params:":
            in_params = True
        elif re.match(r"    \S", line):
            in_params = False
        elif in_params and (m := re.fullmatch(r"      ([^\s-][^:]*):", line)):
            paths.add(m.group(1))
    return {os.path.normpath(directory / p) for p in paths}


def raise_error(error: OSError):
    raise error


def fingerprint(path: Path):
    """内容が変わったかを見る指紋（大きさ、mtime、inode。ディレクトリは中のファイルすべて）。無ければ None"""
    try:
        st = path.stat()
    except FileNotFoundError:
        return None
    if not stat.S_ISDIR(st.st_mode):
        return (st.st_size, st.st_mtime_ns, st.st_ino)
    items = []
    # onerror: os.walk は読めないディレクトリを黙って飛ばし、前後の指紋が同じに見える
    for dirpath, dirnames, filenames in os.walk(path, onerror=raise_error):
        dirnames.sort()
        for name in sorted(filenames):
            file = Path(dirpath) / name
            s = file.stat()
            items.append(
                (str(file.relative_to(path)), s.st_size, s.st_mtime_ns, s.st_ino)
            )
    return tuple(items)


def fingerprints(paths, root: Path) -> dict:
    return {p: fingerprint(root / p) for p in paths}


def changed_paths(before: dict, root: Path) -> list[str]:
    return sorted(p for p, fp in before.items() if fingerprint(root / p) != fp)


def read_record(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    record = json.loads(path.read_text())
    if not isinstance(record, dict):
        raise RunnerError(f"{path} is not a map from stage to CPU seconds")
    for stage, seconds in record.items():
        # bool は int の子クラスなので先に外す
        if (
            isinstance(seconds, bool)
            or not isinstance(seconds, (int, float))
            or seconds < 0
        ):
            raise RunnerError(
                f"{path}: CPU seconds of {stage} is not a non-negative number: {seconds!r}"
            )
    return record


def write_record(path: Path, record: dict[str, float]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(dict(sorted(record.items())), indent=2) + "\n")
    tmp.replace(path)


def stop_on_sigterm() -> None:
    """SIGTERM を SystemExit にする。既定の動作は後始末をせずに終わり、自分のセッションに置いたステージが残る"""

    def handler(signum, frame):
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, handler)


@dataclass
class Finished:
    name: str
    threads: int
    elapsed: float
    cpu: float


def execute(queue, cmds, cwd, budget, finish, log=print) -> list[str]:
    """queue の順にステージを投入し、終わったもの（成功）を finish に渡す。失敗したステージの名前を返す。

    失敗が出たら新しいステージは投入せず、走っているものは終わらせる。finish が例外を出したら、走っているものを止めて
    例外を出し直す。
    """
    queue = list(queue)
    free = budget
    running = {}  # pid → (名前, スレッドの数, Popen, 始めた時刻)
    failures: list[str] = []
    done: list[Finished] = []
    try:
        while True:
            if not failures:
                for name, threads in launchable(queue, free):
                    cmd = with_threads(cmds[name], threads)
                    # 自分のセッションに置き、止めるときにシェルの子ごと止める
                    proc = subprocess.Popen(
                        cmd, shell=True, cwd=cwd, start_new_session=True
                    )
                    running[proc.pid] = (name, threads, proc, time.monotonic())
                    queue.remove((name, threads))
                    free -= threads
                    log(f"start {name} threads={threads}")
            # 空いたスレッドに次を投入してから commit する（commit の数秒の間、コアを空けない）
            if done:
                finish(done)
                done = []
            if not running:
                if failures and queue:
                    log(
                        f"not launched after the failure: {', '.join(n for n, _ in queue)}"
                    )
                elif queue:
                    raise RunnerError(
                        f"stages need more threads than the budget {budget}: {queue}"
                    )
                return failures
            reaped = []
            while not reaped:
                # 待つのは自分の子だけ（wait4(-1) はほかの子も回収しうる）
                for pid in list(running):
                    got, status, usage = os.wait4(pid, os.WNOHANG)
                    if got:
                        reaped.append((running.pop(pid), status, usage))
                if not reaped:
                    time.sleep(0.1)
            for (name, threads, proc, start), status, usage in reaped:
                proc.returncode = os.waitstatus_to_exitcode(status)
                free += threads
                elapsed = time.monotonic() - start
                # wait4 の rusage は、シェルが待った子（xs）を含む
                cpu = usage.ru_utime + usage.ru_stime
                if proc.returncode != 0:
                    failures.append(name)
                    log(f"FAILED {name} exit={proc.returncode} elapsed={elapsed:.1f}")
                else:
                    done.append(Finished(name, threads, elapsed, cpu))
                    log(
                        f"done {name} threads={threads} elapsed={elapsed:.1f} cpu={cpu:.1f}"
                    )
    except BaseException:
        for _, _, proc, _ in running.values():
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for _, _, proc, _ in running.values():
            proc.wait()
        raise


def dvc(*args: str, capture: bool = False) -> str:
    # subprocess.run は例外（SIGTERM からの SystemExit を含む）で子を SIGKILL し、dvc が .dvc/tmp/rwlock と
    # dvc.lock を書きかけで残す。SIGTERM で止めて、終わるのを待つ
    proc = subprocess.Popen(
        ["dvc", *args], text=True, stdout=subprocess.PIPE if capture else None
    )
    try:
        out, _ = proc.communicate()
    except BaseException:
        if proc.poll() is None:
            proc.terminate()
        proc.wait()
        raise
    if proc.returncode != 0:
        raise RunnerError(
            f"dvc {' '.join(args)} failed with exit code {proc.returncode}"
        )
    return out if capture else ""


def run_group(dvc_yaml: Path, group: str, budget: int, log) -> int:
    if dvc_yaml == Path("dvc.yaml"):
        # DVC はリポジトリの直下の dvc.yaml のステージを `dvc.yaml:` を付けずに示すので、組を見分けられない
        raise RunnerError(
            "the runner does not take the dvc.yaml at the root of the repository"
        )
    if re.search(r"^\s*wdir:", dvc_yaml.read_text(), re.MULTILINE):
        raise RunnerError(
            f"{dvc_yaml} has wdir; the runner runs commands in the directory of dvc.yaml"
        )
    target = f"{dvc_yaml}:{group}"
    status = json.loads(dvc("status", "--json", "--with-deps", target, capture=True))
    mine, others = split_status(status, dvc_yaml, group)
    if others:
        log(f"upstream changed, running dvc repro {' '.join(others)}")
        dvc("repro", *others)
        status = json.loads(
            dvc("status", "--json", "--with-deps", target, capture=True)
        )
        mine, others = split_status(status, dvc_yaml, group)
        if others:
            raise RunnerError(f"upstream still changed after dvc repro: {others}")
    if not mine:
        log(f"{target}: nothing changed")
        return 0
    keys = sorted(mine)
    cmds = parse_dry_run(dvc("repro", "--dry", "-f", "-s", *keys, capture=True))
    if set(cmds) != set(keys):
        raise RunnerError(f"dvc repro --dry listed {sorted(cmds)}, expected {keys}")
    prefix = f"{dvc_yaml}:"
    lock = dvc_yaml.parent / "dvc.lock"
    lock_text = lock.read_text() if lock.exists() else ""
    watched = {str(dvc_yaml)} | changed_dep_paths(status)
    for key in keys:
        watched |= (
            lock_dep_paths(lock_text, key.removeprefix(prefix), dvc_yaml.parent)
            or set()
        )
    snapshot = fingerprints(watched, Path())
    record_path = dvc_yaml.parent / RECORD_NAME
    record = read_record(record_path)
    queue = plan([k.removeprefix(prefix) for k in keys], record, budget)
    log(
        f"budget={budget} stages={len(queue)}: "
        + ", ".join(f"{n}:{t}" for n, t in queue)
    )

    def finish(done: list[Finished]) -> None:
        names = [f.name for f in done]
        if changed := changed_paths(snapshot, Path()):
            raise RunnerError(
                f"deps changed during the run: {changed}. Not committed: {names}"
            )
        dvc("commit", "-f", "-q", *(prefix + n for n in names))
        text = lock.read_text()
        for name in names:
            paths = lock_dep_paths(text, name, dvc_yaml.parent)
            if paths is None or not paths <= snapshot.keys():
                raise RunnerError(
                    f"dvc.lock records deps of {name} the runner did not watch: "
                    f"{sorted((paths or set()) - snapshot.keys())}. Do not commit dvc.lock; rerun the stage"
                )
        if changed := changed_paths(snapshot, Path()):
            raise RunnerError(
                f"deps changed while committing {names}: {changed}. Do not commit dvc.lock; rerun the stages"
            )
        for f in done:
            record[f.name] = round(f.cpu, 1)
        write_record(record_path, record)
        log(f"committed {', '.join(names)}")

    failures = execute(
        queue,
        {k.removeprefix(prefix): v for k, v in cmds.items()},
        dvc_yaml.parent,
        budget,
        finish,
        log,
    )
    if failures:
        log(f"ERROR: failed stages (not committed): {', '.join(failures)}")
        return 1
    return 0


def positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be positive: {value}")
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dvc_yaml", type=Path)
    parser.add_argument("group", help="foreach のステージの名前（例: xs-nuclide）")
    parser.add_argument(
        "--threads",
        type=positive_int,
        help="同時に使うスレッドの合計（既定: 物理コアの数）",
    )
    args = parser.parse_args(argv)
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
    )
    dvc_yaml = args.dvc_yaml.resolve().relative_to(root)
    os.chdir(root)
    budget = physical_cores() if args.threads is None else args.threads
    stop_on_sigterm()
    start = time.monotonic()

    def log(message: str) -> None:
        print(f"[{time.monotonic() - start:7.1f}] {message}", flush=True)

    try:
        code = run_group(dvc_yaml, args.group, budget, log)
    except RunnerError as e:
        log(f"ERROR: {e}")
        return 1
    log(f"runner elapsed={time.monotonic() - start:.1f}")
    return code


if __name__ == "__main__":
    sys.exit(main())
