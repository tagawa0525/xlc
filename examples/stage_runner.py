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
import shutil
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
    sysfs: Path = Path("/sys/devices/system/cpu"), allowed: set[int] | None = None
) -> int:
    """このプロセスに許された CPU の中の物理コアの数（ADR-0067）。

    xs の既定（物理コアの数を許された CPU の数で抑える。ADR-0189）は、Slurm の `-c 16` で 16 を返すが、許された 16 の
    ハードウェアスレッドは 8 コアにあたる。runner は各ステージに `--threads` を渡すので、許された CPU のコアを数える
    """
    if allowed is None:
        allowed = os.sched_getaffinity(0)
    cores = set()
    for cpu in allowed:
        topo = sysfs / f"cpu{cpu}" / "topology"
        cores.add(
            (
                (topo / "physical_package_id").read_text().strip(),
                (topo / "core_id").read_text().strip(),
            )
        )
    if not cores:
        raise RunnerError(f"no CPU topology under {sysfs}")
    return len(cores)


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


def launchable(
    queue: list[tuple[str, int]], free: int, budget: int
) -> list[tuple[str, int]]:
    """queue の順に、空きに収まるものを取る（収まらないものは飛ばして後ろを見る）。

    収まらないものは、走っているものが終わるのを待つ（空きの少ないスレッドで始めると、そのステージが最後まで遅い）。
    何も走っておらず（free == budget）1 つも収まらなければ、先頭を予算の全部で取る
    """
    taken = []
    for name, threads in queue:
        if threads <= free:
            taken.append((name, threads))
            free -= threads
    if not taken and queue and free == budget:
        taken.append((queue[0][0], budget))
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
    # followlinks: deps のディレクトリの中の symlink の先も DVC は読む
    for dirpath, dirnames, filenames in os.walk(
        path, onerror=raise_error, followlinks=True
    ):
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
            or not math.isfinite(seconds)
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


def unblock_sigterm() -> None:
    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGTERM})


@dataclass
class Finished:
    name: str
    cpu: float


def execute(
    queue,
    cmds,
    cwd,
    budget,
    finish,
    log=print,
    kill_grace: float = 10.0,
    outputs: dict | None = None,
) -> list[str]:
    """queue の順にステージを投入し、終わったもの（成功）を finish に渡す。失敗したステージの名前を返す。

    finish は commit できなかったステージの名前を返し、それも失敗に数える。失敗が出たら新しいステージは投入せず、
    走っているものは終わらせる。finish が例外を出したら、走っているものを止めて例外を出し直す。
    """
    queue = list(queue)
    free = budget
    running = {}  # pid → (名前, スレッドの数, Popen, 始めた時刻)
    failures: list[str] = []
    done: list[Finished] = []
    try:
        while True:
            if not failures:
                for name, threads in launchable(queue, free, budget):
                    cmd = with_threads(cmds[name], threads)
                    # dvc repro と同じく、実行の直前に出力を消す（投入しないステージの出力は残す）
                    remove_outputs((outputs or {}).get(name, ()), Path())
                    # 自分のセッションに置き、止めるときにシェルの子ごと止める。投入から running に入れるまでの間に
                    # SIGTERM の SystemExit が来ると、後始末がその子を知らずに残すので、その間は SIGTERM を保留する
                    signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
                    try:
                        # 子は保留の mask を受け継ぐので、exec の前に外す（runner は 1 スレッドなので preexec_fn でよい）
                        proc = subprocess.Popen(
                            cmd,
                            shell=True,
                            cwd=cwd,
                            start_new_session=True,
                            preexec_fn=unblock_sigterm,  # noqa: PLW1509
                        )
                        running[proc.pid] = (name, threads, proc, time.monotonic())
                    finally:
                        unblock_sigterm()
                    queue = [item for item in queue if item[0] != name]
                    free -= threads
                    log(f"start {name} threads={threads}")
            # 空いたスレッドに次を投入してから commit する（commit の数秒の間、コアを空けない）
            if done:
                if rejected := finish(done):
                    failures.extend(rejected)
                done = []
            if not running:
                if failures and queue:
                    log(
                        f"not launched after the failure: {', '.join(n for n, _ in queue)}"
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
                # Popen に回収済みと知らせる（wait4 で回収したので、Popen が後で回収を試みないように）
                proc.returncode = os.waitstatus_to_exitcode(status)
                free += threads
                elapsed = time.monotonic() - start
                # wait4 の rusage は、シェルが待った子（xs）を含む
                cpu = usage.ru_utime + usage.ru_stime
                if proc.returncode != 0:
                    failures.append(name)
                    log(f"FAILED {name} exit={proc.returncode} elapsed={elapsed:.1f}")
                else:
                    done.append(Finished(name, cpu))
                    log(
                        f"done {name} threads={threads} elapsed={elapsed:.1f} cpu={cpu:.1f}"
                    )
    except BaseException:
        # 後始末の間の 2 度目の SIGTERM で、残りのステージを止めずに抜けないように
        previous = signal.signal(signal.SIGTERM, signal.SIG_IGN)
        try:
            for _, _, proc, _ in running.values():
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
            deadline = time.monotonic() + kill_grace
            for _, _, proc, _ in running.values():
                try:
                    proc.wait(timeout=max(0.0, deadline - time.monotonic()))
                except subprocess.TimeoutExpired:
                    # SIGTERM を無視するか、止まれない子は、猶予の後にセッションごと SIGKILL する
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    proc.wait()
        finally:
            signal.signal(signal.SIGTERM, previous)
        raise


def dvc(*args: str, capture: bool = False) -> str:
    # subprocess.run は例外（SIGTERM からの SystemExit を含む）で子を SIGKILL し、dvc が .dvc/tmp/rwlock と
    # dvc.lock を書きかけで残す。SIGTERM で止めて、終わるのを待つ
    # 呼び手が SIGTERM を保留していても、dvc は受け取れるように mask を外す（runner は 1 スレッド）
    proc = subprocess.Popen(
        ["dvc", *args],
        text=True,
        stdout=subprocess.PIPE if capture else None,
        preexec_fn=unblock_sigterm,  # noqa: PLW1509
    )
    try:
        out, _ = proc.communicate()
    except BaseException:
        previous = signal.signal(signal.SIGTERM, signal.SIG_IGN)
        try:
            if proc.poll() is None:
                proc.terminate()
            proc.wait()
        finally:
            signal.signal(signal.SIGTERM, previous)
        raise
    if proc.returncode != 0:
        raise RunnerError(
            f"dvc {' '.join(args)} failed with exit code {proc.returncode}"
        )
    return out if capture else ""


def settle_upstream(status, repro, dvc_yaml, group: str, log) -> dict:
    """組の外（上流）が changed でなくなるまで `dvc repro` し、最後の `dvc status` を返す。

    dvc status は、上流が出力を書き換えた後に changed になる下流を示さない（build の後の xs-range）ので、1 回では足りない。
    同じ上流が続けて changed なら（always_changed など）止まる
    """
    previous = None
    while True:
        current = status()
        _, others = split_status(current, dvc_yaml, group)
        if not others:
            return current
        if others == previous:
            raise RunnerError(f"upstream still changed after dvc repro: {others}")
        log(f"upstream changed, running dvc repro {' '.join(others)}")
        repro(others)
        previous = others


def stages_to_run(group_status: dict) -> list[str]:
    """deps か cmd が変わったステージ。出力だけが lock と違うステージがあれば止まる（dvc repro なら run-cache から戻す）"""
    outputs_only = []
    for stage, entries in group_status.items():
        if not any(not isinstance(e, dict) or "changed deps" in e for e in entries):
            outputs_only.append(stage)
    if outputs_only:
        raise RunnerError(
            "outputs differ from dvc.lock with the deps unchanged; restore them "
            f"(scripts/worktree-restore-outputs.sh) or find out why: {outputs_only}"
        )
    return sorted(group_status)


def stage_outputs(
    lock_text: str, stage: str, directory: Path, status_entries: list
) -> set[str]:
    """ステージの出力（リポジトリの直下から）: dvc.lock が記録したものと、dvc status が changed outs に挙げたもの"""
    outs = set()
    if f"  {stage}:" in lock_text.split("\n"):
        outs = {
            os.path.normpath(directory / path)
            for path, _ in record_provenance.parse_stage(lock_text, stage)["outs"]
        }
    for entry in status_entries:
        if isinstance(entry, dict):
            outs.update(entry.get("changed outs", {}))
    return outs


def remove_outputs(paths, root: Path) -> None:
    """dvc repro と同じく、実行の前に出力を消す（xs が書かずに終わったとき、古い出力を新しい deps で commit しない）"""
    for path in paths:
        target = root / path
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink(missing_ok=True)


def stage_names(stage_list: str, dvc_yaml) -> set[str]:
    """`dvc stage list --name-only` の出力から、`<dvc.yaml>:` を外したステージの名前（記録のキーの形）"""
    prefix = f"{dvc_yaml}:"
    names = set()
    for line in stage_list.split():
        if not line.startswith(prefix):
            raise RunnerError(
                f"dvc stage list printed a stage outside {dvc_yaml}: {line}"
            )
        names.add(line.removeprefix(prefix))
    return names


def prune_record(record: dict, group: str, stages: set[str]) -> dict:
    """組の記録のうち、今の組に無いステージ（xs.toml から消した核種）を消す。記録の無いステージの見積もりに効くため"""
    return {
        k: v
        for k, v in record.items()
        if not ("@" in k and group_of(k) == group) or k in stages
    }


def commit_stages(
    names, commit, lock: Path, directory: Path, snapshot, root: Path
) -> list[str]:
    """names を commit し、commit できなかったステージの名前を返す。

    commit の前に deps の指紋を比べ、変わっていれば commit しない。commit の後に、dvc.lock が記録した deps が見たパスに
    入り、指紋が変わっていないことを確かめ、成り立たなければ dvc.lock を commit の前に戻して止まる。
    """
    if changed := changed_paths(snapshot, root):
        raise RunnerError(
            f"deps changed during the run: {changed}. Not committed: {names}"
        )
    before = lock.read_text() if lock.exists() else None
    failed = []
    try:
        commit(names)
        committed = list(names)
    except RunnerError:
        if len(names) == 1:
            return list(names)
        # どのステージが commit できないか（出力が無い、など）を 1 つずつ見る
        committed = []
        for name in names:
            try:
                commit([name])
                committed.append(name)
            except RunnerError:
                failed.append(name)
    text = lock.read_text()
    problem = None
    for name in committed:
        paths = lock_dep_paths(text, name, directory)
        if paths is None or not paths <= snapshot.keys():
            problem = (
                f"dvc.lock records deps of {name} the runner did not watch: "
                f"{sorted((paths or set()) - snapshot.keys())}"
            )
            break
    else:
        if changed := changed_paths(snapshot, root):
            problem = f"deps changed while committing {committed}: {changed}"
    if problem:
        if before is None:
            lock.unlink()
        else:
            tmp = lock.with_name(lock.name + ".tmp")
            tmp.write_text(before)
            tmp.replace(lock)
        raise RunnerError(
            f"{problem}. Restored dvc.lock to before this commit; rerun the stages"
        )
    return failed


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
    status = settle_upstream(
        lambda: json.loads(
            dvc("status", "--json", "--with-deps", target, capture=True)
        ),
        lambda keys: dvc("repro", *keys),
        dvc_yaml,
        group,
        log,
    )
    mine, _ = split_status(status, dvc_yaml, group)
    if not mine:
        log(f"{target}: nothing changed")
        return 0
    keys = stages_to_run(mine)
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
    # cmd は指紋の後に解決する（その間に dvc.yaml や xs.toml が変われば、commit の前の比較で見つかる）
    cmds = parse_dry_run(dvc("repro", "--dry", "-f", "-s", *keys, capture=True))
    if set(cmds) != set(keys):
        raise RunnerError(f"dvc repro --dry listed {sorted(cmds)}, expected {keys}")
    record_path = dvc_yaml.parent / RECORD_NAME
    names_in_yaml = stage_names(
        dvc("stage", "list", "--name-only", str(dvc_yaml), capture=True), dvc_yaml
    )
    record = prune_record(read_record(record_path), group, names_in_yaml)
    queue = plan([k.removeprefix(prefix) for k in keys], record, budget)
    log(
        f"budget={budget} stages={len(queue)}: "
        + ", ".join(f"{n}:{t}" for n, t in queue)
    )

    def finish(done: list[Finished]) -> list[str]:
        names = [f.name for f in done]
        # commit、確認、戻しの途中で SIGTERM の SystemExit が来ると、確かめていない commit や書きかけの dvc.lock が残る。
        # その間は SIGTERM を保留する（dvc の子は mask を外して受け取る）
        # 記録と結果の行を書き終えるまで保留する（commit したステージの時間と、commit できなかったステージを残す）
        signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGTERM})
        try:
            failed = commit_stages(
                names,
                lambda batch: dvc("commit", "-f", "-q", *(prefix + n for n in batch)),
                lock,
                dvc_yaml.parent,
                snapshot,
                Path(),
            )
            for name in failed:
                log(f"FAILED {name}: dvc commit failed")
            committed = [f for f in done if f.name not in failed]
            for f in committed:
                record[f.name] = round(f.cpu, 1)
            write_record(record_path, record)
            if committed:
                log(f"committed {', '.join(f.name for f in committed)}")
        finally:
            unblock_sigterm()
        return failed

    failures = execute(
        queue,
        {k.removeprefix(prefix): v for k, v in cmds.items()},
        dvc_yaml.parent,
        budget,
        finish,
        log,
        outputs={
            k.removeprefix(prefix): stage_outputs(
                lock_text, k.removeprefix(prefix), dvc_yaml.parent, mine[k]
            )
            for k in keys
        },
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
    start = time.monotonic()

    def log(message: str) -> None:
        print(f"[{time.monotonic() - start:7.1f}] {message}", flush=True)

    try:
        git = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
            capture_output=True,
            check=False,
        )
        if git.returncode != 0:
            raise RunnerError(f"not in a git repository: {git.stderr.strip()}")
        root = Path(git.stdout.strip()).resolve()
        dvc_yaml = args.dvc_yaml.resolve().relative_to(root)
        os.chdir(root)
        budget = physical_cores() if args.threads is None else args.threads
        stop_on_sigterm()
        code = run_group(dvc_yaml, args.group, budget, log)
    except (
        RunnerError,
        record_provenance.ProvenanceError,
        OSError,
        ValueError,
    ) as e:
        # 誤りの行で終わる（ValueError は relative_to と JSON の読み取り、OSError は指紋とファイルの読み書き）
        log(f"ERROR: {type(e).__name__}: {e}")
        return 1
    log(f"runner elapsed={time.monotonic() - start:.1f}")
    return code


if __name__ == "__main__":
    sys.exit(main())
