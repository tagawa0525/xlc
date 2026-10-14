"""openmc_depletion.py の燃焼を、燃焼度の点（ステップ）ごとの区間に分けて流し、中断のあとにステップの境界から再開する（ADR-0174）。

使い方（references の devShell で）:

    python3 openmc_depletion_resume.py run [--restart] [--threads N] <入力> <作業ディレクトリ>
    python3 openmc_depletion_resume.py extract <入力> <作業ディレクトリ> <出力ディレクトリ>

<入力> は openmc_depletion.py と同じ TOML。`--threads` は入力の `threads` を上書きする。

OpenMC の燃焼は輸送の途中から再開できず、`depletion_results.h5` はステップごとに追記するので、書いている最中の中断で壊れうる。
そこで 1 ステップを 1 区間として別のプロセスで流す。区間 k は `step-<k>/` で、前の区間の `step-<k-1>/depletion_results.h5`
を写して `prev_results` に渡し（組成と時刻はそこから始まる）、残りの 1 ステップだけを進める。最後の区間以外は
`integrate(final_step=False)` で、終わりの輸送をしない。次の区間は先頭の輸送を、分けない走行と同じく自分でやり直す
（`BosTransport`）: OpenMC の再開は保存した反応率から始めるが、核分裂収率（`chain.fission_yields`）は直前の輸送のスペクトルで
決まって保存されないため、最初の半ステップが分けない走行と違う収率で解かれる（粒子 1000 の試行で、2 点目の k が 0.008 ずれた。
ADR-0174）。輸送の乱数の種と組成が同じなので先頭の輸送の結果は同じで、輸送の回数も分けない走行と変わらない（6 点で 13 回）。
区間の最後に完了の印 `complete`（内容は `<k>/<全ステップ数>`）を置き、`complete` のない区間は再開で捨てる。
元の区間は触らないので、失うのは多くて 1 ステップ。

各区間の `transport-seconds.json` は前の区間までの時間を含む。extract は最後の区間（`complete` の内容が `<k>/<k>`）の
`depletion_results.h5` と `transport-seconds.json` を openmc_depletion.py の extract に渡して
`<出力ディレクトリ>/<入力の名前>.msgpack` を書く。再開した VM は S3 から最後の区間だけを戻すので、前の区間がなくても書ける。
AWS の機種はまちまちなので、この時間は xlc の時間との比較に使わない（ADR-0174）。

`solve` の設定の組み立て（TimedOperator と CECMIntegrator の引数）は openmc_depletion.py の `solve` と同じ。そちらを変えると
`openmc-1C` の重い段階が古くなるので、変えずにここへ写している（#0132）。
"""

import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

STEP = re.compile(r"step-(\d+)$")
COMPLETE = re.compile(r"(\d+)/(\d+)\n$")


def step_numbers(work: Path) -> list[int]:
    """`work` の下の `step-<k>/` の k を昇順に（`complete` の有無によらず）。"""
    if not work.is_dir():
        return []
    return sorted(
        int(m.group(1))
        for p in work.iterdir()
        if p.is_dir() and (m := STEP.fullmatch(p.name))
    )


def complete_steps(work: Path) -> list[int]:
    """`complete` の印がある区間の番号を昇順に。"""
    return [
        k for k in step_numbers(work) if (work / f"step-{k}" / "complete").is_file()
    ]


def write_complete(step: Path, k: int, total: int) -> None:
    """区間の最後に、完了の印を置く（これより前に、区間のファイルはすべて書き終えておく）。"""
    (step / "complete").write_text(f"{k}/{total}\n")


def read_complete(step: Path) -> tuple[int, int]:
    """完了の印の `(k, 全ステップ数)`。"""
    text = (step / "complete").read_text()
    m = COMPLETE.fullmatch(text)
    if not m:
        raise SystemExit(f"{step / 'complete'}: unexpected content {text!r}")
    return int(m.group(1)), int(m.group(2))


def first_step(work: Path, *, restart: bool) -> int:
    """次に流す区間の番号。再開では、`complete` のない区間を捨てる。

    VM は再開で S3 から最後の区間だけを戻すので、それより前の区間は作業ディレクトリになくてよい。
    """
    steps = step_numbers(work)
    if not restart:
        if steps:
            raise SystemExit(
                f"{work} already has steps {steps}; pass --restart to continue them"
            )
        return 1
    done = complete_steps(work)
    if not done:
        raise SystemExit(f"{work} has no complete step to restart from")
    for k in steps:
        if k not in done:
            shutil.rmtree(work / f"step-{k}")
    return done[-1] + 1


def bos_transport_integrator():
    """先頭の輸送を再開でもやり直す `CECMIntegrator`（`openmc` の import を遅らせるため関数にしてある）。"""
    import openmc.deplete

    class BosTransport(openmc.deplete.CECMIntegrator):
        def _get_bos_data_from_restart(self, source_rate, bos_conc):
            conc = list(self.operator.prev_res[-1].data[0])
            return self._get_bos_data_from_operator(0, source_rate, conc)

    return BosTransport


def solve_step(source: Path, work: Path, k: int, threads: int | None) -> None:
    """区間 k を流す（別のプロセスで呼ばれる。`openmc.lib` を区間ごとに初期化し直さないため）。"""
    import openmc.deplete
    from openmc_depletion import (
        TRANSPORT_SECONDS,
        YIELD_CUTOFF,
        TimedOperator,
        build,
        read,
    )
    from openmc_model import FUEL_TEMPERATURE

    run = read(source)
    if threads is not None:
        run["threads"] = threads
    total = len(run["burnup"])
    step = work / f"step-{k}"
    if step.exists():
        raise SystemExit(f"{step} already exists")
    step.mkdir(parents=True)
    prev = None
    if k > 1:
        shutil.copy2(work / f"step-{k - 1}" / "depletion_results.h5", step)
        prev = openmc.deplete.Results(step / "depletion_results.h5")
    os.chdir(step)
    os.environ["OMP_NUM_THREADS"] = str(run["threads"])
    openmc.config["cross_sections"] = run["cross-sections"]
    model, _ = build(run)
    TimedOperator.fuel_temperature = FUEL_TEMPERATURE[run["problem"]]
    operator = TimedOperator(
        model,
        chain_file=run["chain"],
        prev_results=prev,
        normalization_mode="fission-q",
        fission_yield_mode="cutoff",
        fission_yield_opts={"cutoff": YIELD_CUTOFF},
        reaction_rate_mode="direct",
    )
    burnup = run["burnup"]
    steps = [b - a for a, b in itertools.pairwise([0.0, *burnup])]
    integrator = bos_transport_integrator()(
        operator,
        [steps[k - 1]],
        power_density=run["power"],
        timestep_units="MWd/kg",
        solver="cram48",
    )
    start = time.perf_counter()
    integrator.integrate(final_step=k == total)
    seconds = time.perf_counter() - start
    if k < total:
        # final_step=False でも終わりの `operator(...)` は呼ばれ（輸送はせず反応率 0 を返す）、時間の記録に 1 回分が入る
        TRANSPORT_SECONDS.pop()
    # 前の区間までの時間を引き継ぎ、extract が最後の区間だけから書けるようにする（再開した VM に前の区間はない）
    before = {"transport": [], "total": 0.0}
    if k > 1:
        before = json.loads(
            (work / f"step-{k - 1}" / "transport-seconds.json").read_text()
        )
    (step / "transport-seconds.json").write_text(
        json.dumps(
            {
                "transport": before["transport"] + TRANSPORT_SECONDS,
                "total": before["total"] + seconds,
                "transportable": TimedOperator.transportable,
            }
        )
        + "\n"
    )
    # 再開にも extract にも使わない（結果は depletion_results.h5）。区間ごとに 20 MB を S3 に置かない
    for statepoint in step.glob("statepoint.*.h5"):
        statepoint.unlink()
    write_complete(step, k, total)


def solve(source: Path, work: Path, *, restart: bool, threads: int | None) -> None:
    from openmc_depletion import read

    run = read(source)
    total = len(run["burnup"])
    threads = run["threads"] if threads is None else threads
    work.mkdir(parents=True, exist_ok=True)
    # OMP_NUM_THREADS は OpenMC を import する前の環境に要る（import したあとに設定しても、読み込み済みのライブラリには効かない）
    env = {**os.environ, "OMP_NUM_THREADS": str(threads)}
    for k in range(first_step(work, restart=restart), total + 1):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "interval",
            "--threads",
            str(threads),
        ]
        subprocess.run([*command, str(source), str(work), str(k)], check=True, env=env)
        if k not in complete_steps(work):
            raise SystemExit(f"step {k} ended without the complete marker")


def extract(source: Path, work: Path, output: Path) -> None:
    from openmc_depletion import extract as extract_one

    done = complete_steps(work)
    if not done:
        raise SystemExit(f"{work} has no complete step")
    last = work / f"step-{done[-1]}"
    k, total = read_complete(last)
    if k != total:
        raise SystemExit(f"{last} is step {k} of {total}; the run is not finished")
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        merged = Path(tmp)
        for name in ("depletion_results.h5", "transport-seconds.json"):
            (merged / name).symlink_to(last / name)
        extract_one(source, merged, output / f"{source.stem}.msgpack")


def main(argv: list[str]) -> None:
    restart = False
    threads = None
    args = []
    it = iter(argv)
    for a in it:
        if a == "--restart":
            restart = True
        elif a == "--threads":
            threads = int(next(it))
        else:
            args.append(a)
    match args:
        case ["run", source, work]:
            solve(
                Path(source).resolve(),
                Path(work).resolve(),
                restart=restart,
                threads=threads,
            )
        case ["interval", source, work, k]:
            solve_step(Path(source).resolve(), Path(work).resolve(), int(k), threads)
        case ["extract", source, work, output]:
            extract(
                Path(source).resolve(), Path(work).resolve(), Path(output).resolve()
            )
        case _:
            raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
