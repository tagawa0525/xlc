"""Problem 3 と 5-2D の参照解（vera の openmc_model_3d.py、openmc_core_2d.py）の、途中の statepoint と再開。

長い実行を、中断されうる計算機（スポットの VM）で流すための部品。openmc_run.py と分けるのは、ステージ（DVC）の deps に
入る openmc_run.py を変えず、problem 1、2 のステージを古くしないため。

- 入力の `statepoint-interval`（世代数）を書くと、その間隔と最後の世代で statepoint を書く。書かなければ最後の世代だけ
  （OpenMC の動作）。古い statepoint は、最新の `KEEP` 個を残して消す
- `run --statepoint-interval N` は、入力の `statepoint-interval` を上書きする（中断が続くときに、起動し直しだけ間隔を縮めるため。
  入力のファイルは変えない。ADR-0168）
- `run --restart` は、作業ディレクトリの最新の statepoint から `openmc -r` で再開する。入力は同じものを与える（OpenMC は
  乱数の系列と tally の積算を statepoint から引き継ぐ）。再開した世代から先の statepoint が同じ名前で書かれる
- `run --threads N` は、入力の `threads` を上書きする（計算機の物理コア数に合わせるため。値はコマンドの引数に残る）
- extract は、run が作業ディレクトリに書いた model.xml と入力から作ったモデルを `check_model` で比べ、違えば止まる（run は
  ステージでなく、入力やモデルのコードを変えても statepoint は古いまま残るため）
"""

import re
import sys
import tempfile
import threading
import xml.etree.ElementTree as ET
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import openmc

KEEP = 2
POLL_SECONDS = 30
OPTIONAL = {"statepoint-interval"}


@dataclass(frozen=True)
class Options:
    restart: bool = False
    threads: int | None = None
    statepoint_interval: int | None = None


def statepoints(work: Path) -> list[Path]:
    """作業ディレクトリの statepoint（`statepoint.<世代>.h5`）を、世代の順に。"""
    found = [
        (int(m.group(1)), p)
        for p in work.glob("statepoint.*.h5")
        if (m := re.fullmatch(r"statepoint\.(\d+)\.h5", p.name))
    ]
    return [p for _, p in sorted(found)]


def latest_readable(found: list[Path]) -> Path:
    """最新の読める statepoint。書いている途中で中断されて壊れたものは消し、その前のものを使う。"""
    while found:
        try:
            with openmc.StatePoint(found[-1]):
                return found[-1].resolve()
        except (
            OSError,
            KeyError,
            ValueError,
        ) as error:  # h5py が壊れたファイルで出す例外
            print(f"removing unreadable {found[-1].name}: {error}", flush=True)
            found.pop().unlink()
    raise SystemExit("no readable statepoint to restart from")


def set_statepoints(
    settings: openmc.Settings, run: dict[str, Any], options: Options
) -> None:
    """`statepoint-interval`（`options` の指定が入力より先）ごと（と最後）の世代に statepoint を書く。"""
    interval = options.statepoint_interval
    if interval is None:
        interval = run.get("statepoint-interval")
    if interval is None:
        return
    if interval <= 0:
        raise SystemExit("statepoint-interval must be positive")
    last = settings.batches
    settings.statepoint = {"batches": [*range(interval, last, interval), last]}


def prune(work: Path) -> None:
    """最新の `KEEP` 個を除いて statepoint を消す（最新のものは書いている途中でありうるので、古い方だけ）。"""
    for old in statepoints(work)[:-KEEP]:
        old.unlink(missing_ok=True)


def run_model(
    model: openmc.Model, run: dict[str, Any], work: Path, options: Options
) -> None:
    """`model` を `work` で解く。`options.restart` なら最新の statepoint から再開する。"""
    # `model.run(cwd=work)` はプロセスの作業ディレクトリを `work` に移すので、相対パスのままでは別のスレッドの
    # `prune` が `work/work` を探す
    work = work.resolve()
    restart_file = None
    found = statepoints(work)
    if options.restart:
        # model.run は model.xml を書き換えるので、前の区切りのモデルと比べるのはここが最後（extract は最後の区切りの model.xml しか見ない）
        check_model(model, work)
        restart_file = latest_readable(found)
        print(f"restart from {restart_file}", flush=True)
    elif found:
        raise SystemExit(
            f"{work} already has statepoints ({found[-1].name}); add --restart to continue, or remove them to start over"
        )
    stop = threading.Event()

    def prune_while_running() -> None:
        while not stop.wait(POLL_SECONDS):
            prune(work)

    pruner = threading.Thread(target=prune_while_running, daemon=True)
    pruner.start()
    try:
        model.run(
            cwd=work,
            threads=options.threads or run["threads"],
            output=False,
            restart_file=restart_file,
        )
    finally:
        stop.set()
        pruner.join()
    prune(work)


# run の環境と引数で変わり、モデルを変えない要素: ライブラリのパス（実行した worktree の絶対パス）と、statepoint を書く世代
# （`--statepoint-interval` で上書きできる。ADR-0168）
RUN_ONLY = (("materials", "cross_sections"), ("settings", "state_point"))


def model_form(path: Path) -> list[tuple[str, str]]:
    """model.xml の最上位の要素ごとの正規形（`RUN_ONLY` を除き、空白を無視した C14N）。"""
    root = ET.parse(path).getroot()
    for parent, tag in RUN_ONLY:
        for node in root.findall(parent):
            for child in node.findall(tag):
                node.remove(child)
    return [
        (node.tag, ET.canonicalize(ET.tostring(node), strip_text=True)) for node in root
    ]


def check_model(model: openmc.Model, work: Path) -> None:
    """`work` の model.xml（run が書いた、statepoint を出したモデル）が `model` と `RUN_ONLY` のほかで同じことを確かめる。"""
    ran = work / "model.xml"
    if not ran.is_file():
        raise SystemExit(
            f"{ran} is missing; cannot tell which model the statepoints came from"
        )
    with tempfile.TemporaryDirectory() as tmp:
        expected = Path(tmp) / "model.xml"
        model.export_to_model_xml(str(expected))
        want = model_form(expected)
    got = model_form(ran)
    if got != want:
        differ = sorted({tag for tag, _ in set(got) ^ set(want)})
        where = (
            f"in {differ}"
            if differ
            else "in the order or the count of its top-level elements"
        )
        raise SystemExit(
            f"{ran} differs from the model of the input {where}; run it again or fix the input"
        )


def run_or_extract(
    solve: Callable[[Path, Path, Options], None],
    extract: Callable[[Path, Path, Path], None],
    usage: str,
) -> None:
    """コマンドの引数で、run（`--restart`、`--threads N`、`--statepoint-interval N` つき）か extract を呼ぶ（ADR-0123、ADR-0168）。"""
    args = sys.argv[1:]
    if args[:1] == ["run"]:
        rest = args[1:]
        restart = False
        threads = None
        interval = None
        while rest and rest[0].startswith("--"):
            match rest:
                case ["--restart", *tail]:
                    restart = True
                case ["--threads", n, *tail] if (
                    n.isascii() and n.isdigit() and int(n) > 0
                ):
                    threads = int(n)
                case ["--statepoint-interval", n, *tail] if (
                    n.isascii() and n.isdigit() and int(n) > 0
                ):
                    interval = int(n)
                case _:
                    raise SystemExit(usage)
            rest = tail
        if len(rest) != 2:
            raise SystemExit(usage)
        solve(Path(rest[0]), Path(rest[1]), Options(restart, threads, interval))
    elif len(args) == 4 and args[0] == "extract":
        extract(Path(args[1]), Path(args[2]), Path(args[3]))
    else:
        raise SystemExit(usage)
