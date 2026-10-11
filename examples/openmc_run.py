"""VERA の参照解のモデルのスクリプト（examples/vera の openmc_model.py、openmc_model_3d.py、openmc_core_2d.py）で共通に
使う、計算の流れの部品: 入力の読み込み、run と extract のコマンド、初期の中性子源、メッシュの tally の取り出し、k と
ピン出力の書き出し。

ライブラリの段階（openmc_library.py）が deps にする openmc_common.py と分けるのは、モデルの側の変更で、核種ごとの
ライブラリの段階を DVC が古いとみなさないようにするため（ADR-0130）。
"""

import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import msgpack
import numpy as np
import openmc
from openmc_common import input_path, read_input


def box_settings(
    particles: int,
    batches: int,
    inactive: int,
    seed: int,
    lower: tuple[float, float, float],
    upper: tuple[float, float, float],
) -> openmc.Settings:
    """世代の数と乱数の種、箱（角 `lower`、`upper`）の中の燃料に一様な初期の中性子源。"""
    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box(lower, upper),
        constraints={"fissionable": True},
    )
    return settings


# 参照解の入力（p1/openmc-1B.toml の形）のキー。計算の条件はすべて入力に書く（ADR-0130）
RUN_KEYS = {
    "problem",
    "cross-sections",
    "particles",
    "batches",
    "inactive",
    "seed",
    "threads",
}


def read_run(
    source: Path, problems: list[str], optional: set[str] | None = None
) -> dict[str, Any]:
    """参照解の入力を読み、問題が `problems` にあることを確かめ、`cross-sections` を絶対パスにする。"""
    run = read_input(source, RUN_KEYS, optional or set())
    if run["problem"] not in problems:
        raise SystemExit(
            f"unknown problem {run['problem']}; expected one of {problems}"
        )
    # 絶対パスにする。OpenMC も相対パスを Python の作業ディレクトリから絶対パスに直す（openmc.config の resolve_paths）
    # ので、outs の model.xml には worktree の絶対パスが入る。statepoint が日時と計算時間を含み outs はバイト単位では
    # 再現しないので、それを避ける書き方はしない
    run["cross-sections"] = input_path(source, run["cross-sections"])
    return run


def write_k(path: Path, problem: str, k) -> None:
    """k（statepoint の `keff`）を JSON で書き、表示する。"""
    path.write_text(
        json.dumps({"problem": problem, "k": k.nominal_value, "k-std": k.std_dev})
        + "\n"
    )
    print(f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}")


def fuel_mean(rows: np.ndarray) -> float:
    """核分裂のあるセル（値が正のセル）の平均。ピン出力は、これで割って燃料のピンの平均で 1 にする。"""
    return rows[rows > 0].mean()


def write_pin_powers(path: Path, rows: np.ndarray) -> None:
    """セルごとの値（[行（上から）][列]）を、lattice の `pin-powers` と同じ形式（`xlc-pin-powers`）で書く。"""
    path.write_bytes(msgpack.packb({"format": "xlc-pin-powers", "rows": rows.tolist()}))


def run_or_extract(
    solve: Callable[[Path, Path], None],
    extract: Callable[[Path, Path, Path], None],
    usage: str,
) -> None:
    """コマンドの引数で、run（OpenMC を解く）か extract（statepoint から結果を抜く）を呼ぶ（ADR-0123）。"""
    match sys.argv[1:]:
        case ["run", source, work]:
            solve(Path(source), Path(work))
        case ["extract", source, work, results]:
            extract(Path(source), Path(work), Path(results))
        case _:
            raise SystemExit(usage)


def mesh_values(
    sp: openmc.StatePoint, name: str, shape: tuple[int, ...]
) -> tuple[np.ndarray, np.ndarray]:
    """メッシュの tally `name` の核分裂率の平均と統計誤差を、形 `shape`（最後の 2 つが y と x）にし、y を上から並べる
    （メッシュは x が最も速く、y は下から）。"""
    tally = sp.get_tally(name=name)
    return tuple(
        np.flip(tally.get_values(scores=["fission"], value=v).reshape(shape), axis=-2)
        for v in ("mean", "std_dev")
    )
