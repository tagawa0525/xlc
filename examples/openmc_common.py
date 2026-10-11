"""OpenMC の参照解のスクリプト（examples/vera の openmc_model.py、openmc_model_3d.py、openmc_core_2d.py と examples/mox の
openmc_model.py）で共通に使う部品。

入力の TOML の読み込み、run と extract のコマンド、k とピン出力の書き出し、物質、固有値計算の設定、燃料の核種ごとの実効断面積（燃料の領域の反応率 / (個数密度 × 中性子束)）の
tally と、その statepoint からの取り出し。
"""

import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import msgpack
import numpy as np
import openmc
import tomllib

# tally に書く score の名前のうち、statepoint で別の名前になるもの（弾性散乱は "(n,elastic)" になる）と、
# 書き出す名前と score の名前が違うもの（捕獲は "(n,gamma)"）
STORED = {"elastic": "(n,elastic)", "capture": "(n,gamma)"}


def check_keys(
    table: dict[str, Any], where: str, required: set[str], optional: set[str]
) -> None:
    """`table` に `required` のキーがすべてあり、`required` と `optional` のほかのキーがないことを確かめる。

    計算の条件は入力に書き、既定値を持たない（ADR-0122、ADR-0130）ので、足りないキーも綴りを誤ったキーも止める。
    """
    missing = required - table.keys()
    unknown = table.keys() - required - optional
    if missing or unknown:
        raise SystemExit(
            f"{where}: missing keys {sorted(missing)}, unknown keys {sorted(unknown)}"
        )


def read_input(path: Path, required: set[str], optional: set[str]) -> dict[str, Any]:
    """入力の TOML を読み、最上位のキーを `check_keys` で確かめる。"""
    table = tomllib.loads(path.read_text())
    check_keys(table, str(path), required, optional)
    return table


def input_path(input_file: Path, value: str) -> Path:
    """入力に書いたパス（入力のファイルからの相対）を、絶対パスにする。"""
    return (input_file.parent / value).resolve()


def material(
    name: str, densities: dict[str, float], temperature: float
) -> openmc.Material:
    """個数密度 (/(barn cm)) の辞書から物質を作る。"""
    m = openmc.Material(name=name, temperature=temperature)
    for nuclide, density in densities.items():
        m.add_nuclide(nuclide, density, "ao")
    m.set_density("sum")
    return m


def eigenvalue_settings(
    particles: int, batches: int, inactive: int, seed: int, half: float
) -> openmc.Settings:
    """世代の数と乱数の種、燃料に一様な初期の中性子源（`half` は外周の箱の半分の幅）。"""
    return box_settings(
        particles, batches, inactive, seed, (-half, -half, -1), (half, half, 1)
    )


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


def group_filter(library: Path) -> openmc.EnergyFilter:
    """xlc のライブラリ（msgpack）の群構造のエネルギーのフィルタ（OpenMC の境界は昇順）。"""
    bounds = msgpack.unpackb(library.read_bytes())["group_bounds"]
    return openmc.EnergyFilter(sorted(bounds))


def fuel_tallies(
    name: str, filters: list[openmc.Filter], nuclides: list[str], scores: list[str]
) -> list[openmc.Tally]:
    """燃料の核種ごとの反応率の tally（名前 `name`）と、中性子束の tally（`name` + " flux"）。"""
    rates = openmc.Tally(name=name)
    rates.filters = filters
    rates.nuclides = nuclides
    rates.scores = [STORED.get(score, score) for score in scores]
    flux = openmc.Tally(name=f"{name} flux")
    flux.filters = filters
    flux.scores = ["flux"]
    return [rates, flux]


def effective_xs(
    sp: openmc.StatePoint,
    name: str,
    densities: dict[str, float],
    scores: list[str],
    arrange: Callable[[np.ndarray], np.ndarray],
) -> dict[str, dict[str, np.ndarray]]:
    """`fuel_tallies` の tally から、核種と score ごとの実効断面積 (barn) を求める。

    `arrange` は tally の値（フィルタの順）を出力の並びに直す。中性子束が 0 のところは 0 にする。
    """
    phi = arrange(sp.get_tally(name=f"{name} flux").get_values(scores=["flux"]))
    rates = sp.get_tally(name=name)
    xs = {}
    for nuclide, density in densities.items():
        xs[nuclide] = {}
        for score in scores:
            rate = arrange(
                rates.get_values(scores=[STORED.get(score, score)], nuclides=[nuclide])
            )
            xs[nuclide][score] = np.divide(
                rate, density * phi, out=np.zeros_like(rate), where=phi > 0
            )
    return xs
