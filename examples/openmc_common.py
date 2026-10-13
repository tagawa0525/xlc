"""OpenMC の参照解のスクリプト（examples/vera と examples/mox の openmc_model.py）で共通に使う部品。

物質、固有値計算の設定、燃料の核種ごとの実効断面積（燃料の領域の反応率 / (個数密度 × 中性子束)）の
tally と、その statepoint からの取り出し。
"""

from collections.abc import Callable
from pathlib import Path

import msgpack
import numpy as np
import openmc

# tally に書く score の名前のうち、statepoint で別の名前になるもの（弾性散乱は "(n,elastic)" になる）と、
# 書き出す名前と score の名前が違うもの（捕獲は "(n,gamma)"）
STORED = {"elastic": "(n,elastic)", "capture": "(n,gamma)"}


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
    """世代の数と乱数のシード、燃料に一様な初期の中性子源（`half` は外周の箱の半分の幅）。"""
    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.seed = seed
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
        constraints={"fissionable": True},
    )
    return settings


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
