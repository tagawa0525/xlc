"""核種ごとの xs の出力（xs-nuclides/<核種>.msgpack）の燃料の温度の点で、ライブラリの断面積が √T の滑らかな曲線から外れる
大きさを見る（#0274）。

使い方: python3 library_residuals.py <核種の出力>...

燃料の温度の軸（ADR-0239）の 7 点ごとに、ポイントワイズの断面積（ポイントワイズの範囲のうち 0.625 eV より上。範囲の下端が
0.625 eV より上なら下端から）を、狭い共鳴の
近似の中性子束 1/(E (σt + σ0)) で平均した値（σ0 は 1e10 と 50 b）と、Bondarenko の表の 0.625 eV より上の群の値のレサジーの重みの
和（σ0 は無限希釈と 100 b）を、√T の 3 次式に最小 2 乗で当てはめた残差（点の値の平均で割った値、1e-5 単位）を出す。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

FUEL = [560.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0]


def residual(y):
    g = np.sqrt(FUEL)
    return (y - np.polyval(np.polyfit(g, y, 3), g)) / y.mean() * 1e5


def trapezoid(e, f):
    return np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(e))


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    print("fuel (K): " + " ".join(f"{t:>8g}" for t in FUEL))
    for path in sys.argv[1:]:
        d = msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)
        (nuclide,) = d["nuclides"]
        tables = {t["temperature"]: t for t in nuclide["temperatures"]}
        missing = [t for t in FUEL if t not in tables]
        if missing:
            raise SystemExit(f"{path}: no table at {missing} K")
        low, high = d["pointwise_range"]
        bounds = np.array(d["group_bounds"])
        fast = bounds[1:] >= 0.625
        lethargy = np.log(bounds[:-1] / bounds[1:])[fast]
        print(f"{nuclide['name']}（ポイントワイズ {max(low, 0.625):g}〜{high:g} eV）")
        for reaction in ("capture", "fission"):
            for sigma0 in (1e10, 50.0):
                y = []
                for t in FUEL:
                    p = tables[t]["pointwise"]
                    e = np.array(p["energy"])
                    m = (e >= max(low, 0.625)) & (e <= high)
                    e = e[m]
                    flux = 1 / (e * (np.array(p["total"])[m] + sigma0))
                    y.append(
                        trapezoid(e, np.array(p[reaction])[m] * flux)
                        / trapezoid(e, flux)
                    )
                y = np.array(y)
                if y.mean() == 0:
                    continue
                print(
                    f"  ポイントワイズの {reaction}、σ0 = {sigma0:g} b: "
                    + " ".join(f"{x:+8.2f}" for x in residual(y))
                )
            for index in (0, d["sigma0"].index(100.0)):
                y = np.array(
                    [
                        np.sum(
                            np.array(tables[t]["shielded"][index][reaction])[fast]
                            * lethargy
                        )
                        for t in FUEL
                    ]
                )
                if y.mean() == 0:
                    continue
                print(
                    f"  Bondarenko の {reaction}、σ0 = {d['sigma0'][index]:g} b: "
                    + " ".join(f"{x:+8.2f}" for x in residual(y))
                )


if __name__ == "__main__":
    main()
