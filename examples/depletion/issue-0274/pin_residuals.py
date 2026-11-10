"""燃料の温度の軸だけを持つテーブルで、2 群の定数と k∞ が √T の滑らかな曲線から外れる大きさを見る（#0274）。

使い方: python3 pin_residuals.py <テーブル> [<比べるテーブル>]

燃料の温度の点ごとに、次を出す:

- 高速群と熱群の吸収と νΣf、k∞ を √T の 3 次式に最小 2 乗で当てはめた残差（点の値の平均で割った値、1e-5 単位。点が 7 より
  少なければ 2 次式）
- その点を除いた点から √T の線形と、近い 4 点のラグランジュで内挿した k∞ の反応度の誤差（pcm。両端の点は除く）
- 比べるテーブルを渡したときは、温度ごとの 400 kT（分離共鳴を持つ核種の自由ガスの既定の上限。ADR-0204、ADR-0209）と、
  高速群の吸収の相対差（1e-5）と、比べるテーブルの反応度からテーブルの反応度を引いた差 1/k(テーブル) − 1/k(比べる)（pcm）
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

BOLTZMANN = 8.617333262e-5  # eV/K


def load(path):
    t = msgpack.unpackb(Path(path).read_bytes())
    axes = {a["name"]: a["points"] for a in t["axes"]}
    if any(len(p) != 1 for name, p in axes.items() if name != "fuel"):
        raise SystemExit(f"{path}: every axis but fuel must hold one point")
    return np.array(axes["fuel"]), t["states"]


def lagrange(x, y, x0):
    return sum(
        y[i] * np.prod([(x0 - x[j]) / (x[i] - x[j]) for j in range(len(x)) if j != i])
        for i in range(len(x))
    )


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    temps, states = load(sys.argv[1])
    g = np.sqrt(temps)
    degree = 3 if len(temps) >= 7 else 2
    print("fuel (K): " + " ".join(f"{t:>8g}" for t in temps))
    rows = [(f"absorption {n + 1} 群", "absorption", n) for n in range(2)]
    rows += [(f"nu_fission {n + 1} 群", "nu_fission", n) for n in range(2)]
    rows.append(("k_inf", "k_inf", None))
    for label, key, group in rows:
        y = np.array([s[key] if group is None else s[key][group] for s in states])
        r = (y - np.polyval(np.polyfit(g, y, degree), g)) / y.mean() * 1e5
        print(
            f"{label}（√T の {degree} 次式からの残差、1e-5）: "
            + " ".join(f"{x:+8.2f}" for x in r)
        )
    k = np.array([s["k_inf"] for s in states])
    linear, four = [], []
    for i in range(1, len(temps) - 1):
        rest = [j for j in range(len(temps)) if j != i]
        linear.append((1 / np.interp(g[i], g[rest], k[rest]) - 1 / k[i]) * 1e5)
        near = sorted(rest, key=lambda j: abs(g[j] - g[i]))[:4]
        four.append((1 / lagrange(g[near], k[near], g[i]) - 1 / k[i]) * 1e5)
    print(
        "1 点を除いた √T の線形の内挿の誤差（pcm、内側の点）: "
        + " ".join(f"{x:+6.1f}" for x in linear)
    )
    print(
        "1 点を除いた √T の 4 点のラグランジュの内挿の誤差（pcm、内側の点）: "
        + " ".join(f"{x:+6.1f}" for x in four)
    )
    if len(sys.argv) == 3:
        other_temps, other = load(sys.argv[2])
        if not np.array_equal(temps, other_temps):
            raise SystemExit("the two tables must hold the same fuel temperatures")
        print("400 kT (eV): " + " ".join(f"{400 * BOLTZMANN * t:8.2f}" for t in temps))
        a = np.array([s["absorption"][0] for s in states])
        b = np.array([s["absorption"][0] for s in other])
        print(
            "高速群の吸収の差 (比べる − テーブル)/テーブル（1e-5）: "
            + " ".join(f"{x:+8.2f}" for x in (b - a) / a * 1e5)
        )
        kb = np.array([s["k_inf"] for s in other])
        print(
            "反応度の差 1/k(テーブル) − 1/k(比べる)（pcm）: "
            + " ".join(f"{x:+8.2f}" for x in (1 / k - 1 / kb) * 1e5)
        )


if __name__ == "__main__":
    main()
