"""燃料の温度 900 K の 1 点で、分離共鳴を持つ核種の自由ガスの上限（kT の倍数。ADR-0204）を変えたピンセルを、静止した標的と
比べる（#0274）。

使い方: python3 read_caps.py <静止した標的のテーブル> <上限のテーブル>...

上限のテーブルの名前は h-<kT の倍数か inf>.msgpack（`h-*.msgpack` のように静止した標的のテーブルを含めて渡してもよく、静止した
標的のテーブルは上限の並びから除く。同じ倍数を 2 度渡すと止める）。テーブルごとに、上限のエネルギー（倍数 × 900 K の kT）、静止した標的からの
反応度の差 1/k(静止) − 1/k(上限)（pcm）、高速群の吸収の相対差 (上限 − 静止)/静止（1e-5）を、上限のエネルギーの順に出す。
"""

import math
import sys
from pathlib import Path

import msgpack

BOLTZMANN = 8.617333262e-5  # eV/K
TEMPERATURE = 900.0


def state(path):
    t = msgpack.unpackb(Path(path).read_bytes())
    if len(t["states"]) != 1:
        raise SystemExit(f"{path}: the table must hold one state")
    return t["states"][0]


def multiple(path):
    name = Path(path).name
    if not (name.startswith("h-") and name.endswith(".msgpack")):
        raise SystemExit(f"{path}: the name must be h-<multiple of kT or inf>.msgpack")
    value = name[2 : -len(".msgpack")]
    if value == "inf":
        return math.inf
    try:
        return float(value)
    except ValueError:
        raise SystemExit(
            f"{path}: the name must be h-<multiple of kT or inf>.msgpack"
        ) from None


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    static = state(sys.argv[1])
    paths = [
        p for p in sys.argv[2:] if Path(p).resolve() != Path(sys.argv[1]).resolve()
    ]
    multiples = [multiple(p) for p in paths]
    if len(set(multiples)) != len(multiples):
        raise SystemExit("each multiple of kT must be given once")
    rows = sorted(zip(multiples, map(state, paths)), key=lambda row: row[0])
    print("上限 (kT)  上限 (eV)  反応度の差 (pcm)  高速群の吸収の差 (1e-5)")
    for n, s in rows:
        energy = n * BOLTZMANN * TEMPERATURE
        rho = (1 / static["k_inf"] - 1 / s["k_inf"]) * 1e5
        a = (s["absorption"][0] - static["absorption"][0]) / static["absorption"][0]
        print(f"{n:>9g}  {energy:9.2f}  {rho:+16.2f}  {a * 1e5:+23.2f}")


if __name__ == "__main__":
    main()
