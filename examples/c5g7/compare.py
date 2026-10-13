"""C5G7 の 2D の結果を、ベンチマークの参照解（MCNP、NEA/NSC/DOC(2003)16 の Table 3〜5）と比べる。

使い方: python3 compare.py <lattice の出力> <pins.msgpack> [<OpenMC の pins.msgpack>]

lattice の出力は標準出力を受けたファイル（`k <値>` の行から k を読む）。pins.msgpack は lattice のセルごとの核分裂率
（51 × 51）。参照解のピン出力は燃料ピン 1056 本（案内管と核分裂計数管を除く）の平均が 1 なので、同じく燃料ピンだけで
規格化し直す。OpenMC の pins.msgpack（openmc_model.py extract の出力）を足すと、燃料ピンごとの相対差
（lattice / OpenMC − 1）の二乗平均と、絶対値が最大のものも表示する。
"""

import math
import sys
from pathlib import Path

import msgpack

REFERENCE_K = 1.186550
REFERENCE_MAX = 2.498
# Table 4 の値は有効数字2桁なので、最小ピン出力の差には ±2% の丸めが含まれる
REFERENCE_MIN = 0.23
# 集合体ごとのピン出力の和（内側の UO2、MOX、外側の UO2）
REFERENCE_ASSEMBLIES = {"inner UO2": 492.8, "MOX": 211.7, "outer UO2": 139.8}

ASSEMBLY = 17
# 集合体の中の案内管と核分裂計数管（中心）の位置（行, 列）、0 始まり
NON_FUEL = {
    (r, c)
    for r, cols in {
        2: (5, 8, 11),
        3: (3, 13),
        5: (2, 5, 8, 11, 14),
        8: (2, 5, 8, 11, 14),
        11: (2, 5, 8, 11, 14),
        13: (3, 13),
        14: (5, 8, 11),
    }.items()
    for c in cols
}


def lattice_k(path: Path) -> float:
    """lattice の出力の `k <値>` の行の値。行がちょうど 1 つでなければ止める。"""
    values = [
        float(line.split()[1])
        for line in path.read_text().splitlines()
        if line.startswith("k ")
    ]
    if len(values) != 1:
        raise SystemExit(f"{path}: expected one 'k' line, found {len(values)}")
    return values[0]


def read_rates(path: str) -> list[list[float]]:
    with open(path, "rb") as f:
        data = msgpack.unpack(f)
    if data["format"] != "xlc-pin-powers":
        raise SystemExit(f"not a pin power file: {data['format']}")
    rates = data["rows"]
    if len(rates) != 3 * ASSEMBLY or any(len(row) != 3 * ASSEMBLY for row in rates):
        raise SystemExit("expected a 51 x 51 table")
    return rates


def pins(rates: list[list[float]], block_row: int, block_col: int) -> list[float]:
    """集合体 (block_row, block_col) の燃料ピンの値（案内管と核分裂計数管を除く）。"""
    return [
        rates[block_row * ASSEMBLY + r][block_col * ASSEMBLY + c]
        for r in range(ASSEMBLY)
        for c in range(ASSEMBLY)
        if (r, c) not in NON_FUEL
    ]


def fuel_pins(rates: list[list[float]]) -> list[float]:
    """燃料ピンの値を、燃料ピンの平均で 1 にして並べる。"""
    fuel = [
        p for block in [(0, 0), (0, 1), (1, 0), (1, 1)] for p in pins(rates, *block)
    ]
    mean = sum(fuel) / len(fuel)
    return [p / mean for p in fuel]


def main() -> None:
    if len(sys.argv) not in (3, 4):
        raise SystemExit(__doc__)
    k = lattice_k(Path(sys.argv[1]))
    rates = read_rates(sys.argv[2])
    blocks = {"inner UO2": [(0, 0)], "MOX": [(0, 1), (1, 0)], "outer UO2": [(1, 1)]}
    fuel_all = [
        p for block in [(0, 0), (0, 1), (1, 0), (1, 1)] for p in pins(rates, *block)
    ]
    mean = sum(fuel_all) / len(fuel_all)
    fuel = [p / mean for p in fuel_all]

    print(
        f"k          {k:.6f}  reference {REFERENCE_K:.6f}  diff {(k - REFERENCE_K) * 1e5:+.0f} pcm"
    )
    print(
        f"max pin    {max(fuel):.3f}     reference {REFERENCE_MAX:.3f}     diff {(max(fuel) / REFERENCE_MAX - 1) * 100:+.2f} %"
    )
    print(
        f"min pin    {min(fuel):.3f}     reference {REFERENCE_MIN:.3f}     diff {(min(fuel) / REFERENCE_MIN - 1) * 100:+.2f} %"
    )
    for name, positions in blocks.items():
        total = sum(sum(pins(rates, *b)) / mean for b in positions) / len(positions)
        reference = REFERENCE_ASSEMBLIES[name]
        print(
            f"{name:<10} {total:.1f}     reference {reference:.1f}     diff {(total / reference - 1) * 100:+.2f} %"
        )
    if len(sys.argv) == 4:
        reference = fuel_pins(read_rates(sys.argv[3]))
        diffs = [a / b - 1 for a, b in zip(fuel, reference)]
        rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
        worst = max(diffs, key=abs)
        print(
            f"vs OpenMC  rms {rms * 100:.2f} %  max {worst * 100:+.2f} %  ({len(diffs)} fuel pins)"
        )


if __name__ == "__main__":
    main()
