"""C5G7 の 2D の結果を、ベンチマークの参照解（MCNP、NEA/NSC/DOC(2003)16 の Table 3〜5）と比べる。

使い方: python3 compare.py <k> <pins.txt>

pins.txt は lattice のセルごとの核分裂率（51 × 51）。参照解のピン出力は燃料ピン 1056 本
（案内管と核分裂計数管を除く）の平均が 1 なので、同じく燃料ピンだけで規格化し直す。
"""

import sys

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


def main() -> None:
    k = float(sys.argv[1])
    with open(sys.argv[2]) as f:
        rates = [[float(v) for v in line.split()] for line in f if line.strip()]
    if len(rates) != 3 * ASSEMBLY or any(len(row) != 3 * ASSEMBLY for row in rates):
        raise SystemExit("expected a 51 x 51 table")

    def pins(block_row: int, block_col: int) -> list[float]:
        return [
            rates[block_row * ASSEMBLY + r][block_col * ASSEMBLY + c]
            for r in range(ASSEMBLY)
            for c in range(ASSEMBLY)
            if (r, c) not in NON_FUEL
        ]

    blocks = {"inner UO2": [(0, 0)], "MOX": [(0, 1), (1, 0)], "outer UO2": [(1, 1)]}
    fuel = [p for block in [(0, 0), (0, 1), (1, 0), (1, 1)] for p in pins(*block)]
    mean = sum(fuel) / len(fuel)
    fuel = [p / mean for p in fuel]

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
        total = sum(sum(pins(*b)) / mean for b in positions) / len(positions)
        reference = REFERENCE_ASSEMBLIES[name]
        print(
            f"{name:<10} {total:.1f}     reference {reference:.1f}     diff {(total / reference - 1) * 100:+.2f} %"
        )


if __name__ == "__main__":
    main()
