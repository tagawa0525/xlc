"""xlc と OpenMC の、集合体（VERA 2B〜2D）の燃料ピンごとの実効断面積を比べる（自己遮蔽の位置による誤差）。

使い方: python3 compare_pin_xs.py <xlc の pin-xs.msgpack> <OpenMC の pin-xs.msgpack> [核種]

xlc 側は crates/lattice/examples/pin_xs.rs、OpenMC 側は openmc_model.py で作る。核種（省略すると U238）の
非弾性（全 − 弾性。共鳴の群では吸収に等しい）の相対差 (xlc / OpenMC − 1) を、共鳴の群（WIMS 69 群の 15〜27 群）
ごとに、案内管・計装管の隣（上下左右）のピンとそれ以外のピンの平均で表示する。
"""

import sys
from pathlib import Path

import msgpack

# WIMS 69 群の分離共鳴の群（1 始まり）
RESONANCE = range(15, 28)


def read(path: str, nuclide: str) -> list[list[list[float]]]:
    rates = msgpack.unpackb(Path(path).read_bytes())["nuclides"][nuclide]
    return [
        [[t - e for t, e in zip(total, elastic)] for total, elastic in zip(*rows)]
        for rows in zip(rates["total"], rates["elastic"])
    ]


def main() -> None:
    nuclide = sys.argv[3] if len(sys.argv) > 3 else "U238"
    ours, reference = read(sys.argv[1], nuclide), read(sys.argv[2], nuclide)
    n = len(reference)
    if len(ours) != n or any(len(row) != n for row in ours + reference):
        raise SystemExit("the two maps have different shapes")
    if {len(g) for row in ours + reference for g in row} != {len(reference[0][0])}:
        raise SystemExit("the two maps have different numbers of groups")
    fuel = {(i, j) for i in range(n) for j in range(n) if any(reference[i][j])}
    if fuel != {(i, j) for i in range(n) for j in range(n) if any(ours[i][j])}:
        raise SystemExit("the two maps have fuel in different cells")
    holes = {(i, j) for i in range(n) for j in range(n)} - fuel
    beside = {
        (i, j)
        for (i, j) in fuel
        if any(
            (i + di, j + dj) in holes for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1))
        )
    }
    print(f"{nuclide} nonelastic, xlc / OpenMC - 1 (%)")
    print(
        f"group  beside tubes ({len(beside)})  others ({len(fuel) - len(beside)})  worst"
    )
    for group in RESONANCE:
        g = group - 1

        def diff(i: int, j: int, g: int = g) -> float:
            return ours[i][j][g] / reference[i][j][g] - 1

        near = [diff(i, j) for i, j in beside]
        far = [diff(i, j) for i, j in fuel - beside]
        worst = max(fuel, key=lambda p: abs(diff(*p)))
        print(
            f"{group:5d}  {100 * sum(near) / len(near):+19.2f}  "
            f"{100 * sum(far) / len(far):+11.2f}  {100 * diff(*worst):+.2f} at {worst}"
        )


if __name__ == "__main__":
    main()
