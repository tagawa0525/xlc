"""2 つのライブラリの核種の散乱行列を群ごとに比べる（ADR-0133 の自由ガスの上限を決める試行）。

使い方: python3 compare_rows.py <ライブラリ A> <ライブラリ B> <核種> <温度 K>

群ごとに、P0 の行の和、群から出る割合（1 − σ(g→g)/行の和）、上の群へ出る割合、P1 の行の和 / P0 の行の和（平均の
余弦）を A と B で出し、B − A の差を並べる。
"""

import sys
from pathlib import Path

import msgpack


def table(path: str, name: str, temperature: float) -> tuple[list, list]:
    library = msgpack.unpackb(Path(path).read_bytes())
    (nuclide,) = [n for n in library["nuclides"] if n["name"] == name]
    (t,) = [
        t for t in nuclide["temperatures"] if abs(t["temperature"] - temperature) < 1e-6
    ]
    return library["group_bounds"], t["scatter"]


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    bounds, a = table(sys.argv[1], sys.argv[3], float(sys.argv[4]))
    bounds_b, b = table(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    if bounds != bounds_b:
        raise SystemExit("group bounds differ")
    print(
        f"{'g':>3} {'low eV':>10} {'P0 sum':>10} {'out A':>9} {'Δout':>10} {'Δup':>10} {'mu A':>9} {'Δmu':>10}"
    )
    for g in range(len(bounds) - 1):
        rows = []
        for m in (a, b):
            p0 = m[0][g]
            total = sum(p0)
            out = 1.0 - p0[g] / total
            up = sum(p0[:g]) / total
            mu = sum(m[1][g]) / total
            rows.append((total, out, up, mu))
        (ta, oa, ua, ma), (_, ob, ub, mb) = rows
        print(
            f"{g + 1:3d} {bounds[g + 1]:10.4g} {ta:10.5g} {oa:9.5f} {ob - oa:+10.2e} {ub - ua:+10.2e} "
            f"{ma:9.5f} {mb - ma:+10.2e}"
        )


if __name__ == "__main__":
    main()
