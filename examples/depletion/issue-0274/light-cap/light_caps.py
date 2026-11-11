"""分離共鳴を持たない核種の自由ガスの上限を変えたピンセルを、基準のテーブル（上限なし）と比べる（ADR-0246、#0274）。

使い方: python3 light_caps.py <基準のテーブル> <テーブル>...

テーブルは少数群定数の出力で、状態を 1 つだけ持ち、基準と同じ群に分けたものであること。テーブルごとに、反応度の差
1/k(テーブル) − 1/k(基準)（pcm）、最後の群（熱群）より上の群から熱群への移行（P0 の散乱行列の和）の相対差と、群ごとの吸収の相対差
（どちらも (テーブル − 基準)/基準、1e-5）を出す。名前が l-[<条件>-]<kT の倍数>.msgpack のテーブルは、上限の倍数の順に並べ、
その倍数の燃料（900 K）と減速材（580 K）の上限のエネルギー (eV) も出す（ほかの名前は倍数の後に名前の順）。
"""

import math
import re
import sys
from pathlib import Path

import msgpack

# eV/K（CODATA 2018。lattice の `free_gas::BOLTZMANN` と同じ）
BOLTZMANN = 8.617333262e-5
FUEL = 900.0
MODERATOR = 580.0


def state(path: str) -> dict:
    t = msgpack.unpackb(Path(path).read_bytes())
    if len(t["states"]) != 1:
        raise SystemExit(f"{path}: the table must hold one state")
    return t["states"][0]


def multiple(path: str) -> float | None:
    """l-[<条件>-]<倍数>.msgpack の倍数（inf は無限大）。ほかの名前は None"""
    m = re.fullmatch(
        r"l-(?:[a-z0-9]+-)?([0-9]+(?:\.[0-9]+)?|inf)\.msgpack", Path(path).name
    )
    if m is None:
        return None
    return math.inf if m.group(1) == "inf" else float(m.group(1))


def into_thermal(table: dict) -> float:
    """熱群（最後の群）より上の群から熱群への移行の和。散乱行列は [出る群][入る群]"""
    scatter = table["scatter"]
    return sum(row[-1] for row in scatter[:-1])


def compare(reference: dict, table: dict) -> tuple[float, float, list[float]]:
    rho = (1 / table["k_inf"] - 1 / reference["k_inf"]) * 1e5
    if len(table["absorption"]) != len(reference["absorption"]):
        raise SystemExit("the table and the reference must have the same groups")
    absorption = [
        (a - b) / b * 1e5
        for a, b in zip(table["absorption"], reference["absorption"], strict=True)
    ]
    thermal = (into_thermal(table) / into_thermal(reference) - 1) * 1e5
    return rho, thermal, absorption


def rows(reference_path: str, paths: list[str]) -> list[str]:
    reference = state(reference_path)

    def key(path: str):
        n = multiple(path)
        return (n is None, n if n is not None else 0.0, Path(path).name)

    out = []
    for path in sorted(paths, key=key):
        rho, thermal, absorption = compare(reference, state(path))
        n = multiple(path)
        caps = (
            f"{n * BOLTZMANN * FUEL:9.2f} {n * BOLTZMANN * MODERATOR:9.2f}"
            if n is not None
            else f"{'':9} {'':9}"
        )
        cells = " ".join(f"{a:+9.2f}" for a in absorption)
        out.append(f"{Path(path).name:28} {caps} {rho:+9.2f} {thermal:+9.2f}  {cells}")
    return out


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    print(
        "テーブル  燃料の上限 (eV)  減速材の上限 (eV)  反応度の差 (pcm)  熱群への移行の差 (1e-5)  "
        "群ごとの吸収の差 (1e-5)"
    )
    print("\n".join(rows(sys.argv[1], sys.argv[2:])))


if __name__ == "__main__":
    main()
