"""NJOY PURR の出力から、エネルギーごとの直接サンプリングの Bondarenko の表を、平均に対する比で書く。

使い方: python3 purr_table.py <PURR の output>
出力: energy sigma0 capture_ratio total_ratio（タブ区切り。比は σ0 = 1e10 の行に対する）

σ0 の数と値は出力の前置き（sigma zero values）から読み、表の行の数、列の名前、σ0 がそれと合わなければ止める。
温度は 1 つだけを扱う。
"""

import re
import sys

# 列の見出し（データの列は temp、sig0、p0 total、elastic、fission、capture、p1 total の 7 つ）
HEADER = ["temp", "sig0", "p0", "total", "elastic", "fission", "capture", "p1", "total"]
INFINITY = 1.0e10  # PURR は σ0 の表の infinity を 1e10 と書く


def fail(message: str) -> None:
    sys.exit(f"purr_table.py: {message}")


def after_leader(line: str) -> list[str]:
    """点の並び（.....）のあとの値"""
    return re.sub(r"^.*?\.{3,}\s*", "", line).split()


def settings(lines: list[str]) -> list[float]:
    """前置きの σ0 の値（infinity は INFINITY）。温度が 1 つでなければ止める"""
    for n, line in enumerate(lines):
        if line.startswith(" temperatures ...") and len(after_leader(line)) != 1:
            fail("only one temperature is supported")
        if line.startswith(" sigma zero values ..."):
            values = after_leader(line)
            for row in lines[n + 1 :]:
                if row.lstrip() == row or not row.strip() or "...." in row:
                    break
                values.append(row.split()[0])
            return [INFINITY if v == "infinity" else float(v) for v in values]
    fail("no 'sigma zero values' in the output")
    return []


def main() -> None:
    with open(sys.argv[1]) as f:
        lines = f.read().splitlines()
    sigma0 = settings(lines)
    energy = None
    print("energy\tsigma0\tcapture_ratio\ttotal_ratio")
    for n, line in enumerate(lines):
        if line.startswith("e= "):
            energy = float(line.split()[1])
        elif (
            "bondarenko cross sections by direct sampling" in line
            and energy is not None
        ):
            if lines[n + 1].split() != HEADER:
                fail(f"unexpected columns at {energy:e} eV: {lines[n + 1].split()}")
            rows = [
                [float(x) for x in row.split()]
                for row in lines[n + 2 : n + 2 + len(sigma0)]
            ]
            if (
                any(len(v) != 7 for v in rows)
                or [v[1] for v in rows] != sigma0
                or (lines[n + 2 + len(sigma0)].strip() != "")
            ):
                fail(
                    f"the table at {energy:e} eV does not have the sigma zeros {sigma0}"
                )
            mean = rows[0]
            for v in rows:
                print(
                    f"{energy:e}\t{v[1]:e}\t{v[5] / mean[5]:.6f}\t{v[2] / mean[2]:.6f}"
                )
            energy = None


main()
