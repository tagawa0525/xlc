"""区間に分けた燃焼（openmc_depletion_resume.py）が、分けない燃焼（openmc_depletion.py）と統計の範囲で同じ k と数密度を出すことの確認（ADR-0174）。

使い方（references の devShell で。medium.toml の種を 2 に変えた入力でもう 1 本の分けない実行を流して比べる）:

    python3 compare_split.py <分けない実行の msgpack> <区間に分けた実行の msgpack> <種を変えた分けない実行の msgpack>
"""

import math
import statistics
import sys

import msgpack


def load(path: str) -> dict:
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read())


def last_density(result: dict, nuclide: str) -> float:
    return result["rings"][0]["densities"][nuclide][-1]


def main(whole_path: str, split_path: str, other_path: str) -> None:
    whole, split, other = load(whole_path), load(split_path), load(other_path)
    for i in range(len(whole["k"])):
        sigma_split = math.hypot(whole["k-std"][i], split["k-std"][i])
        sigma_other = math.hypot(whole["k-std"][i], other["k-std"][i])
        d_split = split["k"][i] - whole["k"][i]
        d_other = other["k"][i] - whole["k"][i]
        print(
            f"point {i}: whole {whole['k'][i]:.5f} +- {whole['k-std'][i]:.5f}, "
            f"split-whole {d_split:+.5f} ({d_split / sigma_split:+.2f} sigma), "
            f"other seed-whole {d_other:+.5f} ({d_other / sigma_other:+.2f} sigma)"
        )
    rows = []
    for nuclide in whole["rings"][0]["densities"]:
        n = last_density(whole, nuclide)
        if n > 1e-9:
            rows.append(
                (
                    nuclide,
                    n,
                    last_density(split, nuclide) / n - 1,
                    last_density(other, nuclide) / n - 1,
                )
            )
    rows.sort(key=lambda r: -r[1])
    print(f"{len(rows)} nuclides above 1e-9 /(b cm) at the last point")
    for nuclide, n, d_split, d_other in rows[:8]:
        print(f"{nuclide:8s} {n:.4e} split {d_split:+.2e} other seed {d_other:+.2e}")
    for label, column in (("split", 2), ("other seed", 3)):
        values = [abs(r[column]) for r in rows]
        print(
            f"|{label}/whole-1|: median {statistics.median(values):.2e}, "
            f"max {max(values):.2e}"
        )


if __name__ == "__main__":
    main(*sys.argv[1:4])
