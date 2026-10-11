"""xlc の格子計算の k∞（lattice の標準出力の `k <値>` の行）と、OpenMC の k（openmc_model.py extract の k.json）を比べる。

使い方: python3 compare_k.py <lattice の出力> <k.json>

差は検証の記録と同じく (k_xlc − k_OpenMC) × 1e5 pcm で、± は OpenMC の標準偏差 × 1e5。
"""

import json
import sys
from pathlib import Path


def lattice_k(path: Path) -> float:
    """lattice の出力の `k <値>` の行の値。行がちょうど 1 つでなければ止める（状態ごとの出力は対象にしない）。"""
    values = [
        float(line.split()[1])
        for line in path.read_text().splitlines()
        if line.startswith("k ")
    ]
    if len(values) != 1:
        raise SystemExit(f"{path}: expected one 'k' line, found {len(values)}")
    return values[0]


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    ours = lattice_k(Path(sys.argv[1]))
    reference = json.loads(Path(sys.argv[2]).read_text())
    k, sigma = reference["k"], reference["k-std"]
    print(
        f"{reference['problem']} xlc {ours:.6f} OpenMC {k:.6f} +- {sigma:.6f} "
        f"diff {(ours - k) * 1e5:+.0f} +- {sigma * 1e5:.0f} pcm"
    )


if __name__ == "__main__":
    main()
