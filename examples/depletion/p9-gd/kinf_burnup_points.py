"""燃焼の履歴の k∞（lattice の標準出力の `burnup <B> k <k>` の行）を、ADR-0199 の 27 点から ADR-0192 の 3 次 Hermite で内挿し、
履歴のほかの点の Δρ = (1/k_履歴 − 1/k_内挿) × 1e5 pcm を |Δρ| の大きい順に出す（#0176。基準条件の燃焼点が Gd の燃え尽きの付近で
足りるかの見込み。測定でない）。

使い方（リポジトリの直下から）: python3 examples/depletion/p9-gd/kinf_burnup_points.py <lattice の標準出力>

k∞ の内挿で、core が内挿する少数群定数の内挿の誤差ではない（ADR-0192 の比べ方には、履歴のすべての点で解いた基準条件の
テーブルが要る）。
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from compare_burnup_interpolation import interpolate

# ADR-0199 の 27 点（../p9-2A/branch-2A.toml の `burnup`）
POINTS = [
    0, 0.1, 0.25, 0.5, 0.75, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 24,
    28, 32, 36, 40, 48, 56, 60,
]  # fmt: skip


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    burnup, k = [], []
    for line in Path(sys.argv[1]).read_text().splitlines():
        words = line.split()
        if len(words) != 4 or words[0] != "burnup" or words[2] != "k":
            raise SystemExit(f"not a lattice k∞ line: {line!r}")
        burnup.append(float(words[1]))
        k.append(float(words[3]))
    burnup, k = np.array(burnup), np.array(k)
    chosen = [i for i, b in enumerate(burnup) if b in POINTS]
    if len(chosen) != len(POINTS):
        missing = sorted(set(POINTS) - set(burnup.tolist()))
        raise SystemExit(f"the history has no burnups {missing}")
    rest = [i for i in range(len(burnup)) if i not in chosen]
    values = interpolate("hermite", burnup[chosen], k[chosen][:, None], burnup[rest])
    rho = (1 / k[rest] - 1 / values[:, 0]) * 1e5
    for n in np.argsort(-abs(rho)):
        print(f"{burnup[rest][n]:g} {rho[n]:+.1f}")


if __name__ == "__main__":
    main()
