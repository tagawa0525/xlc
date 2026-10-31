"""試行（#0190）: 途中で条件を変えた燃焼の 3 つの core の、参照との差と数密度の効き（compare_switch.py と同じ式）。

- R を w で足した前の core（main f4c5a4f97 の core-2A-residual-switch.msgpack、md5 14c85ae0。DVC のキャッシュから読む）
- U を w、L を熱群の νΣf の補正の比、Δσ を w で足した core（core-2A-residual-switch.msgpack）
- 同じで、Δσ を R の点で a / w 倍したテーブルの core（make_table.py、core-rate-index.toml）

    python3 compare.py   # ワークツリーの直下の references/dvc-cache を読む
"""

import sys

H = "../../examples/depletion/history"
sys.path.insert(0, H)
from compare_core_history import Table, single_node
from compare_residual import coefficients, read, rho
from compare_switch import at, composition, composition_nuclides, reference_states

RUNS = (
    (
        "R を w",
        "../../references/dvc-cache/files/md5/14/c85ae05f529abf591f0f810d6a6332",
    ),
    ("U を w、L を指標", f"{H}/core-2A-residual-switch.msgpack"),
    ("同じで Δσ を指標", "core-rate-index.msgpack"),
)


def main():
    table = read(f"{H}/table-residual-nominal-histories.msgpack")
    reference = read(f"{H}/table-switch-density-branch.msgpack")
    t = Table(table)
    nuclides = table["micro"]["nuclides"]
    for label, path in RUNS:
        nodes = single_node(read(path), label, nuclides)
        parts = []
        for u, b, p in reference_states(reference)[1:]:
            n = at(nodes, u, label)
            state = reference["states"][p]
            point, _ = t.point(u, n["moderator"], label)
            coeff = coefficients(table["micro"]["states"][point])
            n_ref = reference["micro"]["densities"][b]
            difference = rho(n["k"]) - rho(state["k_inf"])
            effect = composition(state, coeff, n_ref, n)
            parts.append(f"{u:g}: {difference:+.3f}（数密度 {effect:+.3f}）")
            last = composition_nuclides(state, coeff, n_ref, n, nuclides)
        print(f"{label}: " + "、".join(parts))
        print(
            "  60 MWd/kg の核種: "
            + "、".join(f"{x['nuclide']} {x['pcm']:+.2f}" for x in last)
        )


if __name__ == "__main__":
    main()
