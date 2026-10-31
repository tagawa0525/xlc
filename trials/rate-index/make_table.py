"""試行（#0190）: 途中で条件を変えた燃焼の core で、反応率の Δσ を w でなく L と同じ指標の比で足したときの効きを見るため、
R を持つ公称のテーブルの Δσ を R の点ごとに a / w 倍したテーブルを書く。a と w は、Δσ を w で足した core の途中で条件を変えた燃焼
（core-2A-residual-switch.msgpack）のその点の値（a は熱群の νΣf の線形の補正の比）。w が 0 の点は変えない。core が Δσ を指標の比で
足すのと、この燃焼の R の点でだけ同じ（点の間の内挿と、組成が変われば a も変わることは含まない）。

    python3 make_table.py   # table-rate-index.msgpack を書く
    ../../target/release/core core-rate-index.toml > core-rate-index.out 2> core-rate-index.log
"""

import sys

import msgpack

H = "../../examples/depletion/history"
sys.path.insert(0, H)
from compare_core_history import Table, production, single_node
from compare_residual import coefficients, read


def scale(x, s):
    if isinstance(x, list):
        return [scale(v, s) for v in x]
    if isinstance(x, dict):
        return {k: (v if k == "mt" else scale(v, s)) for k, v in x.items()}
    return None if x is None else x * s


def main():
    table = read(f"{H}/table-residual-nominal-histories.msgpack")
    t = Table(table)
    micro = table["micro"]
    nodes = single_node(
        read(f"{H}/core-2A-residual-switch.msgpack"), "switch", micro["nuclides"]
    )
    rod = table["histories"][0]["rods"][0]
    for n in nodes:
        k = t.residual_index(n["burnup"])
        w = t.weight(n["history"])
        if k is None or w == 0.0:
            continue
        p, b = t.point(n["burnup"], n["moderator"], "switch")
        factors = coefficients(micro["states"][p])
        own = [x - y for x, y in zip(n["densities"], micro["densities"][b])]
        branch = [x - y for x, y in zip(rod["densities"][k], micro["densities"][b])]
        s = production(factors, own) / production(factors, branch) / w
        print(f"{rod['burnup'][k]:g} MWd/kg: a / w = {s:.4f}")
        rod["micro"][k] = [scale(r, s) for r in rod["micro"][k]]
    with open("table-rate-index.msgpack", "wb") as f:
        msgpack.pack(table, f)


if __name__ == "__main__":
    main()
