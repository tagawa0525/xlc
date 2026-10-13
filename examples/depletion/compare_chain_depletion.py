"""xlc の全核種の鎖と OpenMC の鎖（同じ ENDF/B-VIII.0 から作ったもの。ADR-0164）に同じ 1 群の反応率を与え、1 ステップを
次数 48 の CRAM で解いて比べる。鎖の核種の並び、崩壊と反応の遷移、収率が燃焼行列に入る向きを、行列の外から確かめる。

使い方（references の devShell で。xlc の側は先にビルドした crates/depletion の例を呼ぶ）:

    cargo build --release -p xlc-depletion --example deplete_chain
    python3 examples/depletion/compare_chain_depletion.py examples/depletion/chain-cases.toml \\
        examples/depletion/chain.msgpack examples/depletion/openmc-chain.xml target/release/examples/deplete_chain

差は、絶対値の最大を初期の数密度の最大で割ったもの（ADR-0161 の「確認」と同じ）と、初期の最大の 1e-12 倍より大きな
核種の相対差。
"""

import subprocess
import sys

import msgpack
import numpy as np
import tomllib
from compare_chain import REACTIONS, openmc_name
from openmc.deplete import Chain, ReactionRates
from openmc.deplete.cram import CRAM48

# 経路の核種（表に値を出す）
WATCH = [
    "U-239",
    "Np-239",
    "Pu-239",
    "I-135",
    "Xe-135",
    "Pm-149",
    "Sm-149",
    "Cs-135",
    "Xe-136",
    "Sm-150",
]
THERMAL = 0.0253


def xlc_chain(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False)["nuclides"]


def run_xlc(binary, chain_path, xlc_nuclides, case):
    by_name = {n["name"]: n for n in xlc_nuclides}
    lines = ["order 48", f"t {case['t']!r}"]
    lines += [f"n0 {name} {value!r}" for name, value in case["n0"].items()]
    lines += [
        f"reaction {name} {mt} {rate!r}" for name, mt, rate in case.get("reaction", [])
    ]
    for name, rate in case.get("fission", []):
        energies = by_name[name]["fission"]["energies"]
        weights = [1.0 if e == THERMAL else 0.0 for e in energies]
        if sum(weights) != 1.0:
            sys.exit(f"{name} has no {THERMAL} eV yields")
        lines.append(f"fission {name} {rate!r} " + " ".join(repr(w) for w in weights))
    out = subprocess.run(
        [binary, chain_path],
        input="\n".join(lines),
        capture_output=True,
        text=True,
        check=True,
    )
    return {
        name: float(v) for name, v in (line.split() for line in out.stdout.splitlines())
    }


def run_openmc(chain, case):
    names = [n.name for n in chain.nuclides]
    index = {name: i for i, name in enumerate(names)}
    rates = ReactionRates(["1"], names, chain.reactions)
    for name, mt, rate in case.get("reaction", []):
        rates[0, index[openmc_name(name)], rates.index_rx[REACTIONS[mt]]] = rate
    yields = chain.get_default_fission_yields()
    for name, rate in case.get("fission", []):
        o = openmc_name(name)
        rates[0, index[o], rates.index_rx["fission"]] = rate
        yields[o] = dict(chain[o].yield_data[THERMAL].items())
    a = chain.form_matrix(rates[0], yields)
    n0 = np.zeros(len(names))
    for name, value in case["n0"].items():
        n0[index[openmc_name(name)]] = value
    n = CRAM48(a, n0, case["t"])
    return {name: v for name, v in zip(names, n)}


def main():
    cases_path, chain_path, openmc_path, binary = sys.argv[1:5]
    with open(cases_path, "rb") as f:
        cases = tomllib.load(f)["case"]
    nuclides = xlc_chain(chain_path)
    chain = Chain.from_xml(openmc_path)
    for case in cases:
        xlc = run_xlc(binary, chain_path, nuclides, case)
        omc = run_openmc(chain, case)
        scale = max(case["n0"].values())
        pairs = [(name, v, omc[openmc_name(name)]) for name, v in xlc.items()]
        scaled = max(abs(x - o) for _, x, o in pairs) / scale
        big = [
            (abs(x - o) / abs(o), name, x, o)
            for name, x, o in pairs
            if abs(o) > 1e-12 * scale
        ]
        big.sort(reverse=True)
        print(f"## {case['name']}（t = {case['t']:g} s）")
        print(
            f"差の最大 / 初期の最大 {scaled:.3e}、初期の最大の 1e-12 倍を超える核種 {len(big)}、その相対差の最大 {big[0][0]:.3e}（{big[0][1]}）"
        )
        print("| 核種 | xlc | OpenMC | 相対差 |")
        print("| --- | --- | --- | --- |")
        shown = [name for name in WATCH if abs(omc[openmc_name(name)]) > 1e-12 * scale]
        shown += [b[1] for b in big[:5] if b[1] not in shown]
        for name in shown:
            x, o = xlc[name], omc[openmc_name(name)]
            print(f"| {name} | {x:.10e} | {o:.10e} | {abs(x - o) / abs(o):.2e} |")
        print()


if __name__ == "__main__":
    main()
