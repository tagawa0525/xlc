"""#0133 の試行: openmc_pin_rates.py（と openmc_pin_rates_2a.py）が解いた 1C のピンセルの OpenMC の入力（作業ディレクトリの model.xml）に、
燃料の核種ごとの群ごとの (n,2n)、(n,3n)、(n,4n) の反応率の tally を足して解き直し、`swap_absorption` の例の `--excess` が読む形で書く。

使い方（references の devShell で）: python3 openmc_pin_excess.py <openmc-1C-rings1.toml> <openmc_pin_rates.py の作業ディレクトリ>
    <作業ディレクトリ> <スレッド数>

粒子数、世代数、シードは model.xml のまま（tally は乱数を使わないので、同じ OpenMC とライブラリなら元の実行と同じ粒子の歴史になる。
`same_histories` は k が元の実行とビットで同じか（同じ歴史でも、スレッドの tally の和の順で k が 1e-13 ほど違いうるので、false なら
元の k との差を見る）。歴史が違っても、反応率と中性子束は同じ実行のものなので、その比の断面積は使える）。燃料の物質のフィルタ、群のフィルタ、核種は元の `rates` の tally と同じものを使う（群のフィルタは、448 の版が 10 MeV より
上を別の区間にし、459 の版が群 1 の上端を広げたもの。別の区間は群 1 に足し、標準偏差は 2 つの区間を独立とみて足すので、区間の間の
相関の分だけ小さく出うる）。出力（作業ディレクトリの openmc-excess.msgpack）は `format`、`k`、`k_std`、`same_histories`、`bounds`、`flux`（群ごとの燃料の中性子束、
降順の群）、`densities`、`rates`（xlc の核種名から、反応ごとの `mean` と `std`）で、`bounds`、`densities` は元の
openmc-rates.msgpack の値、`flux` はこの実行の値。
"""

import os
import re
import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
import openmc_depletion

SCORES = ["(n,2n)", "(n,3n)", "(n,4n)"]


def xlc_name(name):
    symbol, mass, isomer = re.fullmatch(r"([A-Z][a-z]?)(\d+)(_m1)?", name).groups()
    return f"{symbol}-{mass}" + ("m" if isomer else "")


def high_first(values, bins_above):
    """昇順の区間の値を、高いエネルギーからの群の順にし、10 MeV より上の区間（`bins_above` 個。0 なら群 1 の上端を広げた tally）を
    群 1 に足す"""
    v = np.asarray(values).ravel()[::-1]
    return np.concatenate([[v[: bins_above + 1].sum()], v[bins_above + 1 :]]).tolist()


def high_first_std(values, bins_above):
    v = np.asarray(values).ravel()[::-1]
    return np.concatenate(
        [[np.sqrt((v[: bins_above + 1] ** 2).sum())], v[bins_above + 1 :]]
    ).tolist()


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    source, rates_dir, work = sys.argv[1:4]
    threads = int(sys.argv[4])
    run = openmc_depletion.read(Path(source).resolve())
    rates_dir = Path(rates_dir).resolve()
    original = msgpack.unpackb((rates_dir / "openmc-rates.msgpack").read_bytes())
    model = openmc.Model.from_model_xml(rates_dir / "model.xml")
    by_name = {t.name: t for t in model.tallies}
    rates = by_name["rates"]
    excess = openmc.Tally(name="excess")
    excess.filters = rates.filters
    excess.nuclides = rates.nuclides
    excess.scores = SCORES
    model.tallies.append(excess)
    # openmc_pin_rates.py の版で、10 MeV より上を 1 つの区間に足した tally（Slurm 448）と、群 1 の上端を広げた tally（459）がある
    (energy,) = [f for f in rates.filters if isinstance(f, openmc.EnergyFilter)]
    edges = list(energy.values)
    bins_above = len(edges) - len(original["bounds"])
    if bins_above not in (0, 1) or edges[-1] <= original["bounds"][0]:
        raise SystemExit(f"the energy filter {edges[-3:]} is not that of the groups")
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    openmc.config["cross_sections"] = run["cross-sections"]
    sp_path = model.run(cwd=work, output=False)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        phi = sp.get_tally(name="flux").get_values(scores=["flux"]).ravel()
        t = sp.get_tally(name="excess")
        out = {
            "format": "xlc-openmc-excess",
            "k": k.nominal_value,
            "k_std": k.std_dev,
            "same_histories": k.nominal_value == original["k"],
            "flux": high_first(phi, bins_above),
            "bounds": original["bounds"],
            "densities": original["densities"],
            "rates": {
                xlc_name(name): {
                    s: {
                        "mean": high_first(
                            t.get_values(scores=[s], nuclides=[name]), bins_above
                        ),
                        "std": high_first_std(
                            t.get_values(scores=[s], nuclides=[name], value="std_dev"),
                            bins_above,
                        ),
                    }
                    for s in SCORES
                }
                for name in rates.nuclides
            },
        }
    (work / "openmc-excess.msgpack").write_bytes(msgpack.packb(out))
    print(
        f"k {k} (the original run {original['k']:.6f}, same histories {out['same_histories']})"
    )


main()
