"""#0133 の試行: 1C のピンセルを OpenMC で解き、燃料の核種の核分裂で生まれる中性子の群ごとの割合（放出のスペクトル χ）を、
xlc の多群ライブラリの χ（ADR-0186 の 1 本のベクトル）と比べ、`swap_absorption` の例の `--chi` が読む形で書く。

使い方（references の devShell で）: python3 openmc_chi.py <openmc-1C-rings1.toml> <xlc のライブラリ> <作業ディレクトリ>
    <粒子数> <世代数> <非アクティブ> <スレッド数> <lattice.toml> <出力 (msgpack)>

燃料の組成を lattice.toml の `[material.fuel]` の `nuclides` に置き換え（openmc_pin_transfer.py と同じ）、組成のアクチニドの
すべてを数える。
出力は `format`、`k`、`k_std`、`bounds`（群の境界 (eV)、降順。群 1 の上端は 10 MeV より上に広げたもの）、`chi`（xlc の核種名から
群ごとの割合）、`chi_std`（その標準偏差）、`share`（数えた核種の生成に占める割合）。生成の無い核種は書かない（xlc は
ライブラリの χ のまま）。
"""

import os
import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc
import tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
sys.path.insert(0, str(ROOT / "examples" / "vera"))
import openmc_depletion
from compare_depletion import openmc_name
from openmc_common import group_filter

ACTINIDES = {"Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf"}
FUEL_TEMPERATURE = 900.0
ABOVE = 1.0e9


def main():
    if len(sys.argv) != 10:
        raise SystemExit(__doc__)
    source, library, work = sys.argv[1:4]
    particles, batches, inactive, threads = map(int, sys.argv[4:8])
    lattice, out_path = Path(sys.argv[8]), Path(sys.argv[9])
    run = openmc_depletion.read(Path(source).resolve())
    run.update(particles=particles, batches=batches, inactive=inactive, seed=1)
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    model, fuels = openmc_depletion.build(run)
    if run["rings"] != 1:
        raise SystemExit(f"the model must have one fuel ring, got {run['rings']}")
    composition = tomllib.loads(lattice.read_text())["material"]["fuel"]["nuclides"]
    (fuel,) = fuels
    for name in list(fuel.get_nuclides()):
        fuel.remove_nuclide(name)
    for name, n in composition.items():
        fuel.add_nuclide(openmc_name(name), n)
    fuel.set_density("sum")
    names = [n for n in composition if n.split("-")[0] in ACTINIDES]
    nuclides = [openmc_name(n) for n in names]
    groups = group_filter(Path(library).resolve())
    # 群構造の上端（10 MeV）より上に生まれた中性子は、xlc の χ と同じく群 1 に数える（ライブラリの χ は群構造の外に出る分を
    # 端の群に入れる。outgoing_group_integrals）。群 1 の上端を、評価済みデータの上限より上に広げる
    out_bounds = list(groups.values)[:-1] + [ABOVE]
    out_filter = openmc.EnergyoutFilter(out_bounds)
    tally = openmc.Tally(name="chi")
    tally.filters = [openmc.MaterialFilter([fuels[0]]), out_filter]
    tally.nuclides = nuclides
    tally.scores = ["nu-fission"]
    model.tallies = openmc.Tallies([tally])
    os.environ["OMP_NUM_THREADS"] = str(threads)
    openmc.config["cross_sections"] = run["cross-sections"]
    # 同じ粒子数と世代数で解いた statepoint があれば、解き直さずに読む（集計の名前 chi と核種は読むときに確かめる）
    done = work / f"statepoint.{batches}.h5"
    sp_path = done if done.exists() else model.run(cwd=work, output=False)
    lib = msgpack.unpackb(Path(library).read_bytes(), strict_map_key=False)
    bounds = lib["group_bounds"]
    with openmc.StatePoint(sp_path) as sp:
        if sp.n_particles != particles or sp.n_batches != batches:
            raise SystemExit(
                f"{sp_path}: {sp.n_particles} particles x {sp.n_batches} batches, "
                f"not {particles} x {batches}"
            )
        t = sp.get_tally(name="chi")
        out_bins = t.find_filter(openmc.EnergyoutFilter).values
        if not np.allclose(out_bins, out_bounds, rtol=1e-12):
            raise SystemExit(f"{sp_path}: the chi tally has other group bounds")
        if list(t.nuclides) != nuclides:
            raise SystemExit(f"{sp_path}: the chi tally counts {list(t.nuclides)}")
        k = sp.keff
        print(f"k {k}")
        produced = {}
        rows = {}
        stds = {}
        for name, nuclide in zip(names, nuclides):
            v = t.get_values(nuclides=[nuclide]).ravel()[::-1]
            s = t.get_values(nuclides=[nuclide], value="std_dev").ravel()[::-1]
            if v.sum() > 0:
                produced[name] = v.sum()
                rows[name] = v / v.sum()
                stds[name] = s / v.sum()
            else:
                print(f"{name}: no neutrons from fission (the library chi stays)")
        total = sum(produced.values())
        chis = {}
        for name in rows:
            entry = next(n for n in lib["nuclides"] if n["name"] == name)
            # 1C の燃料の温度（900 K）そのものの χ。無ければ止める
            chis[name] = np.array(
                next(
                    x
                    for x in entry["temperatures"]
                    if x["temperature"] == FUEL_TEMPERATURE
                )["chi"]
            )
        e = np.sqrt(np.array(bounds[:-1]) * np.maximum(np.array(bounds[1:]), 1e-5))
        fast = [g for g in range(len(bounds) - 1) if bounds[g + 1] >= 1.35e6]
        for name, row in rows.items():
            print(
                f"{name}: 生成の割合 {produced[name] / total:.5f}; 1.353 MeV より上（群 1〜4）の割合 OpenMC "
                f"{row[fast].sum():.5f}, xlc {chis[name][fast].sum():.5f}; 平均の対数エネルギー OpenMC "
                f"{np.exp((row * np.log(e)).sum()):.4g}, xlc {np.exp((chis[name] * np.log(e)).sum()):.4g} eV"
            )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(
            msgpack.packb(
                {
                    "format": "xlc-openmc-chi",
                    "k": k.nominal_value,
                    "k_std": k.std_dev,
                    "bounds": sorted(out_bounds, reverse=True),
                    "chi": {n: rows[n].tolist() for n in rows},
                    "chi_std": {n: stds[n].tolist() for n in rows},
                    "share": {n: produced[n] / total for n in rows},
                }
            )
        )


main()
