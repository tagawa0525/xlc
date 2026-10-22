"""#0133 の試行: 新しい燃料の 1C のピンセルを OpenMC で解き、燃料の U-235 と U-238 の核分裂で生まれる中性子の群ごとの割合
（放出のスペクトル χ）を、xlc の多群ライブラリの χ（ADR-0042 の 1 本のベクトル）と比べる。

使い方（references の devShell で）: python3 openmc_chi.py <openmc-1C-rings1.toml> <xlc のライブラリ> <作業ディレクトリ>
    <粒子数> <世代数> <非アクティブ> <スレッド数>
"""

import os
import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "examples"))
sys.path.insert(0, str(ROOT / "examples" / "depletion"))
sys.path.insert(0, str(ROOT / "examples" / "vera"))
import openmc_depletion
from openmc_common import group_filter

NUCLIDES = ["U235", "U238"]
FUEL_TEMPERATURE = 900.0


def main():
    source, library, work = sys.argv[1:4]
    particles, batches, inactive, threads = map(int, sys.argv[4:8])
    run = openmc_depletion.read(Path(source).resolve())
    run.update(particles=particles, batches=batches, inactive=inactive, seed=1)
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    model, fuels = openmc_depletion.build(run)
    groups = group_filter(Path(library).resolve())
    out_filter = openmc.EnergyoutFilter(groups.values)
    tally = openmc.Tally(name="chi")
    tally.filters = [openmc.MaterialFilter([fuels[0]]), out_filter]
    tally.nuclides = NUCLIDES
    tally.scores = ["nu-fission"]
    model.tallies = openmc.Tallies([tally])
    os.environ["OMP_NUM_THREADS"] = str(threads)
    openmc.config["cross_sections"] = run["cross-sections"]
    # 同じ粒子数と世代数で解いた statepoint があれば、解き直さずに読む（集計の名前 chi は読むときに確かめる）
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
        if not np.allclose(out_bins, groups.values, rtol=1e-12):
            raise SystemExit(f"{sp_path}: the chi tally has other group bounds")
        print(f"k {sp.keff}")
        print(
            "群 | 上端 (eV) | "
            + " | ".join(f"{n} OpenMC | xlc | 差 (%)" for n in NUCLIDES)
        )
        rows = {}
        for name in NUCLIDES:
            v = t.get_values(nuclides=[name]).ravel()[::-1]
            rows[name] = v / v.sum()
        chis = {}
        for name in NUCLIDES:
            entry = next(
                n for n in lib["nuclides"] if n["name"] == f"{name[:1]}-{name[1:]}"
            )
            temps = entry["temperatures"]
            # 1C の燃料の温度（900 K）そのものの χ。無ければ止める
            chis[name] = np.array(
                next(x for x in temps if x["temperature"] == FUEL_TEMPERATURE)["chi"]
            )
        for g in range(len(bounds) - 1):
            if bounds[g] < 5e4:
                break
            cells = []
            for name in NUCLIDES:
                o, x = rows[name][g], chis[name][g]
                cells.append(f"{o:.5f} | {x:.5f} | {(x / o - 1) * 100:+.2f}")
            print(f"{g + 1} | {bounds[g]:.4g} | " + " | ".join(cells))
        fast = [g for g in range(len(bounds) - 1) if bounds[g + 1] >= 1.35e6]
        for name in NUCLIDES:
            e = np.sqrt(np.array(bounds[:-1]) * np.maximum(np.array(bounds[1:]), 1e-5))
            print(
                f"{name}: 1.353 MeV より上（群 1〜4）の割合 OpenMC {rows[name][fast].sum():.5f}, xlc {chis[name][fast].sum():.5f}; "
                f"平均の対数エネルギー OpenMC {np.exp((rows[name] * np.log(e)).sum()):.4g}, xlc {np.exp((chis[name] * np.log(e)).sum()):.4g} eV"
            )


main()
