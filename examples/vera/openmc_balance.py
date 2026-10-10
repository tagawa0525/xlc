"""VERA の 1B〜1D を OpenMC の連続エネルギーモードで解き、物質ごと・群ごとの中性子束、正味の吸収、生成を書く。

使い方: python3 openmc_balance.py <cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [非アクティブバッチ数] [乱数のシード]

xlc の `balance` の例（crates/lattice/examples/balance.rs）と同じ形式の balance.msgpack を書く。群構造は xlc の
ライブラリのもので、高いエネルギーから並べる。正味の吸収は、吸収（中性子が出ない反応と核分裂）から、(n,2n)、(n,3n)、
(n,4n) で増える中性子の数を引いたもの。値は生成の和が 1 になるように規格化する（外周が鏡面反射なので、正味の吸収の和が
1/k になる）。物質の名前は xlc の 1B の入力にそろえる。燃料の核種ごとの値も同じ形で `nuclides` に書く（物質の名前/核種）。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_common import eigenvalue_settings, group_filter, material
from openmc_model import (
    COLD,
    FUEL,
    FUEL_TEMPERATURE,
    MODERATOR,
    ZIRCALOY,
    geometry_of,
)

# xlc の 1B の入力の物質の名前
NAMES = {"fuel": "fuel", "zircaloy-4": "zircaloy", "moderator": "water"}
# 正味の吸収で引く、反応ごとに増える中性子の数
MULTIPLYING = {"(n,2n)": 1.0, "(n,3n)": 2.0, "(n,4n)": 3.0}
SCORES = ["absorption", "nu-fission", *MULTIPLYING]


def main() -> None:
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    cross_sections = Path(sys.argv[1]).resolve()
    problem = sys.argv[2]
    if not problem.startswith("1") or problem not in FUEL_TEMPERATURE:
        raise SystemExit(f"{problem}: only the pin cells 1B〜1D")
    work = Path(sys.argv[3]).resolve()
    library = Path(sys.argv[4]).resolve()
    particles = int(sys.argv[5]) if len(sys.argv) > 5 else 1_000_000
    batches = int(sys.argv[6]) if len(sys.argv) > 6 else 250
    inactive = int(sys.argv[7]) if len(sys.argv) > 7 else 50
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    work.mkdir(parents=True, exist_ok=True)

    fuel = material("fuel", FUEL, FUEL_TEMPERATURE[problem])
    clad = material("zircaloy-4", ZIRCALOY, COLD)
    water = material("moderator", MODERATOR, COLD)
    water.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([fuel, clad, water])
    materials.cross_sections = str(cross_sections)
    geometry, half = geometry_of(problem, fuel, clad, water)
    settings = eigenvalue_settings(particles, batches, inactive, seed, half)
    settings.temperature = {"method": "nearest", "tolerance": 1.0}

    by_material = openmc.MaterialFilter([fuel, clad, water])
    energy = group_filter(library)
    # 燃料の核種ごとと、物質ごとの和（"total"）。全核種を数えると、出来事ごとの断面積の計算が重く何倍も遅くなる
    nuclides = sorted(FUEL) + ["total"]
    rates = openmc.Tally(name="rates")
    rates.filters = [by_material, energy]
    rates.nuclides = nuclides
    rates.scores = SCORES
    flux = openmc.Tally(name="flux")
    flux.filters = [by_material, energy]
    flux.scores = ["flux"]
    model = openmc.Model(geometry, materials, settings, openmc.Tallies([rates, flux]))
    statepoint = model.run(cwd=work, output=False)

    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff
        tally = sp.get_tally(name="rates")
        phi = sp.get_tally(name="flux").get_values(scores=["flux"]).reshape(3, -1)
        groups = phi.shape[1]

        def values(score: str, nuclide: str) -> np.ndarray:
            v = tally.get_values(scores=[score], nuclides=[nuclide])
            return v.reshape(3, groups)

        absorption = {n: values("absorption", n) for n in nuclides}
        for n in nuclides:
            for score, extra in MULTIPLYING.items():
                absorption[n] = absorption[n] - extra * values(score, n)
        production = {n: values("nu-fission", n) for n in nuclides}
        total = production["total"].sum()
        order = [fuel, clad, water]
        out = {
            "format": "xlc-balance",
            "k": k.nominal_value,
            "materials": {},
            "nuclides": {},
        }
        for i, m in enumerate(order):
            name = NAMES[m.name]
            # 群は高いエネルギーから
            out["materials"][name] = {
                "flux": (phi[i][::-1] / total).tolist(),
                "absorption": (absorption["total"][i][::-1] / total).tolist(),
                "production": (production["total"][i][::-1] / total).tolist(),
            }
            for n in FUEL:
                if n in m.get_nuclides():
                    out["nuclides"][f"{name}/{n}"] = {
                        "absorption": (absorption[n][i][::-1] / total).tolist(),
                        "production": (production[n][i][::-1] / total).tolist(),
                    }
        (work / "balance.msgpack").write_bytes(msgpack.packb(out))
        absorbed = sum(sum(v["absorption"]) for v in out["materials"].values())
        print(
            f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f} 1/absorption {1 / absorbed:.6f}"
        )


if __name__ == "__main__":
    main()
