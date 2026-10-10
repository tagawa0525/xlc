"""VERA の 1B〜1D を OpenMC の連続エネルギーモードで解き、物質ごと・群ごとの中性子束、正味の吸収、生成を書く。

使い方: python3 openmc_balance.py <cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

xlc の `balance` の例（crates/lattice/examples/balance.rs）と同じ形式の balance.msgpack を書く。群構造は xlc の
ライブラリのもので、高いエネルギーから並べる。正味の吸収は、xlc（全断面積から、中性子を出すすべての反応の放出の数を
含む散乱行列の行の和を引く）と同じ量にするため、吸収（中性子が出ない反応と核分裂）から、2 個以上の中性子を出す反応の
増える分（放出の数 − 1）を引く。値はすべてのエネルギーの生成の和が 1 になるように規格化し、群の範囲の外（10 MeV より上）の正味の吸収と生成を
`outside` に書く（外周が鏡面反射なので、正味の吸収の和が 1/k に
なる）。各量の標準偏差を `*_std` に書く（規格化の生成の誤差は含めない）。確認のため、全断面積から中性子の放出の数で
数えた散乱（nu-scatter、出来事ごとの推定）を引いた正味の吸収も `absorption_by_scatter` に書く。物質の名前は xlc の 1B の
入力にそろえる。
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
# 2 個以上の中性子を出す反応（ENDF の MT 11、16、17、24、25、30、37、41、42）と、増える中性子の数
MULTIPLYING = {
    "(n,2nd)": 1.0,
    "(n,2n)": 1.0,
    "(n,3n)": 2.0,
    "(n,2na)": 1.0,
    "(n,3na)": 2.0,
    "(n,2n2a)": 1.0,
    "(n,4n)": 3.0,
    "(n,2np)": 1.0,
    "(n,3np)": 2.0,
}
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

    order = [fuel, clad, water]
    filters = [openmc.MaterialFilter(order), group_filter(library)]
    rates = openmc.Tally(name="rates")
    rates.filters = filters
    rates.scores = SCORES
    flux = openmc.Tally(name="flux")
    flux.filters = filters
    flux.scores = ["flux"]
    # 出来事ごとの推定が要る score は別の tally にする（同じ tally だと全体が出来事ごとの推定になる）
    scattering = openmc.Tally(name="scattering")
    scattering.filters = filters
    scattering.scores = ["total", "nu-scatter"]
    # 群のエネルギーの範囲の外（ライブラリの上の境 10 MeV より上など）も含めた和。釣り合いが閉じるかの確認
    everything = openmc.Tally(name="all energies")
    everything.filters = [openmc.MaterialFilter(order)]
    everything.scores = SCORES
    model = openmc.Model(
        geometry,
        materials,
        settings,
        openmc.Tallies([rates, flux, scattering, everything]),
    )
    statepoint = model.run(cwd=work, output=False)

    with openmc.StatePoint(statepoint) as sp:
        k = sp.keff

        def get(name: str, score: str, value: str = "mean") -> np.ndarray:
            v = sp.get_tally(name=name).get_values(scores=[score], value=value)
            return v.reshape(len(order), -1)

        absorption = get("rates", "absorption")
        variance = get("rates", "absorption", "std_dev") ** 2
        for score, extra in MULTIPLYING.items():
            absorption = absorption - extra * get("rates", score)
            variance = variance + (extra * get("rates", score, "std_dev")) ** 2
        production = get("rates", "nu-fission")
        production_std = get("rates", "nu-fission", "std_dev")
        phi = get("flux", "flux")
        phi_std = get("flux", "flux", "std_dev")
        by_scatter = get("scattering", "total") - get("scattering", "nu-scatter")
        by_scatter_std = np.sqrt(
            get("scattering", "total", "std_dev") ** 2
            + get("scattering", "nu-scatter", "std_dev") ** 2
        )

        def whole(score: str) -> float:
            return float(
                sp.get_tally(name="all energies").get_values(scores=[score]).sum()
            )

        all_absorbed = whole("absorption") - sum(
            extra * whole(score) for score, extra in MULTIPLYING.items()
        )
        # すべてのエネルギーの生成で規格化する（群の範囲の外の分は `outside` に書く）
        total = whole("nu-fission")
        out = {"format": "xlc-balance", "k": k.nominal_value, "materials": {}}
        for i, m in enumerate(order):
            # 群は高いエネルギーから
            def put(v: np.ndarray, i: int = i) -> list[float]:
                return (v[i][::-1] / total).tolist()

            out["materials"][NAMES[m.name]] = {
                "flux": put(phi),
                "flux_std": put(phi_std),
                "absorption": put(absorption),
                "absorption_std": put(np.sqrt(variance)),
                "production": put(production),
                "production_std": put(production_std),
                "absorption_by_scatter": put(by_scatter),
                "absorption_by_scatter_std": put(by_scatter_std),
            }
        absorbed = absorption.sum() / total
        absorbed_std = np.sqrt(variance.sum()) / total
        by_scatter_sum = by_scatter.sum() / total
        out["outside"] = {
            "absorption": (all_absorbed - absorption.sum()) / total,
            "production": (total - production.sum()) / total,
        }
        (work / "balance.msgpack").write_bytes(msgpack.packb(out))
        print(
            f"all energies: 1/absorption {total / all_absorbed:.6f}, "
            f"production outside the groups {out['outside']['production']:.2e}"
        )
        print(
            f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}, inside the groups: "
            f"absorption {absorbed:.6f} +- {absorbed_std:.6f}, "
            f"total - nu-scatter {by_scatter_sum:.6f}"
        )


if __name__ == "__main__":
    main()
