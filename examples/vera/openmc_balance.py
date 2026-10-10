"""VERA の 1B〜1D を OpenMC の連続エネルギーモードで解き、物質ごと・群ごとの中性子束、正味の吸収、生成と、群の間の散乱に
よる移行を書く。

使い方: python3 openmc_balance.py <cross_sections.xml> <問題> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

xlc の `balance` の例（crates/lattice/examples/balance.rs）と同じ形式の balance.msgpack を書く。群構造は xlc の
ライブラリのもので、高いエネルギーから並べる。正味の吸収は、xlc（全断面積から、中性子を出すすべての反応の放出の数を
含む散乱行列の行の和を引く）と同じ量にするため、吸収（中性子が出ない反応と核分裂）から、2 個以上の中性子を出す反応の
増える分（放出の数 − 1）を引く。値はすべてのエネルギーの生成の和が 1 になるように規格化する（外周が鏡面反射なので、
正味の吸収の和が 1/k になる）。群の範囲の外（いちばん下の境より下と、いちばん上の境より上）の正味の吸収と生成は
`outside` の `below` と `above` に書く。

統計誤差は、生成、吸収、中性子束が同じ粒子の追跡から出て強く相関するので、世代ごとの値から求める。世代ごとの値
（それぞれの世代のすべてのエネルギーの生成で規格化したもの）を `batches` に書く。移行の行列（nu-scatter、出来事ごとの
推定）は平均と標準偏差を書き、世代ごとには、入射の群ごとに、ほかの群へ出る分（`out_scatter`）と同じ群に残る分
（`within_scatter`）だけを書く。
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
SCORES = ["absorption", "nu-fission", "flux", *MULTIPLYING]
KEYS = ("flux", "absorption", "production")


def net_absorption(values: dict[str, np.ndarray]) -> np.ndarray:
    """score ごとの値から正味の吸収を求める。"""
    net = values["absorption"]
    for score, extra in MULTIPLYING.items():
        net = net - extra * values[score]
    return net


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
    # 世代ごとの値を、続く 2 つの statepoint の和の差から求める
    settings.statepoint = {"batches": list(range(inactive + 1, batches + 1))}

    order = [fuel, clad, water]
    groups = group_filter(library)
    low, high = groups.bins[0][0], groups.bins[-1][1]
    rates = openmc.Tally(name="rates")
    rates.filters = [openmc.MaterialFilter(order), groups]
    rates.scores = SCORES
    # 群の範囲の下、範囲の中、範囲の上
    outside = openmc.Tally(name="outside")
    outside.filters = [
        openmc.MaterialFilter(order),
        openmc.EnergyFilter([0.0, low, high, 1.0e10]),
    ]
    outside.scores = SCORES
    # 出来事ごとの推定が要る score は別の tally にする（同じ tally だと全体が出来事ごとの推定になる）
    transfer = openmc.Tally(name="transfer")
    transfer.filters = [
        openmc.MaterialFilter(order),
        groups,
        openmc.EnergyoutFilter(groups.values),
    ]
    transfer.scores = ["nu-scatter"]
    model = openmc.Model(
        geometry, materials, settings, openmc.Tallies([rates, outside, transfer])
    )
    model.run(cwd=work, output=False)

    n = len(groups.bins)
    width = len(str(batches))
    previous = None
    per_batch = []
    raw_totals = []
    for batch in range(inactive + 1, batches + 1):
        with openmc.StatePoint(work / f"statepoint.{batch:0{width}d}.h5") as sp:
            now = {
                name: {
                    score: sp.get_tally(name=name).get_values(
                        scores=[score], value="sum"
                    )
                    for score in SCORES
                }
                for name in ("rates", "outside")
            }
            now["transfer"] = {
                "nu-scatter": sp.get_tally(name="transfer").get_values(value="sum")
            }
            if batch == batches:
                k = sp.keff
                tally = sp.get_tally(name="transfer")
                mean = tally.get_values(value="mean").reshape(len(order), n, n)
                std = tally.get_values(value="std_dev").reshape(len(order), n, n)
        this = (
            now
            if previous is None
            else {
                name: {s: now[name][s] - previous[name][s] for s in now[name]}
                for name in now
            }
        )
        previous = now
        inside = {s: v.reshape(len(order), n) for s, v in this["rates"].items()}
        edges = {s: v.reshape(len(order), 3) for s, v in this["outside"].items()}
        total = edges["nu-fission"].sum()
        raw_totals.append(total)
        net_edges = net_absorption(edges)
        moved = this["transfer"]["nu-scatter"].reshape(len(order), n, n)
        within = np.diagonal(moved, axis1=1, axis2=2)
        per_batch.append(
            {
                "flux": inside["flux"] / total,
                "absorption": net_absorption(inside) / total,
                "production": inside["nu-fission"] / total,
                "below": net_edges[:, 0].sum() / total,
                "above": net_edges[:, 2].sum() / total,
                "below_production": edges["nu-fission"][:, 0].sum() / total,
                "above_production": edges["nu-fission"][:, 2].sum() / total,
                "out_scatter": (moved.sum(axis=2) - within) / total,
                "within_scatter": within / total,
            }
        )

    def mean_of(key: str) -> np.ndarray:
        return np.mean([b[key] for b in per_batch], axis=0)

    # 移行の行列は 1 世代あたりの平均なので、世代の生成の平均で規格化する
    per_source = np.mean(raw_totals)
    out = {
        "format": "xlc-balance",
        "k": k.nominal_value,
        "k_std": k.std_dev,
        "materials": {},
        "outside": {
            side: {
                "absorption": float(mean_of(side)),
                "production": float(mean_of(f"{side}_production")),
            }
            for side in ("below", "above")
        },
        "batches": [],
    }
    # 群は高いエネルギーから。行列は [入射][出射]
    for i, m in enumerate(order):
        name = NAMES[m.name]
        out["materials"][name] = {key: mean_of(key)[i][::-1].tolist() for key in KEYS}
        out["materials"][name]["transfer"] = (mean[i][::-1, ::-1] / per_source).tolist()
        out["materials"][name]["transfer_std"] = (
            std[i][::-1, ::-1] / per_source
        ).tolist()
    for b in per_batch:
        out["batches"].append(
            {
                "materials": {
                    NAMES[m.name]: {
                        key: b[key][i][::-1].tolist()
                        for key in (*KEYS, "out_scatter", "within_scatter")
                    }
                    for i, m in enumerate(order)
                },
                "below": float(b["below"]),
                "above": float(b["above"]),
            }
        )
    (work / "balance.msgpack").write_bytes(msgpack.packb(out))
    absorbed = np.array(
        [b["absorption"].sum() + b["below"] + b["above"] for b in per_batch]
    )
    inverse = 1 / absorbed.mean()
    inverse_std = absorbed.std(ddof=1) / np.sqrt(len(absorbed)) * inverse**2
    print(
        f"{problem} k {k.nominal_value:.6f} +- {k.std_dev:.6f}, "
        f"1/(net absorption over all energies) {inverse:.6f} +- {inverse_std:.6f}, "
        f"production below/above the groups {out['outside']['below']['production']:.2e}/"
        f"{out['outside']['above']['production']:.2e}"
    )


if __name__ == "__main__":
    main()
