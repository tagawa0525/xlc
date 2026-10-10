"""均質な無限媒質（`homogeneous_medium.py` の入力）を OpenMC の連続エネルギーモードで解き、k と、核種ごと・群ごとの
中性子の釣り合いを書く。

使い方: python3 openmc_homogeneous.py <cross_sections.xml> <媒質 (TOML)> <作業ディレクトリ> <xlc のライブラリ>
        [粒子数/世代] [世代数] [捨てる世代数] [乱数の種]

形状は外周が鏡面反射の立方体 1 つ。H-1 があれば、軽水の H の熱中性子散乱（c_H_in_H2O）を使う（xlc の H-1 と
同じ扱い）。釣り合いは `openmc_balance.py` と同じ形式の balance.msgpack で、物質の代わりに核種（xlc の名前）ごとに、
群ごとの中性子束（どの核種も同じ）、正味の吸収、生成、散乱による移行（nu-scatter）を書く。値はすべてのエネルギーの
生成の和が 1 になるように世代ごとに規格化し、世代ごとの値を `batches` に書く。群の範囲の外の正味の吸収と生成は
`outside` に書く。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import openmc
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from openmc_balance import MULTIPLYING, net_absorption
from openmc_common import eigenvalue_settings, group_filter, material

SCORES = ["absorption", "nu-fission", *MULTIPLYING]


def main() -> None:
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    cross_sections = Path(sys.argv[1]).resolve()
    medium = tomllib.loads(Path(sys.argv[2]).read_text())
    work = Path(sys.argv[3]).resolve()
    library = Path(sys.argv[4]).resolve()
    particles = int(sys.argv[5]) if len(sys.argv) > 5 else 200_000
    batches = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    inactive = int(sys.argv[7]) if len(sys.argv) > 7 else 50
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    work.mkdir(parents=True, exist_ok=True)

    names = list(medium["nuclides"])
    densities = {k.replace("-", ""): v for k, v in medium["nuclides"].items()}
    nuclides = list(densities)
    mix = material("mix", densities, medium["temperature"])
    if "H1" in densities:
        mix.add_s_alpha_beta("c_H_in_H2O")
    materials = openmc.Materials([mix])
    materials.cross_sections = str(cross_sections)
    half = 1.0
    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    top = openmc.ZPlane(z0=half, boundary_type="reflective")
    bottom = openmc.ZPlane(z0=-half, boundary_type="reflective")
    geometry = openmc.Geometry([openmc.Cell(fill=mix, region=-box & -top & +bottom)])
    settings = eigenvalue_settings(particles, batches, inactive, seed, half)
    settings.temperature = {"method": "nearest", "tolerance": 1.0}
    # 世代ごとの値を、続く 2 つの statepoint の和の差から求める
    settings.statepoint = {"batches": list(range(inactive + 1, batches + 1))}

    groups = group_filter(library)
    low, high = groups.bins[0][0], groups.bins[-1][1]
    flux = openmc.Tally(name="flux")
    flux.filters = [groups]
    flux.scores = ["flux"]
    rates = openmc.Tally(name="rates")
    rates.filters = [groups]
    rates.nuclides = nuclides
    rates.scores = SCORES
    # 群の範囲の下、範囲の中、範囲の上
    outside = openmc.Tally(name="outside")
    outside.filters = [openmc.EnergyFilter([0.0, low, high, 1.0e10])]
    outside.nuclides = nuclides
    outside.scores = SCORES
    # 出来事ごとの推定が要る score は別の tally にする（同じ tally だと全体が出来事ごとの推定になる）
    transfer = openmc.Tally(name="transfer")
    transfer.filters = [groups, openmc.EnergyoutFilter(groups.values)]
    transfer.nuclides = nuclides
    transfer.scores = ["nu-scatter"]
    model = openmc.Model(
        geometry,
        materials,
        settings,
        openmc.Tallies([flux, rates, outside, transfer]),
    )
    model.run(cwd=work, output=False)

    n = len(groups.bins)
    m = len(nuclides)
    width = len(str(batches))
    previous = None
    per_batch = []
    for batch in range(inactive + 1, batches + 1):
        with openmc.StatePoint(work / f"statepoint.{batch:0{width}d}.h5") as sp:
            now = {
                "flux": sp.get_tally(name="flux").get_values(value="sum").reshape(n),
                "transfer": sp.get_tally(name="transfer")
                .get_values(value="sum")
                .reshape(n, n, m),
            }
            for name in ("rates", "outside"):
                t = sp.get_tally(name=name)
                for score in SCORES:
                    now[(name, score)] = t.get_values(
                        scores=[score], value="sum"
                    ).reshape(-1, m)
            if batch == batches:
                k = sp.keff
        this = now if previous is None else {s: now[s] - previous[s] for s in now}
        previous = now
        inside = {s: this[("rates", s)].T for s in SCORES}
        edges = {s: this[("outside", s)].T for s in SCORES}
        total = edges["nu-fission"].sum()
        net_edges = net_absorption(edges)
        per_batch.append(
            {
                "flux": this["flux"] / total,
                "absorption": net_absorption(inside) / total,
                "production": inside["nu-fission"] / total,
                "transfer": np.moveaxis(this["transfer"], 2, 0) / total,
                "below": net_edges[:, 0].sum() / total,
                "above": net_edges[:, 2].sum() / total,
                "below_production": edges["nu-fission"][:, 0].sum() / total,
                "above_production": edges["nu-fission"][:, 2].sum() / total,
            }
        )

    def rates_of(b: dict, i: int) -> dict:
        """1 組の値の、核種 i の値（群は高いエネルギーから）"""
        moved = b["transfer"][i]
        within = np.diagonal(moved)
        return {
            "flux": b["flux"][::-1].tolist(),
            "absorption": b["absorption"][i][::-1].tolist(),
            "production": b["production"][i][::-1].tolist(),
            "out_scatter": (moved.sum(axis=1) - within)[::-1].tolist(),
            "within_scatter": within[::-1].tolist(),
        }

    means = {key: np.mean([b[key] for b in per_batch], axis=0) for key in per_batch[0]}
    out = {
        "format": "xlc-balance",
        "k": k.nominal_value,
        "k_std": k.std_dev,
        "materials": {},
        "outside": {
            side: {
                "absorption": float(means[side]),
                "production": float(means[f"{side}_production"]),
            }
            for side in ("below", "above")
        },
        "batches": [
            {
                "materials": {name: rates_of(b, i) for i, name in enumerate(names)},
                "below": float(b["below"]),
                "above": float(b["above"]),
            }
            for b in per_batch
        ],
    }
    # 移行の行列は [入射][出射]
    for i, name in enumerate(names):
        out["materials"][name] = rates_of(means, i)
        out["materials"][name]["transfer"] = means["transfer"][i][::-1, ::-1].tolist()
    (work / "balance.msgpack").write_bytes(msgpack.packb(out))
    absorbed = np.array(
        [b["absorption"].sum() + b["below"] + b["above"] for b in per_batch]
    )
    inverse = 1 / absorbed.mean()
    inverse_std = absorbed.std(ddof=1) / np.sqrt(len(absorbed)) * inverse**2
    print(
        f"{sys.argv[2]} k {k.nominal_value:.6f} +- {k.std_dev:.6f}, "
        f"1/(net absorption over all energies) {inverse:.6f} +- {inverse_std:.6f}"
    )


if __name__ == "__main__":
    main()
