"""#0133 の試行: 同じ組成の 1C のピンセルの燃料で、核種ごとの反応率を生成 1 あたり（中性子 1 個あたり）で xlc と OpenMC で比べ、
エネルギーの範囲（熱 0.625 eV 以下、共鳴 0.625〜9118 eV、高速 9118 eV より上）に分ける。

使い方: python3 compare_rates.py <作業ディレクトリ>（openmc-rates.msgpack、xlc-fuel-xs.msgpack、xlc-balance.msgpack）[核種 ...]

- 反応: 核分裂と、非弾性 − 核分裂（全 − 弾性 − 核分裂。熱と共鳴では捕獲、高速では非弾性散乱と (n,2n) を含む）。xlc は fuel_xs の
  燃料のピンの実効断面積 × balance の燃料の中性子束（生成 1 あたり、体積で積分）× 数密度、OpenMC は燃料の群ごとの反応率 /
  全物質の生成（nu-fission）
- 中性子束の比（xlc/OpenMC − 1）は生成 1 あたりの燃料の中性子束の範囲ごとの和。反応率の差は、断面積の差と中性子束の差に分かれる
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

work = Path(sys.argv[1])
names = sys.argv[2:] or [
    "U-235",
    "U-238",
    "Pu-239",
    "Pu-240",
    "Pu-241",
    "Xe-135",
    "Sm-149",
]
omc = msgpack.unpackb((work / "openmc-rates.msgpack").read_bytes())
xs = msgpack.unpackb((work / "xlc-fuel-xs.msgpack").read_bytes())["nuclides"]
bal = msgpack.unpackb((work / "xlc-balance.msgpack").read_bytes())
bounds = np.array(omc["bounds"])
P_o = sum(np.sum(m["nu-fission"]) for m in omc["materials"].values())
phi_x = np.array(bal["materials"]["fuel"]["flux"])
phi_o = np.array(omc["flux"]) / P_o
low = bounds[1:]
ranges = [
    ("熱", low < 0.6249),
    ("共鳴", (low >= 0.6249) & (low < 9117)),
    ("高速", low >= 9117),
]
print(f"k xlc {bal['k']:.5f}, OpenMC {omc['k']:.5f} ± {omc['k_std']:.5f}")
print(
    "中性子束（生成 1 あたり）| "
    + " | ".join(f"{r}: xlc/OpenMC − 1 (%)" for r, _ in ranges)
)
print(
    "  | "
    + " | ".join(
        f"{(phi_x[m].sum() / phi_o[m].sum() - 1) * 100:+.3f}" for _, m in ranges
    )
)
print(
    "核種 | 反応 | 範囲 | OpenMC（生成 1 あたり ×1e5）± | xlc | 差 ×1e5 | xlc/OpenMC − 1 (%) | うち断面積 (%) | うち中性子束 (%)"
)
for name in names:
    key = name.replace("-", "")
    r = omc["rates"].get(name)
    if r is None or key not in xs:
        raise SystemExit(f"{name}: not in the OpenMC rates or the xlc cross sections")
    n = omc["densities"][name]
    x = xs[key]
    tot, el, fi = (np.array(x[s]) for s in ("total", "elastic", "fission"))
    for label, sx, so, so_std in [
        (
            "核分裂",
            fi,
            np.array(r["fission"]["mean"]),
            np.array(r["fission"]["std"]),
        ),
        (
            "非弾性−核分裂",
            tot - el - fi,
            # 9118 eV より下は (n,gamma)（非弾性散乱と (n,2n) のしきいより下）、上は全 − 弾性 − 核分裂（統計誤差は出さない）
            np.where(
                low < 9117,
                np.array(r["(n,gamma)"]["mean"]),
                np.array(r["total"]["mean"])
                - np.array(r["(n,elastic)"]["mean"])
                - np.array(r["fission"]["mean"]),
            ),
            np.where(low < 9117, np.array(r["(n,gamma)"]["std"]), np.nan),
        ),
    ]:
        if so.sum() == 0.0:
            continue
        for rl, m in ranges + [("計", np.ones_like(low, dtype=bool))]:
            ro = so[m].sum() / P_o
            if ro <= 0:
                continue
            rx = n * (sx[m] * phi_x[m]).sum()
            # 断面積だけの差: xlc の断面積と OpenMC の中性子束
            mix = n * (sx[m] * phi_o[m]).sum()
            std = np.sqrt((so_std[m] ** 2).sum()) / P_o
            print(
                f"{name} | {label} | {rl} | {ro * 1e5:.2f} ± {std * 1e5:.2f} | {rx * 1e5:.2f} | {(rx - ro) * 1e5:+.2f} | "
                f"{(rx / ro - 1) * 100:+.2f} | {(mix / ro - 1) * 100:+.2f} | {(rx / mix - 1) * 100:+.2f}"
            )
