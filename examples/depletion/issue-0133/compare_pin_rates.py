"""#0133 の試行: 同じ組成の 1C のピンセルの燃料で、核種の 1 群の微視的断面積を xlc と OpenMC で比べる。

使い方: python3 compare_pin_rates.py <作業ディレクトリ>（openmc-rates.msgpack、xlc-fuel-xs.msgpack、xlc-balance.msgpack）[群の詳細を出す核種 ...]

xlc の断面積は `fuel_xs` の燃料のピンの実効断面積（全 − 弾性 = 非弾性。捕獲の核種では捕獲）、中性子束は `balance` の燃料の
群の中性子束。OpenMC は燃料の群ごとの反応率 / (N φ)。1 群は
  OpenMC: Σ R_g / (N Σ φ_g)
  xlc:    Σ σ_g φ_g / Σ φ_g（xlc の中性子束）
  混合:   Σ σ_g(xlc) φ_g(OpenMC) / Σ φ_g(OpenMC)（断面積だけの差）
"""

import sys
from pathlib import Path

import msgpack
import numpy as np


def main():
    work = Path(sys.argv[1])
    detail = sys.argv[2:]
    omc = msgpack.unpackb((work / "openmc-rates.msgpack").read_bytes())
    xs = msgpack.unpackb((work / "xlc-fuel-xs.msgpack").read_bytes())["nuclides"]
    bal = msgpack.unpackb((work / "xlc-balance.msgpack").read_bytes())
    phi_x = np.array(bal["materials"]["fuel"]["flux"])
    phi_o = np.array(omc["flux"])
    bounds = omc["bounds"]
    print(f"k OpenMC {omc['k']:.5f} ± {omc['k_std']:.5f}, xlc {bal['k']:.5f}")
    # 中性子束のスペクトル（和で規格化）の比べ: 熱（0.625 eV 以下）の割合
    thermal = [g for g in range(len(bounds) - 1) if bounds[g] <= 0.625 + 1e-9]
    fx, fo = phi_x / phi_x.sum(), phi_o / phi_o.sum()
    print(
        f"熱群の中性子束の割合 xlc {fx[thermal].sum():.5f}, OpenMC {fo[thermal].sum():.5f}"
    )
    print(
        "核種 | 反応 | OpenMC 1 群 (b) ± | xlc 1 群 | xlc/OpenMC − 1 (%) | 混合/OpenMC − 1 (%) | スペクトル（xlc/混合 − 1, %）"
    )
    for name, r in omc["rates"].items():
        n = omc["densities"].get(name)
        key = name.replace("-", "")
        if n is None or key not in xs:
            continue
        x = xs[key]
        for score, ours in [("(n,gamma)", None), ("fission", "fission")]:
            R = np.array(r[score]["mean"])
            if R.sum() == 0.0:
                continue
            if ours is None:
                # 非弾性 − 核分裂（非弾性の散乱は熱と共鳴で 0。高速の群は U-238 などで含むので、ここでは捕獲の核種だけを見る）
                sx = (
                    np.array(x["total"])
                    - np.array(x["elastic"])
                    - np.array(x["fission"])
                )
            else:
                sx = np.array(x["fission"])
            so = R.sum() / (n * phi_o.sum())
            so_std = np.sqrt(np.sum(np.array(r[score]["std"]) ** 2)) / (n * phi_o.sum())
            ox = (sx * phi_x).sum() / phi_x.sum()
            mix = (sx * phi_o).sum() / phi_o.sum()
            print(
                f"{name} | {score} | {so:.5g} ± {so_std / so * 100:.2f}% | {ox:.5g} | {(ox / so - 1) * 100:+.2f} | "
                f"{(mix / so - 1) * 100:+.2f} | {(ox / mix - 1) * 100:+.2f}"
            )
            if name in detail:
                sg_o = np.divide(R, n * phi_o, out=np.zeros_like(R), where=phi_o > 0)
                share = R / R.sum()
                print(
                    "   群 | 上端 (eV) | 反応の割合 | σ OpenMC | σ xlc | 差 (%) | 1 群への寄与 (%) | φ xlc/OpenMC − 1 (%)"
                )
                for g in range(len(R)):
                    if share[g] < 0.005:
                        continue
                    contrib = (sx[g] - sg_o[g]) * phi_o[g] / phi_o.sum() / so * 100
                    print(
                        f"   {g + 1} | {bounds[g]:.4g} | {share[g]:.3f} | {sg_o[g]:.5g} | {sx[g]:.5g} | "
                        f"{(sx[g] / sg_o[g] - 1) * 100:+.2f} | {contrib:+.3f} | {(fx[g] / fo[g] - 1) * 100:+.2f}"
                    )


main()
