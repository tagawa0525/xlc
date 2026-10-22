"""#0133 の試行: 同じ組成の 1C のピンセルで、xlc と OpenMC の生成 1 あたりの吸収（1/k の内訳）の差を、物質とエネルギーの範囲、
燃料の核種に分ける。

使い方: python3 decompose.py <xlc の作業ディレクトリ（xlc-balance.msgpack、xlc-fuel-xs.msgpack）> <OpenMC の openmc-rates.msgpack>

- 物質: xlc は `balance` の正味の吸収（全 − 散乱の行の和）、OpenMC は absorption（(n,2n) などの増倍を引かない。高速の群の差に含まれる）
- 燃料の核種: 9118 eV より下（非弾性散乱のしきいより下とみなす）で、xlc は (全 − 弾性) × 中性子束、OpenMC は absorption
差は (xlc − OpenMC) × 1e5（pcm 相当。Δ(1/k) ≈ −Δk/k²）
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

work = Path(sys.argv[1])
bal = msgpack.unpackb((work / "xlc-balance.msgpack").read_bytes())
xs = msgpack.unpackb((work / "xlc-fuel-xs.msgpack").read_bytes())["nuclides"]
omc = msgpack.unpackb(Path(sys.argv[2]).read_bytes())
bounds = np.array(omc["bounds"])
names = {"fuel": "fuel-0", "zircaloy": "zircaloy-4", "water": "moderator"}
P_o = sum(np.sum(m["nu-fission"]) for m in omc["materials"].values())
k_x = bal["k"]
print(
    f"k xlc {k_x:.5f}, OpenMC {omc['k']:.5f} ± {omc['k_std']:.5f}: 差 {(k_x - omc['k']) * 1e5:+.0f} pcm, "
    f"1/k の差 ×1e5 {(1 / k_x - 1 / omc['k']) * 1e5:+.1f}"
)
ranges = [
    ("> 9118 eV", bounds[:-1] > 9118 * 0.999),
    ("4〜9118 eV", (bounds[:-1] <= 9118 * 1.001) & (bounds[:-1] > 4.0 * 1.001)),
    ("0.625〜4 eV", (bounds[:-1] <= 4.001) & (bounds[:-1] > 0.626)),
    ("< 0.625 eV", bounds[:-1] <= 0.626),
]
total = 0.0
print("物質 | 範囲 | xlc A/P | OpenMC A/P | 差 ×1e5")
for x, o in names.items():
    ax = np.array(bal["materials"][x]["absorption"])
    ao = np.array(omc["materials"][o]["absorption"]) / P_o
    for label, mask in ranges:
        d = (ax[mask].sum() - ao[mask].sum()) * 1e5
        total += d
        print(f"{x} | {label} | {ax[mask].sum():.5f} | {ao[mask].sum():.5f} | {d:+.1f}")
    px = np.sum(bal["materials"][x]["production"])
    po = np.sum(omc["materials"][o]["nu-fission"]) / P_o
print(f"吸収の差の和 ×1e5 {total:+.1f}")
A_o = sum(np.sum(m["absorption"]) for m in omc["materials"].values()) / P_o
print(
    f"OpenMC の Σ A/P − 1/k ×1e5 {(A_o - 1 / omc['k']) * 1e5:+.1f}（(n,2n) などの増倍。高速の群の燃料の差に入る）"
)
phi = np.array(bal["materials"]["fuel"]["flux"])
low = bounds[:-1] <= 9118 * 1.001
rows = []
for name, r in omc["rates"].items():
    key = name.replace("-", "")
    n = omc["densities"][name]
    sx = np.array(xs[key]["total"]) - np.array(xs[key]["elastic"])
    ax = n * sx * phi
    ao = np.array(r["absorption"]["mean"]) / P_o
    rows.append(
        (name, (ax[low].sum() - ao[low].sum()) * 1e5, ao[low].sum(), ax[low].sum())
    )
rows.sort(key=lambda t: -abs(t[1]))
print("燃料の核種（9118 eV より下）| 差 ×1e5 | OpenMC A/P | xlc/OpenMC − 1 (%)")
for name, d, ao, ax in rows[:20]:
    print(f"{name} | {d:+.1f} | {ao:.5f} | {(ax / ao - 1) * 100:+.2f}")
print(f"和（全核種） {sum(t[1] for t in rows):+.1f}")
