"""#0133: 燃焼で増える k の差の候補ごとに、既に測った反応率の相対差と OpenMC の tally の割合から k への寄与の桁を見積もる。

使い方: python3 estimate_candidates.py <openmc-rates.msgpack>...（openmc_pin_rates*.py か openmc_pin_transfer.py の出力）

k∞ = P / A（P は生成、A は吸収。外周が鏡面反射）なので、核種と群の反応率が相対的に r だけ違えば、δk ≈ k (r ΔP/P − r ΔA/A)。
P と A は燃料の全核種の和（被覆管と減速材の吸収を A に足す）。相対差は #0133、#0148、#0151、#0152 の測定（新しい燃料か 2A の
40 MWd/kg の組成）で、状態によらず同じとみる（見積もりの仮定。確かめていない）。範囲のある差は両端を出す。
"""

import sys
from pathlib import Path

import msgpack

# (名前, [(核種, 最初の群, 最後の群, 反応 ("capture" / "fission"), 下の相対差, 上の相対差)])。群は 1 始まりで高いエネルギーから
CANDIDATES = [
    # U-238 の核分裂（15 keV より上）が生成 1 あたり +2.08%（新しい燃料。5 MWd/kg は +2.02%）。断面積 +0.83%、中性子束 +1.24%（#0148）
    ("#0148 U-238 の高速の核分裂", [("U-238", 1, 13, "fission", 0.0202, 0.0208)]),
    # 0.005〜0.1 eV の U-238 の捕獲 +0.10〜+0.58%、0.14〜0.25 eV の Pu-239 +0.39〜+0.83%、0.32〜0.5 eV の Pu-239 −0.60〜−0.09%（#0151）
    (
        "#0151 熱群の重み",
        [
            ("U-238", 57, 68, "capture", 0.0010, 0.0058),
            ("Pu-239", 50, 52, "capture", 0.0039, 0.0083),
            ("Pu-239", 50, 52, "fission", 0.0039, 0.0083),
            ("Pu-239", 46, 48, "capture", -0.0060, -0.0009),
            ("Pu-239", 46, 48, "fission", -0.0060, -0.0009),
        ],
    ),
    # 4〜9.877 eV（27 群）の燃料の中性子束の傾きで、薄い核種の捕獲が U-236 −1.7%、Nd-145 −1.5%、Ag-109 −1.3%、Tc-99 −1.2%、
    # U-238 は +0.8%（#0152。27 群の群平均の比）
    (
        "#0152 4〜10 eV の中性子束の傾き",
        [
            ("U-236", 27, 27, "capture", -0.017, -0.017),
            ("Nd-145", 27, 27, "capture", -0.015, -0.015),
            ("Ag-109", 27, 27, "capture", -0.013, -0.013),
            ("Tc-99", 27, 27, "capture", -0.012, -0.012),
            ("U-238", 27, 27, "capture", 0.008, 0.008),
        ],
    ),
    # 収率の重みを ADR-0183 にする前の数密度の差（Rh-103 −1.6%、Cs-133 −0.5%、I-135 −0.9%、Nd-149 と Pm-149 −0.5%）。今は
    # Rh-103 +0.2%、ほかは ±0.1%（#0133 の経緯）。捕獲の率は数密度に比例する
    (
        "#0133 核分裂生成物の数密度（ADR-0183 の前）",
        [
            ("Rh-103", 1, 69, "capture", -0.016, -0.016),
            ("Cs-133", 1, 69, "capture", -0.005, -0.005),
            ("I-135", 1, 69, "capture", -0.009, -0.009),
            ("Pm-149", 1, 69, "capture", -0.005, -0.005),
        ],
    ),
    ("#0133 核分裂生成物の数密度（今）", [("Rh-103", 1, 69, "capture", 0.002, 0.002)]),
]


def estimate(path: Path) -> list[tuple[str, float, float]]:
    d = msgpack.unpackb(path.read_bytes())
    rates = d["rates"]
    production = sum(sum(v["nu-fission"]["mean"]) for v in rates.values())
    absorption = sum(sum(m["absorption"]) for m in d["materials"].values())
    k = d["k"]
    out = []
    for name, parts in CANDIDATES:
        ends = []
        for end in (4, 5):
            dk = 0.0
            for part in parts:
                nuclide, first, last = part[:3]
                if nuclide not in rates:
                    continue
                r = rates[nuclide]
                g = slice(first - 1, last)
                if part[3] == "capture":
                    dk -= k * part[end] * sum(r["(n,gamma)"]["mean"][g]) / absorption
                else:
                    dk += k * part[end] * sum(r["nu-fission"]["mean"][g]) / production
                    dk -= k * part[end] * sum(r["fission"]["mean"][g]) / absorption
            ends.append(dk * 1e5)
        out.append((name, min(ends), max(ends)))
    return out


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    for path in argv:
        print(path)
        for name, low, high in estimate(Path(path)):
            span = (
                f"{low:+.1f}" if abs(high - low) < 0.05 else f"{low:+.1f}〜{high:+.1f}"
            )
            print(f"  {name}: {span} pcm")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
