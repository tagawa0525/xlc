"""#0133: 燃焼で増える k の差の候補ごとに、既に測った反応率の相対差と OpenMC の tally の割合から k への寄与の桁を見積もる。

使い方: python3 estimate_candidates.py <openmc-rates.msgpack>...（openmc_pin_rates*.py の出力。物質ごとの吸収（`materials` の
`absorption`）を読むので、openmc_pin_transfer.py の出力は受け付けない）

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
            ("Pu-239", 53, 55, "capture", 0.0039, 0.0083),
            ("Pu-239", 53, 55, "fission", 0.0039, 0.0083),
            ("Pu-239", 47, 49, "capture", -0.0060, -0.0009),
            ("Pu-239", 47, 49, "fission", -0.0060, -0.0009),
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


def estimate(
    d: dict, candidates: list = CANDIDATES
) -> list[tuple[str, float, float, list[str]]]:
    """候補ごとの (名前, 下限, 上限, 組成に無くて数えなかった核種)。寄与は測った比ごとに線形なので、和の下限と上限は比ごとの両端の
    小さい方と大きい方の和。同じ核種と群の範囲と比の範囲を持つ部分（捕獲と核分裂）は 1 つの測った比を共有するので、先に係数を
    足してから端を取る（別々に端を取ると、どの入力も許さない組み合わせになる）"""
    rates = d["rates"]
    production = sum(sum(v["nu-fission"]["mean"]) for v in rates.values())
    absorption = sum(sum(m["absorption"]) for m in d["materials"].values())
    k = d["k"]
    out = []
    for name, parts in candidates:
        coefficients: dict[tuple, float] = {}
        missing = []
        for nuclide, first, last, reaction, *ends in parts:
            if reaction not in ("capture", "fission"):
                raise ValueError(f"{name}: unknown reaction `{reaction}`")
            if nuclide not in rates:
                if nuclide not in missing:
                    missing.append(nuclide)
                continue
            r = rates[nuclide]
            g = slice(first - 1, last)
            if reaction == "capture":
                per = -k * sum(r["(n,gamma)"]["mean"][g]) / absorption
            else:
                per = k * (
                    sum(r["nu-fission"]["mean"][g]) / production
                    - sum(r["fission"]["mean"][g]) / absorption
                )
            key = (nuclide, first, last, *ends)
            coefficients[key] = coefficients.get(key, 0.0) + per
        low = high = 0.0
        for (*_, a, b), per in coefficients.items():
            values = [per * e * 1e5 for e in (a, b)]
            low += min(values)
            high += max(values)
        out.append((name, low, high, missing))
    return out


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__, file=sys.stderr)
        return 2
    for path in argv:
        print(path)
        for name, low, high, missing in estimate(
            msgpack.unpackb(Path(path).read_bytes())
        ):
            span = (
                f"{low:+.1f}" if abs(high - low) < 0.05 else f"{low:+.1f}〜{high:+.1f}"
            )
            note = f"（組成に無い: {'、'.join(missing)}）" if missing else ""
            print(f"  {name}: {span} pcm{note}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
