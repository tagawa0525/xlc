"""区間に分けた燃焼（openmc_depletion_resume.py）が、分けない燃焼（openmc_depletion.py）と統計の範囲で同じ k と数密度を出すことの確認（ADR-0174）。

使い方（split_check.sh が流して呼ぶ）:

    python3 compare_split.py <分けない走行の msgpack>... -- <区間に分けた走行の msgpack>...

2 つの列は同じ数で、同じ位置のものが同じ乱数の種の対。種ごとの相対差 (分けた - 分けない) / 分けない の平均と標準誤差 (n-1 の
標準偏差 / sqrt(n)) から、対のある t を、1 点目と 2 点目の k と、数密度が 1e-9 /(b cm) を超える核種（最後の点）について出す。
同じ種の対は 1 点目までの輸送が同じなので、独立な 2 組より差の分散が小さい。n が小さいので t は自由度 n-1 の Student の t で見る
（正規分布の閾値では、外れを多く見積もる）。
"""

import math
import statistics
import sys

import msgpack


def load(path: str) -> dict:
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read())


def last_density(result: dict, nuclide: str) -> float:
    return result["rings"][0]["densities"][nuclide][-1]


def paired(whole: list[dict], split: list[dict], get) -> tuple[float, float]:
    diffs = [(get(s) - get(w)) / get(w) for w, s in zip(whole, split, strict=True)]
    return statistics.mean(diffs), statistics.stdev(diffs) / math.sqrt(len(diffs))


def main(argv: list[str]) -> None:
    mid = argv.index("--")
    whole = [load(p) for p in argv[:mid]]
    split = [load(p) for p in argv[mid + 1 :]]
    n = len(whole)
    print(f"{n} seeds, paired by seed")
    for i in (1, 2):
        mean, se = paired(whole, split, lambda r, i=i: r["k"][i])
        print(f"k[{i}]: split/whole-1 {mean:+.2e} +- {se:.2e}  t {mean / se:+.2f}")
    nuclides = [
        nuc
        for nuc in whole[0]["rings"][0]["densities"]
        if last_density(whole[0], nuc) > 1e-9
    ]
    ts = {}
    for nuc in nuclides:
        mean, se = paired(whole, split, lambda r, nuc=nuc: last_density(r, nuc))
        ts[nuc] = (mean, se, mean / se if se else 0.0)
    for nuc in ("U235", "U238", "Np239", "Pu239", "Xe135", "Sm149"):
        mean, se, t = ts[nuc]
        print(f"{nuc}: split/whole-1 {mean:+.2e} +- {se:.2e}  t {t:+.2f}")
    over = {c: sum(abs(t) > c for _, _, t in ts.values()) for c in (2, 3)}
    print(
        f"{len(ts)} nuclides above 1e-9: |t|>2 {over[2]}, |t|>3 {over[3]}; "
        f"largest {sorted(ts, key=lambda k: -abs(ts[k][2]))[:4]}"
    )
    print(
        "expected for no difference (Student t, 7 degrees of freedom, 151 nuclides): "
        "|t|>2 12.9, |t|>3 3.0"
    )


if __name__ == "__main__":
    main(sys.argv[1:])
