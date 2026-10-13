"""xlc の CRAM と OpenMC の CRAM（openmc.deplete.cram.CRAM16、CRAM48）を、同じ燃焼行列で比べる（ADR-0161）。

使い方（references の devShell で。xlc の側は先にビルドした crates/depletion の例を呼ぶ）:

    cargo build --release -p xlc-depletion --example deplete_text
    nix develop .#references -c python3 examples/depletion/compare_openmc.py \
        examples/depletion/cases.toml target/release/examples/deplete_text

差は、誤差の絶対値の最大を初期の数密度の最大で割ったもの（scaled）で判定し、閾値を超えたら終了コード 1。
"""

import math
import random
import subprocess
import sys

import numpy as np
import scipy.sparse as sp
import tomllib
from openmc.deplete.cram import CRAM16, CRAM48

# ADR-0161 の「確認」: 同じ有理近似を同じ係数で当てるので、差は線形の解法の丸めだけ
SCALED_TOL = 1e-13


def build(case):
    """case から (核種の名前, 行列の項 [(行, 列, 値)], 初期値) を作る。"""
    if "random" in case:
        return build_random(**case["random"])
    names = case["nuclides"]
    index = {name: i for i, name in enumerate(names)}
    entries = []

    def loss_and_gain(parent, products, rate):
        p = index[parent]
        entries.append((p, p, -rate))
        for product, fraction in products:
            entries.append((index[product], p, fraction * rate))

    for parent, product, half_life in case.get("decay", []):
        loss_and_gain(parent, [(product, 1.0)], math.log(2.0) / half_life)
    flux = case.get("flux", 0.0)
    for parent, product, sigma in case.get("reaction", []):
        loss_and_gain(parent, [(product, 1.0)], sigma * 1e-24 * flux)
    for parent, sigma, yields in case.get("fission", []):
        products = list(yields.items())
        products.append(("other", 2.0 - sum(yields.values())))
        loss_and_gain(parent, products, sigma * 1e-24 * flux)
    n0 = np.zeros(len(names))
    for name, value in case["n0"].items():
        n0[index[name]] = value
    return names, entries, n0


def build_random(seed, n, fissile, yields, initial):
    rng = random.Random(seed)
    entries = []
    for col in range(n):
        others = [i for i in range(n) if i != col]
        count = yields if col < fissile else rng.randint(1, 3)
        daughters = rng.sample(others, count)
        rate = 10.0 ** rng.uniform(-12.0, 2.0)
        entries.append((col, col, -rate))
        for d in daughters:
            entries.append((d, col, rate / len(daughters)))
    n0 = np.zeros(n)
    for i in rng.sample(range(n), initial):
        n0[i] = 10.0 ** rng.uniform(-6.0, 0.0)
    return [f"x{i}" for i in range(n)], entries, n0


def run_xlc(binary, n, entries, n0, t, order):
    lines = [
        f"n {n}",
        f"order {order}",
        f"t {t!r}",
        "n0 " + " ".join(repr(float(v)) for v in n0),
    ]
    lines += [f"entry {r} {c} {v!r}" for r, c, v in entries]
    out = subprocess.run(
        [binary], input="\n".join(lines), capture_output=True, text=True, check=True
    )
    return np.array([float(v) for v in out.stdout.split()])


def main():
    cases_path, binary = sys.argv[1], sys.argv[2]
    with open(cases_path, "rb") as f:
        cases = tomllib.load(f)["case"]
    failed = False
    print(
        f"{'case':12} {'order':>5} {'scaled xlc-openmc':>18} {'max rel (n>1e-12 max n0)':>25} {'scaled 16-48 (openmc)':>22}"
    )
    for case in cases:
        names, entries, n0 = build(case)
        n = len(names)
        rows, cols, vals = zip(*entries)
        # coo_matrix の tocsr は同じ位置の値を足す（xlc の from_triplets と同じ）
        a = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
        t = case["t"]
        scale = np.max(np.abs(n0))
        ref = {16: CRAM16(a, n0, t), 48: CRAM48(a, n0, t)}
        order_gap = np.max(np.abs(ref[16] - ref[48])) / scale
        for order in (16, 48):
            got = run_xlc(binary, n, entries, n0, t, order)
            diff = np.abs(got - ref[order])
            scaled = np.max(diff) / scale
            big = np.abs(ref[order]) > 1e-12 * scale
            rel = np.max(diff[big] / np.abs(ref[order][big]))
            ok = bool(
                np.isfinite(got).all()
                and np.isfinite(ref[order]).all()
                and scaled <= SCALED_TOL
            )
            mark = "" if ok else "  FAIL"
            failed |= not ok
            print(
                f"{case['name']:12} {order:>5} {scaled:>18.3e} {rel:>25.3e} {order_gap:>22.3e}{mark}"
            )
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
