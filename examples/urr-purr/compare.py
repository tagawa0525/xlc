"""xlc の遮蔽の比（urr_convergence の出力）を NJOY PURR の直接サンプリングの表（purr_table.py の出力）と並べる。

使い方: python3 compare.py <xlc の tsv> <PURR の tsv>

ADR-0177 の「確認」の表の作り方（565 K、ENDF/B-VIII.0。<核種> は U238、Pu240、Pm149）:

- PURR: <核種>/input を、ENDF を tape20 にした作業ディレクトリで references の devShell の njoy に与え、出力を
  `python3 purr_table.py output > purr-<核種>.tsv`
- xlc: `target/release/examples/urr_convergence <ENDF> 565 512 250 <エネルギー>... > xlc-<核種>.tsv`（エネルギーは各 tsv の値）
出力: 核種のエネルギー、σ0、捕獲の比（xlc、種 1 つあたりの相対標準偏差、PURR、差）、全断面積の比の差と相対標準偏差
"""

import sys


def read(path):
    rows = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "energy")):
                continue
            v = [float(x) for x in line.split()]
            rows[(v[0], v[1])] = v[2:]
    return rows


xlc, purr = read(sys.argv[1]), read(sys.argv[2])
print(
    "energy\tsigma0\tcapture_xlc\trel_std\tcapture_purr\tdiff_%\ttotal_diff_%\ttotal_rel_std"
)
for key, (c, cs, t, ts) in xlc.items():
    if key not in purr or key[1] >= 1e9:
        continue
    pc, pt = purr[key]
    print(
        f"{key[0]:g}\t{key[1]:g}\t{c:.5f}\t{cs:.2e}\t{pc:.5f}\t{100 * (c / pc - 1):+.2f}\t{100 * (t / pt - 1):+.2f}\t{ts:.2e}"
    )
