"""NJOY PURR の出力から、エネルギーごとの直接サンプリングの Bondarenko の表を、平均に対する比で書く。

使い方: python3 purr_table.py <PURR の output>
出力: energy sigma0 capture_ratio total_ratio（タブ区切り。比は σ0 = 1e10 の行に対する）
"""

import sys

with open(sys.argv[1]) as f:
    lines = f.read().splitlines()
energy = None
i = 0
print("energy\tsigma0\tcapture_ratio\ttotal_ratio")
while i < len(lines):
    line = lines[i]
    if line.startswith("e= "):
        energy = float(line.split()[1])
    elif "bondarenko cross sections by direct sampling" in line and energy is not None:
        rows = []
        for row in lines[i + 2 : i + 7]:
            v = [float(x) for x in row.split()]
            rows.append(v)
        mean = rows[0]
        for v in rows:
            print(f"{energy:e}\t{v[1]:e}\t{v[5] / mean[5]:.6f}\t{v[2] / mean[2]:.6f}")
        energy = None
    i += 1
