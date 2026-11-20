"""MF6 の中性子の LAW=6（位相空間分布）を、NJOY GROUPR の getmf6 と同じ入射エネルギーの格子（多重度の最初の点 / 1.1999 から
1.2 倍ずつ、多重度の最後の点まで）で LAW=1 の表にし、入射エネルギーの間を単位基底（INT=22）で補間するコピーを書く。
xs を NJOY の格子で解かせ、NJOY との差が格子によるかを見る。

`cm` は重心系の表（LCT=2）、`lab` は NJOY と同じく格子の点ごとに実験室系に移し、E' = 0 から始めた表（LCT=1）。実験室系の密度は閉じた式
f(E') = K/(4 c √E_max (m + 1)) [(1 − s_a)^(m+1) − (1 − s_b)^(m+1)]、c = √Ec、s = v²/E_max、v は [|√E' − c|, √E' + c] と
[0, √E_max] の重なりの両端（重心系で等方な √E' = v の中性子が実験室系で一様に分布することから）。

使い方: python3 -I tabulate_phase_space.py <ENDF> <MT> <Q eV> <cm|lab> <コピー>
"""

import math
import sys
from pathlib import Path

src, mt, q, frame, dst = (
    sys.argv[1],
    int(sys.argv[2]),
    float(sys.argv[3]),
    sys.argv[4],
    sys.argv[5],
)
assert frame in ("cm", "lab")


def num(s):
    s = s.strip()
    if not s:
        return 0.0
    for i in range(1, len(s)):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def fmt(v):
    if v == 0:
        return " 0.000000+0"
    m, e = f"{v:.6e}".split("e")
    e = int(e)
    return f"{m}{'+' if e >= 0 else '-'}{abs(e)}".rjust(11)


def rows(values, tail):
    out = []
    for k in range(0, len(values), 6):
        chunk = values[k : k + 6]
        out.append("".join(chunk).ljust(66) + tail + "\n")
    return out


lines = Path(src).read_text().splitlines(keepends=True)
idx = [i for i, l in enumerate(lines) if l[70:72] == " 6" and int(l[72:75]) == mt]
first = idx[0]
tail = lines[first][66:75]
head = lines[first]
awr = num(head[11:22])
nk = int(head[44:55])
# 中性子の多重度の TAB1（最初の生成粒子とする）
tab = lines[first + 1]
assert num(tab[0:11]) == 1.0 and int(tab[33:44]) == 6, tab
nr, np_ = int(tab[44:55]), int(tab[55:66])
nbt_lines = (2 * nr + 5) // 6
pt_lines = (2 * np_ + 5) // 6
pts = []
for l in lines[first + 2 + nbt_lines : first + 2 + nbt_lines + pt_lines]:
    pts += [num(l[k * 11 : (k + 1) * 11]) for k in range(6)]
energies_mult = pts[0 : 2 * np_ : 2]
cont = lines[first + 2 + nbt_lines + pt_lines]
apsx, npsx = num(cont[0:11]), int(cont[55:66])
end = first + 2 + nbt_lines + pt_lines + 1
m = 1.5 * npsx - 4
k_norm = {3: 8 / math.pi, 4: 105 / 16, 5: 256 / (7 * math.pi)}[npsx]

grid = []
enow = energies_mult[0] / 1.1999
pspmax = energies_mult[-1]
while enow < pspmax:
    enow = min(enow * 1.2, pspmax)
    grid.append(enow)

out = []
out.append(head[:33] + f"{2 if frame == 'cm' else 1:>11}" + head[44:])
out.append(tab[:33] + f"{1:>11}" + tab[44:])  # LAW=1
out += lines[first + 2 : first + 2 + nbt_lines + pt_lines]
out += rows(
    [fmt(0), fmt(0), f"{1:>11}", f"{2:>11}", f"{1:>11}", f"{len(grid):>11}"], tail
)
out += rows([f"{len(grid):>11}", f"{22:>11}"], tail)
n = 400
for e in grid:
    emax = (apsx - 1) / apsx * (awr / (awr + 1) * e + q)
    if emax <= 0:
        emax = 1.0
    xs, ys = [], []
    if frame == "cm":
        for i in range(n + 1):
            t = 0.5 * math.pi * i / n
            u = math.sin(t) ** 2
            xs.append(u * emax)
            ys.append(k_norm * math.sqrt(u) * max(1 - u, 0.0) ** m / emax)
    else:
        c, top = math.sqrt(e / (awr + 1) ** 2), math.sqrt(emax)
        lo, hi = max(c - top, 0.0), c + top
        roots = sorted(
            {lo + (hi - lo) * i / n for i in range(n + 1)}
            | {r for r in (c, top - c) if lo < r < hi}
        )

        def density(r, c=c, top=top, emax=emax):
            a, b = abs(r - c), min(r + c, top)
            if a >= b:
                return 0.0
            sa, sb = a * a / emax, b * b / emax
            return (
                k_norm
                / (4 * c * top * (m + 1))
                * (max(1 - sa, 0.0) ** (m + 1) - max(1 - sb, 0.0) ** (m + 1))
            )

        # NJOY の cm2lab の表は E' = 0 から始まる（単位基底の範囲が 0 からになる）
        if lo > 0:
            xs.append(0.0)
            ys.append(0.0)
        for r in roots:
            xs.append(r * r)
            ys.append(density(r))
    vals = []
    for x, y in zip(xs, ys):
        vals += [fmt(x), fmt(y)]
    out += rows(
        [
            fmt(0),
            fmt(e),
            f"{0:>11}",
            f"{0:>11}",
            f"{2 * len(xs):>11}",
            f"{len(xs):>11}",
        ],
        tail,
    )
    out += rows(vals, tail)
lines[first:end] = out
Path(dst).write_text("".join(lines))
