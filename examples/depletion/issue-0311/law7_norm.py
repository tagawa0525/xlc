"""MF6 MT16 の中性子の LAW=7 の、入射エネルギーごとの ∫∫ f dμ dE'（ENDF-102 の式 6.2 で 1）を、xs と同じ補間則で求める。

余弦は lin-lin（台形）かヒストグラム（左の表）、E' は lin-lin（台形）かヒストグラム（左の点の値 × 幅）。
使い方: python3 -I law7_norm.py <ENDF>...
"""

import sys
from pathlib import Path


def num(s):
    s = s.strip()
    if not s:
        return 0.0
    for i in range(1, len(s)):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


class Reader:
    def __init__(self, lines):
        self.lines, self.i = lines, 0

    def cont(self):
        line = self.lines[self.i]
        self.i += 1
        return [num(line[k * 11 : (k + 1) * 11]) for k in range(6)]

    def values(self, n):
        out = []
        while len(out) < n:
            out += self.cont()
        return out[:n]


def table_integral(laws, breaks, x, y):
    total = 0.0
    for k in range(len(x) - 1):
        law = laws[next(j for j, b in enumerate(breaks) if b >= k + 2)]
        width = x[k + 1] - x[k]
        total += y[k] * width if law == 1 else 0.5 * (y[k] + y[k + 1]) * width
    return total


for path in sys.argv[1:]:
    lines = [
        line for line in Path(path).read_text().splitlines() if line[70:75] == " 6 16"
    ]
    r = Reader(lines)
    r.cont()
    head = r.cont()
    nr, np_ = int(head[4]), int(head[5])
    r.values(2 * nr)
    r.values(2 * np_)
    assert int(head[3]) == 7, head
    outer = r.cont()
    r.values(2 * int(outer[4]))
    worst = (0.0, None, None)
    for _ in range(int(outer[5])):
        inner = r.cont()
        ranges = [int(v) for v in r.values(2 * int(inner[4]))]
        mu_breaks, mu_laws = ranges[0::2], ranges[1::2]
        mu, g = [], []
        for _ in range(int(inner[5])):
            t = r.cont()
            rr = [int(v) for v in r.values(2 * int(t[4]))]
            pts = r.values(2 * int(t[5]))
            mu.append(t[1])
            g.append(table_integral(rr[1::2], rr[0::2], pts[0::2], pts[1::2]))
        total = table_integral(mu_laws, mu_breaks, mu, g)
        if abs(total - 1) > abs(worst[0]):
            worst = (total - 1, inner[1], total)
    print(
        f"{path}: max |sum - 1| {abs(worst[0]):.3e} at E = {worst[1]:.6e} eV (sum {worst[2]:.9f})"
    )
