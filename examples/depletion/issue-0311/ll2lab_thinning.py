"""NJOY GROUPR の ll2lab（LAW=7 を実験室系の Legendre にする）を、Be-9 の MF6 MT16 の中性子の 1 つの入射エネルギーで真似て、
ENDF のとおりのヒストグラムの g(E') = ∫ f(μ, E') dμ と、群ごとの割合を比べる。

ll2lab の点を間引く判定（前の点が要るか）の補間の重みが入れ替わっている版（NJOY のとおり）と、正しい重みの版を並べる。
使い方: python3 -I ll2lab_thinning.py <ENDF> <入射エネルギー eV>
"""

import sys
from itertools import pairwise
from pathlib import Path

import numpy as np

WIMS = [
    1.0e-5,
    0.005,
    0.01,
    0.015,
    0.02,
    0.025,
    0.03,
    0.035,
    0.042,
    0.05,
    0.058,
    0.067,
    0.08,
    0.1,
    0.14,
    0.18,
    0.22,
    0.25,
    0.28,
    0.3,
    0.32,
    0.35,
    0.4,
    0.5,
    0.625,
    0.78,
    0.85,
    0.91,
    0.95,
    0.972,
    0.996,
    1.02,
    1.045,
    1.071,
    1.097,
    1.123,
    1.15,
    1.3,
    1.5,
    2.1,
    2.6,
    3.3,
    4.0,
    9.877,
    15.968,
    27.7,
    48.052,
    75.501,
    148.728,
    367.262,
    906.898,
    1425.1,
    2239.45,
    3519.1,
    5530.0,
    9118.0,
    1.503e4,
    2.478e4,
    4.085e4,
    6.734e4,
    1.11e5,
    1.83e5,
    3.025e5,
    5.0e5,
    8.21e5,
    1.353e6,
    2.231e6,
    3.679e6,
    6.0655e6,
    1.0e7,
][::-1]


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

    def tab1(self):
        head = self.cont()
        nr, np_ = int(head[4]), int(head[5])
        ranges = [int(v) for v in self.values(2 * nr)]
        points = self.values(2 * np_)
        return head, ranges[1::2], points[0::2], points[1::2]


def law7_at(path, energy):
    lines = [
        l
        for l in Path(path).read_text().splitlines(keepends=True)
        if l[70:75] == " 6 16"
    ]
    r = Reader(lines)
    r.cont()
    head, _, _, _ = r.tab1()
    assert int(head[3]) == 7
    outer = r.cont()
    r.values(2 * int(outer[4]))
    for _ in range(int(outer[5])):
        inner = r.cont()
        codes = [int(v) for v in r.values(2 * int(inner[4]))][1::2]
        assert codes == [2], codes
        tables = []
        for _ in range(int(inner[5])):
            h, laws, x, y = r.tab1()
            assert laws == [1], laws
            tables.append((h[1], np.array(x), np.array(y)))
        if inner[1] == energy:
            return tables
    raise SystemExit(f"no incident energy {energy}")


def hist(x, y, e):
    """ヒストグラムの表の e での値（区間の左の点の値、範囲の外は 0）"""
    if e < x[0] or e >= x[-1]:
        return 0.0
    return y[np.searchsorted(x, e, side="right") - 1]


def exact_fractions(tables):
    mu = np.array([t[0] for t in tables])
    w = np.zeros(len(mu))
    w[:-1] += 0.5 * np.diff(mu)
    w[1:] += 0.5 * np.diff(mu)
    grid = np.unique(np.concatenate([t[1] for t in tables]))
    out = np.zeros(69)
    for a, b in pairwise(grid):
        g = sum(wi * hist(x, y, a) for wi, (_, x, y) in zip(w, tables))
        add_linear(out, a, b, g, g)
    return out


def add_linear(out, a, b, ya, yb):
    """[a, b] で線形な密度を群に振り分ける（範囲の外は端の群）"""
    cuts = sorted({a, b, *[c for c in WIMS if a < c < b]})
    for lo, hi in pairwise(cuts):
        mid = 0.5 * (lo + hi)
        y = lambda e: ya + (yb - ya) * (e - a) / (b - a)
        g = next(
            (k for k in range(69) if WIMS[k + 1] <= mid <= WIMS[k]),
            0 if mid > WIMS[0] else 68,
        )
        out[g] += 0.5 * (hi - lo) * (y(lo) + y(hi))


def ll2lab(tables, swapped, nl=4):
    """ll2lab の P0..P3 の点（ep, term）を、間引きを含めて作る"""
    qp, qw = np.polynomial.legendre.leggauss(8)
    legendre = [np.polynomial.legendre.Legendre.basis(l)(qp) for l in range(nl)]
    mu = np.array([t[0] for t in tables])
    points = []
    ep, emax = 0.0, 1.0e10
    while ep < emax * (1 - 1e-10):
        epnext = emax
        fmu = []
        for _, x, y in tables:
            # terpa: 値と次の点。ヒストグラム（idis=1）は次の点の 0.99999 倍で止まる
            if ep < x[0]:
                f, nxt = 0.0, x[0]
            elif ep >= x[-1] * 1.00001:
                f, nxt = 0.0, 1e12
            elif ep >= x[-1]:
                f, nxt = y[-1], (1.00001**2 * x[-1] if y[-1] > 0 else 1e12)
            else:
                k = np.searchsorted(x, ep, side="right") - 1
                f, nxt = y[k], x[k + 1]
                if ep < 0.99999 * nxt:
                    nxt = 0.99999 * nxt
            fmu.append(f)
            epnext = min(epnext, nxt)
        fmu = np.array(fmu)
        j = np.searchsorted(mu, qp, side="right")
        frac = (qp - mu[j - 1]) / (mu[j] - mu[j - 1])
        ff = (1 - frac) * fmu[j - 1] + frac * fmu[j]
        term = np.array([np.sum(ff * p * qw) for p in legendre])
        total = np.sum(np.diff(mu) * (fmu[1:] + fmu[:-1]) / 2)
        if term[0] != 0:
            term *= total / term[0]
        points.append([ep, *term])
        if len(points) > 3:
            (x0, *y0), (x1, *y1), (x2, *y2) = points[-3:]
            f1 = (x1 - x0) / (x2 - x0)
            keep = False
            for l in range(nl):
                tmid = (
                    f1 * y0[l] + (1 - f1) * y2[l]
                    if swapped
                    else (1 - f1) * y0[l] + f1 * y2[l]
                )
                if abs(tmid - y1[l]) > 0.005 * abs(tmid) + 1e-12:
                    keep = True
            if not keep:
                del points[-2]
        ep = epnext
    return np.array(points)


def fractions(points):
    out = np.zeros(69)
    for (a, ya, *_), (b, yb, *_) in pairwise(points):
        if b > a:
            add_linear(out, a, b, ya, yb)
    return out


tables = law7_at(sys.argv[1], float(sys.argv[2]))
exact = exact_fractions(tables)
for name, swapped in [("NJOY (swapped weights)", True), ("corrected weights", False)]:
    approx = fractions(ll2lab(tables, swapped))
    diff = np.sum(np.abs(approx - exact)) / np.sum(exact)
    mean = lambda p: (
        np.sum(p * np.sqrt(np.array(WIMS[:-1]) * np.array(WIMS[1:]))) / np.sum(p)
    )
    print(
        f"{name}: sum |diff| / sum {diff:.3e}; sum exact {exact.sum():.6f}, approx {approx.sum():.6f}; "
        f"group-mean energy exact {mean(exact):.4e}, approx {mean(approx):.4e} eV"
    )
