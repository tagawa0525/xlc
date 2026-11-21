"""Compare the two forms of the coarse-mesh diagonal-block solve on 69-group blocks (#0241).

The blocks come from the fine-group medium of VERA 2A (main 630459cba): A = V (diag(Σt + leak) − Sᵀ), the shape of
`cmfd::Coarse::homogenize`. The absorption is scaled down to make the blocks worse conditioned. The LU and inverse
forms copy xlc_numeric::lu::Lu (partial pivoting, row-wise inner products) and xlc_cmfd::Coarse::new/invert operation
by operation, so the float results are those of the Rust code up to the order of the sums.
The reference refines the LU solution with residuals computed exactly in fractions.
"""

import sys
from fractions import Fraction
from pathlib import Path

import msgpack
import numpy as np


def lu_factor(a, n):
    a = list(a)
    swaps = []
    for col in range(n):
        best = max(range(col, n), key=lambda i: abs(a[i * n + col]))
        swaps.append(best)
        if best != col:
            for j in range(n):
                a[best * n + j], a[col * n + j] = a[col * n + j], a[best * n + j]
        for row in range(col + 1, n):
            factor = a[row * n + col] / a[col * n + col]
            a[row * n + col] = factor
            for j in range(col + 1, n):
                a[row * n + j] -= factor * a[col * n + j]
    return a, swaps


def lu_solve(lu, swaps, n, b):
    b = list(b)
    for col, row in enumerate(swaps):
        if row != col:
            b[col], b[row] = b[row], b[col]
    for row in range(n):
        s = 0.0
        for j in range(row):
            s += lu[row * n + j] * b[j]
        b[row] -= s
    for row in reversed(range(n)):
        s = 0.0
        for j in range(row + 1, n):
            s += lu[row * n + j] * b[j]
        b[row] = (b[row] - s) / lu[row * n + row]
    return b


def inverse(lu, swaps, n):
    inv = [0.0] * (n * n)
    for j in range(n):
        e = [0.0] * n
        e[j] = 1.0
        col = lu_solve(lu, swaps, n, e)
        for i in range(n):
            inv[i * n + j] = col[i]
    return inv


def apply_inverse(inv, n, b):
    out = []
    for i in range(n):
        s = 0.0
        for j in range(n):
            s += inv[i * n + j] * b[j]
        out.append(s)
    return out


def reference(a, n, b, lu, swaps):
    exact = [Fraction(v) for v in a]
    bx = [Fraction(v) for v in b]
    x = lu_solve(lu, swaps, n, b)
    for _ in range(4):
        xf = [Fraction(v) for v in x]
        r = [bx[i] - sum(exact[i * n + j] * xf[j] for j in range(n)) for i in range(n)]
        d = lu_solve(lu, swaps, n, [float(v) for v in r])
        x = [float(xf[i] + Fraction(d[i])) for i in range(n)]
    return x


def main(path):
    state = msgpack.unpackb(Path(path).read_bytes(), raw=False, strict_map_key=False)[
        "states"
    ][0]
    total = np.array(state["total"])
    scatter = np.array(state["scatter"][0])  # [from][to]
    chi = np.array(state["chi"])
    n = len(total)
    volume = 1.26**2
    out_scatter = scatter.sum(axis=1)
    absorption = total - out_scatter
    rng = np.random.default_rng(0)
    print(f"groups {n}; min absorption/total {np.min(absorption / total):.3e}")
    print(
        "scale leak cond lu_err inv_err (max relative error over the groups, vs refined reference)"
    )
    for scale in [1.0, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6]:
        for leak in [0.0, 1.0, -0.5]:
            t = out_scatter + scale * absorption
            diffusion = 1.0 / (3.0 * t)
            diagonal = t + leak * 4.0 * diffusion / 1.26
            a = volume * (np.diag(diagonal) - scatter.T)
            flat = [float(v) for v in a.flatten()]
            cond = np.linalg.cond(a)
            lu, swaps = lu_factor(flat, n)
            inv = inverse(lu, swaps, n)
            worst = [0.0, 0.0]
            for trial in range(3):
                b = chi + (0.0 if trial == 0 else 0.5) * rng.random(n) * chi.max()
                b = [float(v) for v in b]
                ref = reference(flat, n, b, lu, swaps)
                for k, x in enumerate(
                    [lu_solve(lu, swaps, n, b), apply_inverse(inv, n, b)]
                ):
                    err = max(abs(x[i] - ref[i]) / abs(ref[i]) for i in range(n))
                    worst[k] = max(worst[k], err)
            print(f"{scale:g} {leak:g} {cond:.3e} {worst[0]:.3e} {worst[1]:.3e}")
            sys.stdout.flush()


if __name__ == "__main__":
    main(sys.argv[1])
