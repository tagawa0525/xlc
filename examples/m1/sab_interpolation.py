"""H in H2O の S(α,β) の散乱断面積を、ENDF の補間則（ln S を α, β について線形）と THERMR の補間（thermr.f90 の
`terpq`: ln S を 3 点の 2 次式、差が 2 を超える区間は線形）で求めて比べる（ADR-0133、#0005 の 4〜10 eV の差の候補）。
表の外と S = 0 の扱いは、どちらも THERMR の `sig` と同じ短時間衝突近似（SCT）。比だけを出すので定数は省く。

使い方: python3 sab_interpolation.py <tsl-HinH2O.endf> <温度 K> <入射エネルギー eV ...>
"""

import sys

import endf
import numpy as np

KT0 = 0.0253
BK = 8.617333262e-5
FLAG = -225.0


def load(path: str, temperature: float):
    d = endf.Material(path).section_data[7, 4]
    awr = d["B"][2]
    beta = np.array([b[0]["beta"] for b in d["beta_data"]])
    temps = [t["T"] for t in d["beta_data"][0]]
    k = next(i for i, t in enumerate(temps) if abs(t - temperature) < 1e-6)
    first = d["beta_data"][0][0]["S"]
    alpha = np.asarray(first.x)
    values = lambda entry: np.asarray(entry["S"].y if k == 0 else entry["S"])
    s = np.array([values(b[k]) for b in d["beta_data"]])  # [β][α]
    t_eff = d["Teff"][0](temperature)
    with np.errstate(divide="ignore"):
        ln_s = np.where(s > 0, np.log(np.maximum(s, 1e-300)), FLAG)
    return awr, np.array(alpha), beta, ln_s, t_eff


def sct(a, b, kt, kt_eff):
    """THERMR の SCT の対称な S（α, |β| は物質の温度の値）。"""
    return np.exp(-((a - b) ** 2) * kt / (4 * a * kt_eff) - b / 2) / np.sqrt(
        4 * np.pi * a * kt_eff / kt
    )


def terpq(x1, y1, x2, y2, x3, y3, x):
    """thermr.f90 の terpq（配列）。"""
    lin12 = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    lin23 = y2 + (y3 - y2) * (x - x2) / (x3 - x2)
    b = (y2 - y1) * (x3 - x1) / ((x2 - x1) * (x3 - x2)) - (y3 - y1) * (x2 - x1) / (
        (x3 - x1) * (x3 - x2)
    )
    c = (y3 - y1) / ((x3 - x1) * (x3 - x2)) - (y2 - y1) / ((x2 - x1) * (x3 - x2))
    quad = y1 + b * (x - x1) + c * (x - x1) ** 2
    corner = (np.abs(y1 - y2) > 2) | (np.abs(y2 - y3) > 2)
    y = np.where(corner, np.where(x < x2, lin12, lin23), quad)
    # x < x1: ln-lin（y1 > y2 なら y1）。x > x3: 線形（y3 > y2 なら y3）
    y = np.where(
        x < x1,
        np.where(y1 > y2, y1, y1 + (y2 - y1) * np.log(x / x1) / np.log(x2 / x1)),
        y,
    )
    y = np.where(x > x3, np.where(y3 > y2, y3, lin23), y)
    return np.maximum(y, FLAG)


def ln_s_at(method, alpha, beta, ln_s, a, b):
    """表の中の点 (a, b)（表の kT0 の値）の ln S と、THERMR の条件で SCT にするかどうか。"""
    na, nb = len(alpha), len(beta)
    ia = np.clip(np.searchsorted(alpha, a, side="right") - 1, 0, na - 2)
    ib = np.clip(np.searchsorted(beta, b, side="right") - 1, 0, nb - 2)
    corners = np.stack(
        [ln_s[ib, ia], ln_s[ib, ia + 1], ln_s[ib + 1, ia], ln_s[ib + 1, ia + 1]]
    )
    zero = (corners <= FLAG).any(axis=0)
    if method == "endf":
        u = (b - beta[ib]) / (beta[ib + 1] - beta[ib])
        v = (a - alpha[ia]) / (alpha[ia + 1] - alpha[ia])
        lo = corners[0] + v * (corners[1] - corners[0])
        hi = corners[2] + v * (corners[3] - corners[2])
        y = lo + u * (hi - lo)
        y = np.where(zero, -np.inf, y)
    else:
        ja = np.minimum(ia, na - 3)
        jb = np.minimum(ib, nb - 3)
        rows = []
        for k in range(3):
            r = jb + k
            rows.append(
                terpq(
                    alpha[ja],
                    ln_s[r, ja],
                    alpha[ja + 1],
                    ln_s[r, ja + 1],
                    alpha[ja + 2],
                    ln_s[r, ja + 2],
                    a,
                )
            )
        y = terpq(beta[jb], rows[0], beta[jb + 1], rows[1], beta[jb + 2], rows[2], b)
        y = np.where(y <= FLAG, -np.inf, y)
    return y, zero


def sigma(method, e, awr, alpha, beta, ln_s, kt, kt_eff):
    """∫ dE' e^(−β/2) ∫ S dα（α は物質の温度の値）。"""
    e_out = np.concatenate(
        [np.linspace(1e-6, e, 3001)[1:], np.linspace(e, e + 2.0, 1001)[1:]]
    )
    total = []
    for x in e_out:
        lo = (np.sqrt(e) - np.sqrt(x)) ** 2 / (awr * kt)
        hi = (np.sqrt(e) + np.sqrt(x)) ** 2 / (awr * kt)
        a = np.geomspace(max(lo, 1e-9), hi, 4001)
        bb = (x - e) / kt
        b = abs(bb)
        at, bt = a * kt / KT0, b * kt / KT0
        inside = (at <= alpha[-1]) & (bt <= beta[-1])
        s = sct(a, b, kt, kt_eff) * np.exp(-bb / 2)
        if bt <= beta[-1]:
            y, zero = ln_s_at(
                method,
                alpha,
                beta,
                ln_s,
                np.minimum(at, alpha[-1]),
                np.full_like(at, bt),
            )
            use_sct = zero & ((at * awr >= 30) | (bt >= 30))
            # (kT0/kT) ∫ S dα_表 = ∫ S dα
            s = np.where(inside & ~use_sct, np.exp(y - bb / 2), s)
        total.append(np.trapezoid(s, a))
    return np.trapezoid(total, e_out)


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    temperature = float(sys.argv[2])
    awr, alpha, beta, ln_s, t_eff = load(sys.argv[1], temperature)
    kt, kt_eff = BK * temperature, BK * t_eff
    for e in map(float, sys.argv[3:]):
        ours = sigma("endf", e, awr, alpha, beta, ln_s, kt, kt_eff)
        thermr = sigma("thermr", e, awr, alpha, beta, ln_s, kt, kt_eff)
        print(f"E {e:6.3f} eV: THERMR / ENDF − 1 = {thermr / ours - 1:+.4%}")


if __name__ == "__main__":
    main()
