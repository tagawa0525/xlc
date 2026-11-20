"""cases.rs の Coulomb の透過因子 P_l(η, ρ) = ρ / (F_l² + G_l²) の多倍長の参照値の行を、Rust の配列の要素の形で書き出す。

F_l と G_l は mpmath の coulombf と coulombg。40 桁と 80 桁で求めた値の相対差が 1e-30 を超えたら止める（mpmath の
超幾何級数の桁落ちを、桁を増やした値で確かめる）。mpmath は devShell に無いので、nixpkgs の Python で実行する:

    nix shell --impure --expr 'let p = (builtins.getFlake "nixpkgs").legacyPackages.x86_64-linux;
      in p.python3.withPackages (ps: [ps.mpmath])' -c python3 -I cases.py
"""

import mpmath

# (η, ρ)。荷電粒子の対を持つ評価（Cl-35 の (n,p)、Ca-40 の (n,p) と (n,α)）の共鳴エネルギーでの範囲（η 0.88〜17.9、
# ρ 0.18〜3.9）、Ca-40 の (n,p) で NJOY の近似の式が切り替わる 823 keV の点（η 5.66、ρ 0.567）、障壁の外（ρ が転回点
# η + √(η² + l(l+1)) より大きい）、転回点のまわり、η = 0（剛体球と同じ値）、障壁の深い下（P が 1e-300 の桁）
POINTS = [
    (0.0, 0.5), (0.0, 3.0),
    (0.1, 0.01), (0.1, 1.0),
    (0.883, 2.946), (3.179, 0.819), (3.3, 0.8),
    (8.143, 2.687), (8.6, 2.7), (5.666, 3.861),
    (11.513, 0.279), (17.942, 0.179), (5.66, 0.567),
    (5.0, 9.9), (5.0, 10.1), (5.0, 30.0), (2.0, 100.0),
    (60.0, 0.05), (110.0, 0.02),
]  # fmt: skip
# P が倍精度の最小の正の数より十分に小さく、0 になる点
UNDERFLOW_POINTS = [(300.0, 0.01), (600.0, 0.005)]


def penetrability(l, eta, rho, dps):
    with mpmath.workdps(dps):
        eta, rho = mpmath.mpf(eta), mpmath.mpf(rho)
        f = mpmath.coulombf(l, eta, rho)
        g = mpmath.coulombg(l, eta, rho)
        return rho / (f * f + g * g)


def checked(l, eta, rho):
    p = penetrability(l, eta, rho, 40)
    q = penetrability(l, eta, rho, 80)
    if abs(p / q - 1) > 1e-30:
        raise SystemExit(
            f"l={l} eta={eta} rho={rho}: 40 and 80 digits differ ({p} and {q})"
        )
    return q


if __name__ == "__main__":
    print("// COULOMB_REFERENCE")
    for l in range(5):
        for eta, rho in POINTS:
            print(f"    ({l}, {eta!r}, {rho!r}, {float(checked(l, eta, rho))!r}),")
    print("// COULOMB_UNDERFLOW")
    for l in (0, 4):
        for eta, rho in UNDERFLOW_POINTS:
            # 多倍長の値が最小の正の非正規化数（約 4.9e-324）の半分より小さいことを確かめる
            p = penetrability(l, eta, rho, 40)
            if p >= mpmath.mpf(2) ** -1075:
                raise SystemExit(f"l={l} eta={eta} rho={rho}: P={p} does not underflow")
            print(f"    ({l}, {eta!r}, {rho!r}),")
