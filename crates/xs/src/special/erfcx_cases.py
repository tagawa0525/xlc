"""erfcx_cases.rs と special/tests.rs の多倍長の参照値（x, hi, lo）の行を、Rust の配列の要素の形で書き出す。

hi は 60 桁の値を倍精度に丸めた値、lo は残りを丸めた値。mpmath は devShell に無いので、nixpkgs の Python で実行する:

    nix shell --impure --expr 'let p = (builtins.getFlake "nixpkgs").legacyPackages.x86_64-linux;
      in p.python3.withPackages (ps: [ps.mpmath])' -c python3 -I erfcx_cases.py

ERFCX_BITS は参照値でなく実装のビットなので、ここでは作らない（テストが違うビットを示す）。
"""

import mpmath

mpmath.mp.dps = 60

# 負の値、区間の端（0.46875、4、50）、ADR-0211 の計測で誤差が大きかった点
ERFCX_POINTS = [
    -0.25, -1.0, -2.0, -3.0, -10.0,
    0.0, 0.25, 0.46875, 0.5, 0.7837817412143546, 1.0, 1.4793300144592165, 1.5, 2.0, 2.171084449465437,
    2.4752524646965584, 2.4999999999999996, 2.5, 2.5124651723206233, 2.5846515745793606, 2.6050776195927656, 2.75,
    3.0, 3.5, 3.9999999999999996, 4.0, 4.000000000000001, 4.44892545588937, 4.463089702624895, 5.0, 6.0, 8.0,
    9.57370480820499, 10.0, 20.0, 36.82194217585145, 50.0, 100.0, 602.8467785386605, 1e3, 1e4, 1e6, 1e8, 1e12,
    1e100, 1e200,
]  # fmt: skip
ERF_POINTS = [
    0.25, 0.46875, 0.47280819081272235, 1.0, 2.0, 2.4160745503842005, 2.4999999999999996, 3.0, 5.0, 6.0,
]  # fmt: skip
ERFC_POINTS = [
    0.46875, 1.0, 2.0, 2.0462943283432917, 2.4810855475278304, 2.4999999999999996, 3.0, 5.0,
]  # fmt: skip


def erfcx(x):
    x = mpmath.mpf(x)
    if x > 1e150:
        # mpmath の erfc が扱えない。漸近展開の 2 項（次の項は 3/(4x⁴) の相対の大きさ）
        return (1 - 1 / (2 * x * x)) / (mpmath.sqrt(mpmath.pi) * x)
    return mpmath.erfc(x) * mpmath.exp(x * x)


def rows(name, f, points):
    print(f"// {name}")
    for x in points:
        v = f(x)
        hi = float(v)
        lo = float(v - mpmath.mpf(hi))
        print(f"    ({x!r}, {hi!r}, {lo!r}),")


if __name__ == "__main__":
    rows("ERFCX_REFERENCE", erfcx, ERFCX_POINTS)
    rows("ERF_REFERENCE", lambda x: mpmath.erf(mpmath.mpf(x)), ERF_POINTS)
    rows("erfc（tests.rs）", lambda x: mpmath.erfc(mpmath.mpf(x)), ERFC_POINTS)
