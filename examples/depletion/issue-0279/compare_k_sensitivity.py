"""k の数密度の感度（lattice の `[sensitivity]`。ADR-0247 の段 2）を、同じ入力の直接の摂動の中心差分と核種ごとに比べる（#0279）。

一般化摂動の値 S（`relative`）と、数密度を (1 ± ε) 倍と (1 ± ε/2) 倍した k の中心差分の Richardson の外挿 R（`direct`）の差を、
許容 max(2 e, 1e-4 |R|, 3 τ / ε) と比べる。e = |D(ε/2) − D(ε)| / 3 は D(ε/2) の刻みの誤差の Richardson の見積もり（#0279 の計画の
許容の目安）、3 τ / ε は k の判定 τ の誤差が外挿 R = (4 D(ε/2) − D(ε)) / 3 に入る大きさ（各 k の誤差が τ のとき
(4/3)(2τ/ε) + (1/3)(τ/ε)）。感度の小さい核種では k の変化が判定の誤差に近く、参照の側がこの項で決まる。

使い方: python3 compare_k_sensitivity.py <結果の TOML>...（どれかの核種が許容を超えれば終了コード 1）
"""

import sys
from pathlib import Path

import tomllib

RELATIVE = 1.0e-4


def compare(report: dict) -> list[dict]:
    """核種ごとの比較の行"""
    epsilon = report["epsilon"]
    tau = report["k_tolerance"]
    rows = []
    for n in report["nuclides"]:
        direct = n["direct"]
        coarse, fine = direct["central"]
        reference = direct["richardson"]
        estimate = abs(fine - coarse) / 3.0
        floor = 3.0 * tau / epsilon
        allowed = max(2.0 * estimate, RELATIVE * abs(reference), floor)
        difference = n["relative"] - reference
        rows.append(
            {
                "name": n["name"],
                "adjoint": n["relative"],
                "direct": reference,
                "difference": difference,
                "relative": difference / reference
                if reference != 0.0
                else float("inf"),
                "allowed": allowed,
                "limited_by": max(
                    [
                        ("step", 2.0 * estimate),
                        ("relative", RELATIVE * abs(reference)),
                        ("k tolerance", floor),
                    ],
                    key=lambda p: p[1],
                )[0],
                "ok": abs(difference) <= allowed,
            }
        )
    return rows


def main(paths: list[str]) -> int:
    failed = False
    for path in paths:
        report = tomllib.loads(Path(path).read_text())
        rows = compare(report)
        print(
            f"{path}: k {report['k']:.10f}, adjoint k {report['adjoint_k']:.10f}, eps {report['epsilon']}"
        )
        print(
            f"{'nuclide':8s} {'S (adjoint)':>14s} {'R (direct)':>14s} {'rel diff':>10s} {'allowed':>10s}  limited by"
        )
        for r in rows:
            mark = "" if r["ok"] else "  OUTSIDE"
            print(
                f"{r['name']:8s} {r['adjoint']:14.6e} {r['direct']:14.6e} {r['relative']:10.2e} "
                f"{r['allowed'] / abs(r['direct']):10.2e}  {r['limited_by']}{mark}"
            )
        worst = max(rows, key=lambda r: abs(r["relative"]))
        print(f"largest relative difference: {worst['name']} {worst['relative']:.2e}")
        seconds = report["seconds"]
        print(
            f"seconds: forward {seconds['forward']:.3f}, adjoint {seconds['adjoint']:.3f} "
            f"(ratio {seconds['adjoint'] / seconds['forward']:.2f}), products {seconds['products']:.3f}"
        )
        failed |= not all(r["ok"] for r in rows)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
