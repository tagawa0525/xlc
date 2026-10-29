"""ミクロ燃焼の線形の補正（ADR-0216）が、追う核種の組成の差の輸送の k∞ をどれだけ再現するかを、燃焼度の隣り合う点の組ごとに出す。

    python3 compare_linearization.py --json linearization.json   # このディレクトリの出力から JSON を書く
    python3 compare_linearization.py --report linearization.json  # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは 2 つの分岐のテーブル（`xlc-table`。燃焼度の軸と 1 点の分岐の軸を持ち、燃焼度の点ごとに状態 1 つ）: 2A の履歴の点を解いた `table-nominal.msgpack`
（`micro` を持つ）と、点 k の組成を点 k - 1 の組成のうち追う核種だけを点 k の値にしたもので解いた `table-shifted.msgpack`
（make_shifted.py）。点の組 (k - 1, k) ごとに:

- 輸送: ρ(shifted の点 k の k∞) − ρ(nominal の点 k - 1 の k∞)。追う核種の組成の差だけの効き
- 線形: ρ(予測) − ρ(nominal の点 k - 1 の少数群定数の 2 群の k)。予測は、点 k - 1 の少数群定数（全断面積、散乱行列、νΣf）に
  Σᵢ σᵢ(k - 1)(N̄ᵢ(k) − N̄ᵢ(k - 1)) を足した定数の 2 群の無限媒質の k
- 誤差: 線形 − 輸送（2 つの基準が同じなら ρ(予測) − ρ(shifted の点 k の k∞)。正なら線形の補正は反応度を高く見積もる）

ρ は反応度 1 − 1/k で、差は pcm（Δρ = (1/k_a − 1/k_b) × 1e5、a から b へ）。2 群の k は、少数群定数から無限媒質の固有値として求め、テーブルの
k∞ との差（`self_check_pcm`。nominal のすべての状態の最大）を出す（外周が鏡面反射の格子の少数群定数は格子の k∞ を保つ）。
"""

import argparse
import hashlib
import json
from pathlib import Path

import msgpack

# 燃焼度 0 の組成は 2 つの履歴で同じなので、k∞ は同じ問題の値（同じ入力と実行ファイル）。違えば別の条件で解いたテーブル
FRESH_TOLERANCE_PCM = 1e-6
# 組成の差の相対値を出す核種
NUCLIDES = ("U-235", "Pu-239", "Pu-241", "Xe-135", "Sm-149")
FILES = {"nominal": "table-nominal.msgpack", "shifted": "table-shifted.msgpack"}


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def solve(matrix, vector):
    """小さな連立 1 次方程式（部分ピボットの Gauss の消去）"""
    n = len(vector)
    a = [list(row) + [b] for row, b in zip(matrix, vector)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(a[r][c]))
        a[c], a[p] = a[p], a[c]
        for r in range(c + 1, n):
            f = a[r][c] / a[c][c]
            for j in range(c, n + 1):
                a[r][j] -= f * a[c][j]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - sum(a[r][j] * x[j] for j in range(r + 1, n))) / a[r][r]
    return x


def k_inf(state):
    """少数群定数の無限媒質の k: (diag(Σt) − Sᵀ) φ = χ を解き、k = νΣf · φ（χ に比例する核分裂源は 1 つの固有値だけを持つ）"""
    total, scatter = state["total"], state["scatter"]
    groups = len(total)
    removal = [
        [(total[g] if g == h else 0.0) - scatter[h][g] for h in range(groups)]
        for g in range(groups)
    ]
    flux = solve(removal, state["chi"])
    return sum(n * f for n, f in zip(state["nu_fission"], flux))


def corrected(state, sigmas, deltas):
    """状態の少数群定数に、核種ごとの微視的断面積 `sigmas` (barn) × 数密度の差 `deltas` (1/(barn cm)) を足したもの"""
    out = dict(state)
    for key in ("total", "nu_fission"):
        out[key] = [
            value + sum(s[key][g] * d for s, d in zip(sigmas, deltas))
            for g, value in enumerate(state[key])
        ]
    out["scatter"] = [
        [
            value + sum(s["scatter"][g][h] * d for s, d in zip(sigmas, deltas))
            for h, value in enumerate(row)
        ]
        for g, row in enumerate(state["scatter"])
    ]
    return out


def burnups(table):
    """燃焼度の点。燃焼度の軸を最初に持ち、ほかの軸（分岐の軸）は 1 点で、点ごとに状態 1 つのテーブルだけ"""
    axes = table["axes"]
    single = all(len(a.get("points", a.get("values", []))) == 1 for a in axes[1:])
    if (
        not axes
        or axes[0]["name"] != "burnup"
        or not single
        or len(table["states"]) != len(axes[0]["points"])
    ):
        raise SystemExit(
            "the tables must have a burnup axis, branch axes of one point and a state per burnup"
        )
    return axes[0]["points"]


def rho(k):
    """反応度 1 − 1/k (pcm)"""
    return 1e5 * (1.0 - 1.0 / k)


def compare(nominal, shifted):
    if nominal.get("micro") is None:
        raise SystemExit("the nominal table has no micro data")
    points = burnups(nominal)
    if burnups(shifted) != points:
        raise SystemExit("the tables have different burnup points")
    micro = nominal["micro"]
    names = micro["nuclides"]
    fresh = rho(shifted["states"][0]["k_inf"]) - rho(nominal["states"][0]["k_inf"])
    if abs(fresh) > FRESH_TOLERANCE_PCM:
        raise SystemExit(f"the fresh points differ by {fresh} pcm")
    # 予測に使うのは nominal の少数群定数なので、その状態で確かめる（shifted は k∞ だけを使う）
    self_check = max(abs(rho(k_inf(s)) - rho(s["k_inf"])) for s in nominal["states"])
    rows = []
    for k in range(1, len(points)):
        before, after = micro["densities"][k - 1], micro["densities"][k]
        deltas = [b - a for a, b in zip(before, after)]
        base = nominal["states"][k - 1]
        predicted = k_inf(corrected(base, micro["states"][k - 1], deltas))
        transport = shifted["states"][k]["k_inf"]
        # どちらの効きも自分の基準（テーブルの k∞ と、同じ定数の 2 群の k）からの差なので、基準の収束の差（self_check）は入らない
        tracked = rho(transport) - rho(base["k_inf"])
        linear = rho(predicted) - rho(k_inf(base))
        rows.append(
            {
                "from": points[k - 1],
                "to": points[k],
                "tracked_pcm": tracked,
                "linear_pcm": linear,
                "error_pcm": linear - tracked,
                "density_change_percent": {
                    n: (
                        None
                        if before[names.index(n)] == 0.0
                        else 100.0 * deltas[names.index(n)] / before[names.index(n)]
                    )
                    for n in NUCLIDES
                    if n in names
                },
            }
        )
    return {
        "format": "xlc-compare-micro-linearization",
        "nuclides": len(names),
        "self_check_pcm": self_check,
        "rows": rows,
    }


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def report(summary):
    rows = summary["rows"]
    errors = [r["error_pcm"] for r in rows]
    worst = max(rows, key=lambda r: abs(r["error_pcm"]))
    shares = [
        r["linear_pcm"] / r["tracked_pcm"] for r in rows if abs(r["tracked_pcm"]) > 1.0
    ]
    return (
        f"追う核種 {summary['nuclides']}、点の組 {len(rows)}。線形の補正の誤差（線形 − 輸送）は "
        f"{min(errors):+.2f}〜{max(errors):+.2f} pcm（絶対値の最大は {worst['from']:g}→{worst['to']:g} MWd/kg の "
        f"{worst['error_pcm']:+.2f} pcm、そのときの輸送の効き {worst['tracked_pcm']:+.1f} pcm）。効きが 1 pcm を超える組で、線形の補正は "
        f"輸送の効きの {100 * min(shares):.2f}〜{100 * max(shares):.2f}% を再現した。2 群の k とテーブルの k∞ の差は最大 "
        f"{summary['self_check_pcm']:.2e} pcm"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", help="write the summary to this file")
    parser.add_argument("--report", help="print the report of this summary")
    args = parser.parse_args()
    if args.report:
        with open(args.report, encoding="utf-8") as f:
            print(report(json.load(f)))
        return
    here = Path(__file__).resolve().parent
    tables = {name: read(here / file) for name, file in FILES.items()}
    summary = compare(tables["nominal"], tables["shifted"])
    summary["inputs"] = {name: md5(here / file) for name, file in FILES.items()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
