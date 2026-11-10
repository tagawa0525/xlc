"""直積の分岐のテーブルの格子点の間の内挿の誤差を測る（ADR-0239 の確認）。連続の軸（ホウ素、減速材、燃料の温度）の内側の格子点
ごとに、その軸の両隣の格子点の状態の少数群定数を core と同じ式で内挿し、格子点の状態と比べる。測った誤差から、本番の間隔
（隣り合う格子点の区間）の中点の誤差を見込む。

使い方: python3 compare_grid_interpolation.py --json <要約の出力先> <直積のテーブル>
        python3 compare_grid_interpolation.py --report <要約の JSON>

- 内挿: core の `bracket`（crates/core/src/state.rs）と同じく、軸の内挿の変数 g（テーブルの軸の `transform`。`linear` は値、`sqrt` は
  √値）について隣の 2 点で線形にする。格子点 x_i を両隣 x_{i−1}、x_{i+1} から重み w = (g_i − g_{i−1}) / (g_{i+1} − g_{i−1}) で
  (1 − w) f_{i−1} + w f_{i+1} とする。ほかの軸は格子点のまま（直積の状態）で、燃焼度の点ごとに、ほかの軸のすべての点の組で比べる
- 誤差: Δρ = (1/k_格子点 − 1/k_内挿) × 1e5 pcm（k∞ は少数群定数の無限媒質の固有値。compare_burnup_interpolation.py の `k_inf`）。
  格子点ごと、燃焼度の点ごとに、|Δρ| の最大とその状態を書く（`centers`）。これは両隣の間隔 h1 + h2 の区間の誤差で、本番の
  間隔の誤差ではない
- 見込み（測定でない）: 区間の中で f の g についての 2 階の微分が一定なら、線形の内挿の誤差は (1/2) f''(g − a)(g − b) で、
  格子点 i の誤差は (1/2) f'' h1 h2（h1 = g_i − g_{i−1}、h2 = g_{i+1} − g_i）、幅 h の区間の中点の誤差は f'' h² / 8。比の
  h² / (4 h1 h2) を格子点 i の |Δρ| の最大に掛けて、i を端に持つ区間（幅 h1 と h2）の中点の誤差を見込む（`intervals` の `from`）。
  内側の区間は両端の格子点の見込みを持ち、大きい方を `max_abs_pcm` に書く。軸の端の区間は内側の端の格子点の見込みだけ。Δρ が
  少数群定数の差について線形とみなす
- `inputs` は読んだファイルのパス（実行した cwd から）と md5
"""

import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_branch_cross_terms import table_grid
from compare_burnup_interpolation import k_inf, unvector, vector

TRANSFORMS = {"linear": lambda x: x, "sqrt": math.sqrt}


def axis_transforms(table: dict) -> list:
    """燃焼度の軸を除く軸の内挿の変数（core と同じく、テーブルの軸の `transform`）"""
    out = []
    for a in table["axes"]:
        if a["name"] == "burnup":
            continue
        name = a.get("transform")
        if name not in TRANSFORMS:
            raise SystemExit(
                f"{a['name']}: the axis needs a transform of {sorted(TRANSFORMS)}, not {name}"
            )
        out.append(TRANSFORMS[name])
    return out


def summarize(table: dict) -> dict:
    # table_grid は公称の値を引数に取る（測る量には使わない）。lattice は、入力の物質の状態がどの軸でも格子点にあるときだけ軸に公称の
    # 値を書く（crates/lattice/src/branches.rs）ので、公称の値の無いテーブルはここで止まる。制御棒の軸（`points` が無い）は
    # table_grid が止める
    state_axes = [a for a in table["axes"] if a["name"] != "burnup" and "points" in a]
    if any("nominal" not in a for a in state_axes):
        raise SystemExit("the table needs the nominal value of every state axis")
    grid_of = table_grid("table", table, {a["name"]: a["nominal"] for a in state_axes})
    transforms = axis_transforms(table)
    burnup, names, points = grid_of["burnup"], grid_of["names"], grid_of["points"]
    burnups = burnup if burnup is not None else [None]
    per_point = grid_of["per_point"]
    indices = list(itertools.product(*[range(len(p)) for p in points]))
    flat = {index: n for n, index in enumerate(indices)}
    vectors = []
    truths = []
    for b in range(len(burnups)):
        chunk = table["states"][b * per_point : (b + 1) * per_point]
        vectors.append([vector(s) for s in chunk])
        truths.append([(s, k_inf(s)) for s in chunk])
    centers = []
    gaps = {}
    for a, name in enumerate(names):
        g = [transforms[a](x) for x in points[a]]
        for i in range(1, len(g) - 1):
            w = (g[i] - g[i - 1]) / (g[i + 1] - g[i - 1])
            out = []
            for b in range(len(burnups)):
                worst = None
                for index in indices:
                    if index[a] != i:
                        continue
                    low = flat[index[:a] + (i - 1,) + index[a + 1 :]]
                    high = flat[index[:a] + (i + 1,) + index[a + 1 :]]
                    state, k_grid = truths[b][flat[index]]
                    rebuilt = unvector(
                        (1 - w) * vectors[b][low] + w * vectors[b][high], state
                    )
                    rho = (1 / k_grid - 1 / k_inf(rebuilt)) * 1e5
                    if worst is None or abs(rho) > abs(worst[0]):
                        worst = (rho, index)
                out.append(
                    {
                        "burnup": burnups[b],
                        "max_abs_pcm": abs(worst[0]),
                        "worst_pcm": worst[0],
                        "worst": {
                            n: points[k][worst[1][k]] for k, n in enumerate(names)
                        },
                    }
                )
            centers.append(
                {
                    "axis": name,
                    "point": points[a][i],
                    "neighbours": [points[a][i - 1], points[a][i + 1]],
                    "weight": w,
                    "points": out,
                }
            )
            h1, h2 = g[i] - g[i - 1], g[i + 1] - g[i]
            for left, h in ((i - 1, h1), (i, h2)):
                factor = h * h / (4 * h1 * h2)
                gaps.setdefault((a, left), []).append(
                    {
                        "center": points[a][i],
                        "factor": factor,
                        "max_abs_pcm": [p["max_abs_pcm"] * factor for p in out],
                    }
                )
    intervals = [
        {
            "axis": names[a],
            "interval": [points[a][left], points[a][left + 1]],
            "from": sources,
            "max_abs_pcm": [
                max(s["max_abs_pcm"][b] for s in sources) for b in range(len(burnups))
            ],
        }
        for (a, left), sources in sorted(gaps.items())
    ]
    return {
        "format": "xlc-compare-grid-interpolation",
        "burnup": burnup,
        "axes": {
            n: {"points": p, "transform": a["transform"]}
            for n, p, a in zip(
                names,
                points,
                [a for a in table["axes"] if a["name"] != "burnup"],
                strict=True,
            )
        },
        "table_k_max_abs_pcm": float(
            np.max(
                [abs(1 / s["k_inf"] - 1 / k) * 1e5 for row in truths for s, k in row]
            )
        ),
        "centers": centers,
        "intervals": intervals,
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"k∞ の解き直しとテーブルの k_inf の差の最大 {summary['table_k_max_abs_pcm']:.4f} pcm",
        "格子点を両隣から内挿した |Δρ| の最大（pcm。燃焼度の点ごと）:",
    ]
    for c in summary["centers"]:
        values = "、".join(
            f"{p['worst_pcm']:+.2f}（{p['worst']}）" for p in c["points"]
        )
        lines.append(
            f"  {c['axis']} {c['point']:g}（{c['neighbours'][0]:g} と {c['neighbours'][1]:g}、重み {c['weight']:.3f}）: {values}"
        )
    lines.append("区間の中点の見込み（pcm。燃焼度の点ごと。測定でない）:")
    for i in summary["intervals"]:
        sources = "、".join(
            f"{s['center']:g} から × {s['factor']:.3f} で "
            + "/".join(f"{v:.2f}" for v in s["max_abs_pcm"])
            for s in i["from"]
        )
        lines.append(
            f"  {i['axis']} {i['interval'][0]:g}〜{i['interval'][1]:g}: {sources}"
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 3 and args[0] == "--json":
        path = Path(args[2])
        data = path.read_bytes()
        summary = summarize(msgpack.unpackb(data))
        summary["inputs"] = [{"path": str(path), "md5": hashlib.md5(data).hexdigest()}]
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
