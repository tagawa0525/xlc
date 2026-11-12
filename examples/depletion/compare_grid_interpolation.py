"""直積の分岐のテーブルの格子点の間の内挿の誤差を測る（ADR-0239 の確認）。連続の軸（ホウ素、減速材、燃料の温度）の内側の格子点
ごとに、その軸の両隣の格子点の状態の少数群定数を隣の 2 点の線形（ADR-0097 の式）で内挿し、格子点の状態と比べる。測った誤差から、本番の間隔
（隣り合う格子点の区間）の中点の誤差を見込む。

使い方: python3 compare_grid_interpolation.py --json <要約の出力先> <入力の TOML>
        python3 compare_grid_interpolation.py --report <要約の JSON>

入力の TOML は、直積のテーブル `table`（TOML のディレクトリからのパス）と、比べる内挿の方式 `[[scheme]]`（`axis`、`variable`、
`method`）。方式ごとに、内側の格子点を、その点を除いたその軸の格子点から内挿して比べる（`schemes`）:

- `variable`: 内挿の変数。`value`（軸の値。ppm か K）、`sqrt`（√値）、`log`（ln 値）、`density`（減速材の軸の格子点の密度。ADR-0229 の
  案 T2）
- `method`: `linear`（隣の 2 点の線形。core のホウ素と減速材の軸の式）、`hermite`（傾きを点と両隣の 3 点の放物線の微分にする 3 次 Hermite。core の
  燃焼度の軸（ADR-0192）と燃料の温度の軸（`sqrt`。ADR-0242）と同じ形で、燃焼度 0 の扱いは無い）、`lagrange`（区間の 2 点とその外側の 1 点ずつの 4 点の 3 次。端の区間は端の
  4 点）、`spline`（自然 3 次スプライン。値は軸のすべての点に依る）。どれも値は格子点の値の線形結合で、少数群定数の項目ごとに同じ
  重みで足す。点を除いた格子が方式の要る点の数（`linear` 2、`hermite` と `spline` 3、`lagrange` 4）に足りなければエラー

- 線形の式（`centers` と `intervals`）の内挿: ADR-0097 の式（core のホウ素と減速材の軸の `bracket`、crates/core/src/state.rs。燃料の温度の軸の
  core の式は、`schemes` の `sqrt` の `hermite`。ADR-0242）と同じく、軸の内挿の変数 g（テーブルの軸の `transform`。`linear` は値、`sqrt` は
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
- `schemes` は方式ごとの `centers` と同じ形の誤差（両隣でなく、その点を除いた格子の内挿。`linear` の両隣の重みは `centers` と
  同じ）。見込みは書かない（`intervals` は線形の内挿の誤差の 2 次の見込み）
- `inputs` は読んだファイル（TOML とテーブル）のパス（実行した cwd から）と md5
"""

import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib

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


VARIABLES = {
    "value": lambda axis: list(axis["points"]),
    "sqrt": lambda axis: [math.sqrt(x) for x in axis["points"]],
    "log": lambda axis: [math.log(x) for x in axis["points"]],
    # 密度は温度とともに下がるので、符号を変えて増える列にする（内挿の重みは変数の 1 次の変換で変わらない）
    "density": lambda axis: [-d for d in axis["density"]],
}

MIN_POINTS = {"linear": 2, "hermite": 3, "spline": 3, "lagrange": 4}


def interval_of(g: list[float], q: float) -> int:
    return max(j for j in range(len(g) - 1) if g[j] <= q)


def linear_weights(g: list[float], q: float) -> dict[int, float]:
    j = interval_of(g, q)
    t = (q - g[j]) / (g[j + 1] - g[j])
    return {j: 1 - t, j + 1: t}


def parabola_slope(g: list[float], k: int) -> dict[int, float]:
    """点 k の傾きの係数: k と両隣（端では端の 3 点）を通る放物線の微分（crates/core/src/state.rs の `slope`）"""
    s = min(max(k - 1, 0), len(g) - 3)
    p = g[s : s + 3]
    out = {}
    for m in range(3):
        a, b = [p[n] for n in range(3) if n != m]
        out[s + m] = (2 * g[k] - a - b) / ((p[m] - a) * (p[m] - b))
    return out


def hermite_weights(g: list[float], q: float) -> dict[int, float]:
    """3 次 Hermite の重み（crates/core/src/state.rs の `burnup_terms` の基底と傾き。燃焼度 0 の扱いは無い）"""
    j = interval_of(g, q)
    h = g[j + 1] - g[j]
    t = (q - g[j]) / h
    w = {j: (1 + 2 * t) * (1 - t) ** 2, j + 1: t * t * (3 - 2 * t)}
    for k, hk in ((j, t * (1 - t) ** 2), (j + 1, -t * t * (1 - t))):
        for p, c in parabola_slope(g, k).items():
            w[p] = w.get(p, 0.0) + h * hk * c
    return w


def lagrange_weights(g: list[float], q: float) -> dict[int, float]:
    j = interval_of(g, q)
    low = min(max(j - 1, 0), len(g) - 4)
    window = range(low, low + 4)
    return {
        m: math.prod((q - g[o]) / (g[m] - g[o]) for o in window if o != m)
        for m in window
    }


def spline_weights(g: list[float], q: float) -> dict[int, float]:
    """自然 3 次スプラインの重み: 点ごとに、その点だけ 1 の値のスプラインを q で引く"""
    n = len(g)
    a = np.zeros((n, n))
    a[0, 0] = a[-1, -1] = 1.0
    rhs = np.zeros((n, n))
    for i in range(1, n - 1):
        h0, h1 = g[i] - g[i - 1], g[i + 1] - g[i]
        a[i, i - 1 : i + 2] = [h0 / 6, (h0 + h1) / 3, h1 / 6]
        rhs[i, i - 1 : i + 2] = [1 / h0, -1 / h0 - 1 / h1, 1 / h1]
    # 2 階の微分（列 m は点 m だけ 1 の値のもの）
    second = np.linalg.solve(a, rhs)
    j = interval_of(g, q)
    h = g[j + 1] - g[j]
    left, right = (g[j + 1] - q) / h, (q - g[j]) / h
    out = {}
    for m in range(n):
        w = left * (m == j) + right * (m == j + 1)
        w += (
            ((left**3 - left) * second[j, m] + (right**3 - right) * second[j + 1, m])
            * h
            * h
            / 6
        )
        out[m] = float(w)
    return out


METHODS = {
    "linear": linear_weights,
    "hermite": hermite_weights,
    "lagrange": lagrange_weights,
    "spline": spline_weights,
}


def check_scheme(spec: dict, axes: dict) -> None:
    axis, variable, method = spec["axis"], spec["variable"], spec["method"]
    if axis not in axes:
        raise SystemExit(f"scheme: no state axis {axis!r} (axes {sorted(axes)})")
    if variable not in VARIABLES:
        raise SystemExit(
            f"scheme: variable {variable!r} is not one of {sorted(VARIABLES)}"
        )
    if variable == "density" and "density" not in axes[axis]:
        raise SystemExit(f"scheme: the {axis} axis has no densities")
    if variable == "log" and min(axes[axis]["points"]) <= 0:
        raise SystemExit(f"scheme: the {axis} axis has points at or below 0 for log")
    # 区間の選び方と重みは、変数が格子点の順に増えることを前提にする
    g = VARIABLES[variable](axes[axis])
    if any(b <= a for a, b in itertools.pairwise(g)):
        raise SystemExit(
            f"scheme: the {variable} of the {axis} axis does not increase strictly "
            "with the points (densities must fall with the temperature)"
        )
    if method not in METHODS:
        raise SystemExit(f"scheme: method {method!r} is not one of {sorted(METHODS)}")
    if len(axes[axis]["points"]) - 1 < MIN_POINTS[method]:
        raise SystemExit(
            f"scheme: {method} needs {MIN_POINTS[method]} points without the inner point, "
            f"the {axis} axis has {len(axes[axis]['points'])}"
        )


def summarize(table: dict, schemes: list[dict] = ()) -> dict:
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

    def errors(a: int, i: int, weights: dict[int, float]) -> list[dict]:
        """軸 a の格子点 i を、軸 a の格子点の重み `weights` で内挿した状態の、燃焼度の点ごとの |Δρ| の最大"""
        out = []
        for b in range(len(burnups)):
            worst = None
            for index in indices:
                if index[a] != i:
                    continue
                state, k_grid = truths[b][flat[index]]
                rebuilt = unvector(
                    sum(
                        c * vectors[b][flat[index[:a] + (k,) + index[a + 1 :]]]
                        for k, c in weights.items()
                    ),
                    state,
                )
                rho = (1 / k_grid - 1 / k_inf(rebuilt)) * 1e5
                if worst is None or abs(rho) > abs(worst[0]):
                    worst = (rho, index)
            out.append(
                {
                    "burnup": burnups[b],
                    "max_abs_pcm": abs(worst[0]),
                    "worst_pcm": worst[0],
                    "worst": {n: points[k][worst[1][k]] for k, n in enumerate(names)},
                }
            )
        return out

    centers = []
    gaps = {}
    for a, name in enumerate(names):
        g = [transforms[a](x) for x in points[a]]
        for i in range(1, len(g) - 1):
            w = (g[i] - g[i - 1]) / (g[i + 1] - g[i - 1])
            out = errors(a, i, {i - 1: 1 - w, i + 1: w})
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
    state_axes_by_name = {a["name"]: a for a in table["axes"] if a["name"] != "burnup"}
    scheme_rows = []
    for spec in schemes:
        check_scheme(spec, state_axes_by_name)
        a = names.index(spec["axis"])
        full = VARIABLES[spec["variable"]](state_axes_by_name[spec["axis"]])
        rows = []
        for i in range(1, len(full) - 1):
            rest = [k for k in range(len(full)) if k != i]
            weights = METHODS[spec["method"]]([full[k] for k in rest], full[i])
            rows.append(
                {
                    "point": points[a][i],
                    "points": errors(a, i, {rest[m]: c for m, c in weights.items()}),
                }
            )
        scheme_rows.append(
            {
                "axis": spec["axis"],
                "variable": spec["variable"],
                "method": spec["method"],
                "centers": rows,
            }
        )
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
        "schemes": scheme_rows,
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
    if summary["schemes"]:
        lines.append(
            "方式ごとに、格子点をその点を除いた格子から内挿した |Δρ| の最大（pcm。燃焼度の点ごと）:"
        )
    for x in summary["schemes"]:
        lines.append(f"  {x['axis']}、{x['variable']}、{x['method']}:")
        for c in x["centers"]:
            values = "、".join(
                f"{p['worst_pcm']:+.2f}（{p['worst']}）" for p in c["points"]
            )
            lines.append(f"    {c['point']:g}: {values}")
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 3 and args[0] == "--json":
        config_path = Path(args[2])
        config_bytes = config_path.read_bytes()
        config = tomllib.loads(config_bytes.decode())
        unknown = sorted(set(config) - {"table", "scheme"})
        if unknown:
            raise SystemExit(f"{config_path}: unknown keys {unknown}")
        path = config_path.parent / config["table"]
        data = path.read_bytes()
        for spec in config.get("scheme", []):
            extra = sorted(set(spec) - {"axis", "variable", "method"})
            if extra:
                raise SystemExit(f"{config_path}: unknown scheme keys {extra}")
        summary = summarize(msgpack.unpackb(data), config.get("scheme", []))
        summary["inputs"] = [
            {"path": str(config_path), "md5": hashlib.md5(config_bytes).hexdigest()},
            {"path": str(path), "md5": hashlib.md5(data).hexdigest()},
        ]
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
