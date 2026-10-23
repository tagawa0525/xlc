"""燃焼度の軸の内挿の方式を比べる（ADR-0192）。履歴の全点で作った分岐のテーブル（燃焼度のほかの軸は 1 点）から燃焼度の点を
選び（間引き）、選ばなかった点の k∞ を、方式ごとに選んだ点から内挿した少数群定数で解き、全点のテーブルのその点の状態を
真値として比べる。

使い方: python3 compare_burnup_interpolation.py --json <要約の出力先> <入力（TOML）>
        python3 compare_burnup_interpolation.py --report <要約の JSON>

入力の TOML（パスは TOML のあるディレクトリから）:

    [[table]]
    name = "2A"
    file = "table-2A-all.msgpack"
    selections = [{ name = "17 点", burnup = [0, 0.1, ...] }]

- k∞ は少数群定数の無限媒質の固有値（A φ = (1/k) χ νΣfᵀ φ、A = diag(Σt) − Σsᵀ。Σs は [元の群][先の群]）。集合体 1 体の core
  （反射境界、ノード 1 つ）の k と同じ量（examples/depletion/verification.md の「燃焼度の内挿の方式」）。全点のテーブルの各点で、
  この k∞ とテーブルの `k_inf`（格子計算の k∞）の差の最大を `table_k_max_abs_pcm` に書く
- 誤差は反応度の差 Δρ = (1/k_真値 − 1/k_内挿) × 1e5 pcm（compare_core_depletion.py の core − 履歴と同じ向き）
- 方式（`METHODS`）は、少数群定数の項目ごとに同じ方式で燃焼度について内挿する。`-with-zero` のない方式は、燃焼度 0 の点を
  傾き（曲線）に使わず、0 から次の点までを線形にする（ADR-0192。core の内挿と同じ）
- オーバーシュート: 選ばなかった点で、真値がその点を挟む 2 点の値の範囲 [lo, hi] の中にある項目について、内挿の値が範囲を
  越えた量 / max(|lo|, |hi|) の最大。真値も範囲の外にある項目（区間の中に極値がある）は数えない。chi は数えない（熱群の
  chi は 0 の丸めの値で、速群は 1 からその値を引いたもの）
- `inputs` は読んだファイルのパス（実行した cwd から）と md5
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib
from scipy.interpolate import Akima1DInterpolator, CubicSpline, PchipInterpolator

# 少数群定数のうち内挿して比べる項目（`k_inf` は格子計算の値で、内挿した定数からは解き直す）
ITEMS = [
    "total",
    "absorption",
    "nu_fission",
    "fission",
    "chi",
    "diffusion",
    "scatter",
    "adf",
]


def vector(state: dict) -> np.ndarray:
    """状態の項目を 1 列に並べる（形状関数があれば最後に足す）"""
    parts = [np.ravel(np.asarray(state[i], dtype=float)) for i in ITEMS]
    if state.get("pins"):
        parts.append(np.ravel(np.asarray(state["pins"]["form"], dtype=float)))
    return np.concatenate(parts)


def item_names(state: dict) -> list[str]:
    """`vector` の並びの項目名"""
    names = []
    for i in ITEMS:
        names += [i] * np.asarray(state[i]).size
    if state.get("pins"):
        names += ["form"] * np.asarray(state["pins"]["form"]).size
    return names


def unvector(v: np.ndarray, like: dict) -> dict:
    """`vector` の逆（k∞ を解くのに要る項目だけ）"""
    out = {}
    i = 0
    for name in ITEMS:
        shape = np.asarray(like[name]).shape
        n = int(np.prod(shape))
        out[name] = v[i : i + n].reshape(shape)
        i += n
    return out


def k_inf(state: dict) -> float:
    """少数群定数の無限媒質の k∞"""
    total = np.asarray(state["total"], dtype=float)
    scatter = np.asarray(state["scatter"], dtype=float)
    a = np.diag(total) - scatter.T
    f = np.outer(
        np.asarray(state["chi"], dtype=float),
        np.asarray(state["nu_fission"], dtype=float),
    )
    return float(max(np.linalg.eigvals(np.linalg.solve(a, f)).real))


def parabola_slopes(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """各点の傾き: その点と両隣（端では端の 3 点）を通る放物線の微分（core の `state::slope` と同じ）"""
    n = len(x)
    d = np.empty_like(y)
    for k in range(n):
        s = min(max(k - 1, 0), n - 3)
        p = x[s : s + 3]
        for m in range(3):
            others = [c for c in range(3) if c != m]
            numerator = sum(
                np.prod([x[k] - p[c] for c in others if c != o]) for o in others
            )
            denominator = np.prod([p[m] - p[c] for c in others])
            if m == 0:
                d[k] = y[s + m] * numerator / denominator
            else:
                d[k] += y[s + m] * numerator / denominator
    return d


def hermite(x: np.ndarray, y: np.ndarray, q: np.ndarray) -> np.ndarray:
    """3 点の放物線の傾きの 3 次 Hermite。点が 2 つなら線形"""
    if len(x) == 2:
        return linear(x, y, q)
    d = parabola_slopes(x, y)
    i = np.clip(np.searchsorted(x, q, side="right") - 1, 0, len(x) - 2)
    h = (x[i + 1] - x[i])[:, None]
    t = ((q - x[i]) / (x[i + 1] - x[i]))[:, None]
    h00 = (1 + 2 * t) * (1 - t) ** 2
    h01 = t * t * (3 - 2 * t)
    h10 = t * (1 - t) ** 2
    h11 = -t * t * (1 - t)
    return h00 * y[i] + h01 * y[i + 1] + h * (h10 * d[i] + h11 * d[i + 1])


def linear(x: np.ndarray, y: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.stack([np.interp(q, x, y[:, c]) for c in range(y.shape[1])], axis=1)


def lagrange3(x: np.ndarray, y: np.ndarray, q: np.ndarray) -> np.ndarray:
    """区間の 2 点と、q に近い側の隣の点の 3 点のラグランジュの 2 次（端の区間は内側の隣）"""
    out = np.empty((len(q), y.shape[1]))
    for n, v in enumerate(q):
        i = min(max(int(np.searchsorted(x, v, side="right")) - 1, 0), len(x) - 2)
        if i == 0:
            s = 0
        elif i == len(x) - 2:
            s = len(x) - 3
        else:
            s = i - 1 if v - x[i] < x[i + 1] - v else i
        p = x[s : s + 3]
        out[n] = sum(
            y[s + a] * np.prod([(v - p[b]) / (p[a] - p[b]) for b in range(3) if b != a])
            for a in range(3)
        )
    return out


def apart_from_zero(inner):
    """燃焼度 0 の点を曲線に使わず、0 から次の点までを線形にした方式（ADR-0192）"""

    def method(x: np.ndarray, y: np.ndarray, q: np.ndarray) -> np.ndarray:
        if x[0] != 0.0:
            return inner(x, y, q)
        out = np.empty((len(q), y.shape[1]))
        first = q <= x[1]
        out[first] = linear(x[:2], y[:2], q[first])
        if (~first).any():
            out[~first] = inner(x[1:], y[1:], q[~first])
        return out

    return method


# 方式の名前 → (記録の見出し, 内挿の関数)。名前は JSON のキーと記録の表の列
METHODS = {
    "linear": ("線形", linear),
    "lagrange3": ("3 点の 2 次", apart_from_zero(lagrange3)),
    "natural": (
        "自然 3 次スプライン",
        apart_from_zero(lambda x, y, q: CubicSpline(x, y, bc_type="natural")(q)),
    ),
    "pchip": ("PCHIP", apart_from_zero(lambda x, y, q: PchipInterpolator(x, y)(q))),
    "akima": ("Akima", apart_from_zero(lambda x, y, q: Akima1DInterpolator(x, y)(q))),
    "hermite": ("3 点の放物線の傾きの Hermite", apart_from_zero(hermite)),
    "akima-with-zero": (
        "Akima（0 も曲線に使う）",
        lambda x, y, q: Akima1DInterpolator(x, y)(q),
    ),
    "hermite-with-zero": ("3 点の放物線の傾きの Hermite（0 も曲線に使う）", hermite),
}


def interpolate(method: str, x: np.ndarray, y: np.ndarray, q: np.ndarray) -> np.ndarray:
    """燃焼度の点 x（昇順）の値 y（[点][項目]）を、q の燃焼度へ方式 `method` で内挿した値（[q][項目]）"""
    if method not in METHODS:
        raise SystemExit(f"unknown method {method}; choose from {list(METHODS)}")
    return np.asarray(METHODS[method][1](x, y, q), dtype=float)


def summarize_table(table: dict, selections: list[dict]) -> dict:
    """1 つのテーブルの、選び方ごと方式ごとの誤差（モジュールの説明）"""
    if table.get("format") != "xlc-table" or table.get("split") != 1:
        raise SystemExit(
            "the table must be an xlc-table that does not split the assembly"
        )
    axes = table["axes"]
    if axes[0]["name"] != "burnup":
        raise SystemExit("the first axis of the table must be the burnup")
    others = [a for a in axes[1:] if len(a["points"]) != 1]
    if others:
        raise SystemExit(
            f"the axes {[a['name'] for a in others]} have more than one point; give a table of one state per burnup"
        )
    points = [float(b) for b in axes[0]["points"]]
    states = table["states"]
    if len(states) != len(points):
        raise SystemExit(
            f"the table has {len(states)} states for {len(points)} burnups"
        )
    full = np.array([vector(s) for s in states])
    names = item_names(states[0])
    counted = np.array([n != "chi" for n in names])
    truth_k = np.array([k_inf(s) for s in states])
    table_k = np.array([s["k_inf"] for s in states])
    out = []
    for sel in selections:
        chosen = [float(b) for b in sel["burnup"]]
        if chosen != sorted(set(chosen)) or any(b not in points for b in chosen):
            raise SystemExit(
                f"{sel['name']}: the burnups {chosen} must be increasing points of the table"
            )
        if chosen[0] != points[0] or chosen[-1] != points[-1]:
            raise SystemExit(
                f"{sel['name']}: the burnups must run from the first to the last point of the table"
            )
        index = [points.index(b) for b in chosen]
        removed = [i for i in range(len(points)) if i not in index]
        x = np.array(chosen)
        q = np.array([points[i] for i in removed])
        results = {}
        for method in METHODS:
            values = interpolate(method, x, full[index], q)
            k = np.array([k_inf(unvector(v, states[0])) for v in values])
            rho = (1.0 / truth_k[removed] - 1.0 / k) * 1e5
            overshoot = 0.0
            for v, i in zip(values, removed, strict=True):
                hi_index = int(np.searchsorted(x, points[i]))
                a, b = full[index[hi_index - 1]], full[index[hi_index]]
                lo, hi = np.minimum(a, b), np.maximum(a, b)
                scale = np.maximum(np.abs(lo), np.abs(hi))
                inside = (full[i] >= lo) & (full[i] <= hi) & counted & (scale > 0)
                excess = np.maximum(lo - v, 0.0) + np.maximum(v - hi, 0.0)
                if inside.any():
                    overshoot = max(
                        overshoot, float(np.max(excess[inside] / scale[inside]))
                    )
            results[method] = {
                "rho_pcm": rho.tolist(),
                "max_abs_pcm": float(np.max(np.abs(rho))) if len(rho) else 0.0,
                "rms_pcm": float(np.sqrt(np.mean(rho**2))) if len(rho) else 0.0,
                "overshoot": overshoot,
            }
        out.append(
            {
                "name": sel["name"],
                "burnup": chosen,
                "removed": q.tolist(),
                "results": results,
            }
        )
    return {
        "burnup": points,
        "table_k_max_abs_pcm": float(np.max(np.abs(1.0 / table_k - 1.0 / truth_k)))
        * 1e5,
        "selections": out,
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    for t in summary["tables"]:
        lines.append(
            f"{t['name']}: 燃焼度の点 {len(t['burnup'])}、k∞ の解き直しとテーブルの k_inf の差の最大 {t['table_k_max_abs_pcm']:.4f} pcm"
        )
        for sel in t["selections"]:
            lines.append(
                f"{t['name']} {sel['name']}: 点 {len(sel['burnup'])}、比べた点 {len(sel['removed'])}"
            )
            for method, r in sel["results"].items():
                rho = r["rho_pcm"]
                if not rho:
                    lines.append(f"  {method}: 比べた点がない")
                    continue
                worst = int(np.argmax(np.abs(rho)))
                lines.append(
                    f"  {method}: |Δρ| の最大 {r['max_abs_pcm']:.1f} pcm（{sel['removed'][worst]:g} MWd/kg）、"
                    f"rms {r['rms_pcm']:.1f} pcm、{min(rho):+.1f}〜{max(rho):+.1f} pcm、オーバーシュート {r['overshoot']:.1e}"
                )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 3 and args[0] == "--json":
        config_path = Path(args[2])
        config = tomllib.loads(config_path.read_text())
        files = [config_path]
        tables = []
        for t in config["table"]:
            path = config_path.parent / t["file"]
            files.append(path)
            summary = summarize_table(
                msgpack.unpackb(path.read_bytes()), t["selections"]
            )
            tables.append({"name": t["name"], **summary})
        summary = {
            "format": "xlc-compare-burnup-interpolation",
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            "methods": {name: label for name, (label, _) in METHODS.items()},
            "tables": tables,
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
