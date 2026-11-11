"""直積の分岐のテーブルの燃料の温度の軸で、少数群定数と k∞ の温度への依存を物理のモデルの式と比べ、当てはめ（最小 2 乗）と
内挿の候補を同じ測り方で比べる（#0274 の「燃料の温度への依存のモデルと当てはめの式」）。

使い方: python3 fuel_models.py --json <要約の出力先> <直積のテーブル>
        python3 fuel_models.py --report <要約の JSON>

テーブルの軸は burnup、boron、moderator、fuel の順の直積（../compare_grid_interpolation.py と同じ並び）。燃焼度、ホウ素、減速材の
点の組（状態）ごとに燃料の温度の軸の点を見て、次を書く:

- `exponents`: 項目（`ITEMS`）と 1/k∞ を状態ごとに a + b T^p（p = 0 は a + b ln T）に最小 2 乗で当てはめ、残差の 2 乗和が最小の p
  （`EXPONENTS` の格子の上）の、燃焼度の点ごとの最小、中央値、最大（温度で変わらない状態は除いて `flat` に数え、すべて
  変わらなければ None）。p = 1/2 は共鳴積分の √T の相関、p = 1 は温度の線形
- `fits`: 基底（`BASES`。定数の項はどれにもある）ごとに、軸のすべての点に当てはめた残差。項目の残差は点の値の平均で割った値
  （1e-5）の、状態の中の最小と最大。k∞ は当てはめた少数群定数から解いた k∞ の Δρ。係数の滑らかさ（`roughness`）は、項目ごとに、
  基底の関数を公称の温度からの差にした係数を公称の温度の当てはめの値で割り、ホウ素と減速材の軸の線ごとにその軸の値の 2 次式に
  当てはめた残差の最大を、燃焼度の点の中の係数の幅で割った値（滑らかに動けば 0 に近く、状態ごとに跳べば 1 の桁。幅が係数の 1e-12 倍以下なら 0）
- `designs`: 格子点の選び方（`DESIGNS`）ごとに、選んだ点から方式（`SCHEMES`）で内挿か当てはめをし、選ばなかった点（選んだ点の
  範囲の中）の状態と比べた Δρ。`loo` は内側の点を 1 つずつ除く（../compare_grid_interpolation.py の `schemes` と同じ測り方）。
  方式の要る点の数に足りない選び方は `skipped` に書く。選び方は添字で決めるので、燃料の温度の軸は 7 点（`FUEL_POINTS`）
- `surrogates`（見込み。測定でない）: 状態ごとに軸のすべての点に当てはめた滑らかな曲線（`SURROGATES`）を真値とみなし、その曲線の
  格子点の値から、すべての点か選び方（`DESIGNS` の `loo` はすべての点）の点で方式で引いた区間の中点（温度の中点）の状態と、
  曲線の中点の状態の Δρ。曲線と同じ基底の当てはめは曲線そのものなので `skipped` に書く。格子点の値の当てはめの残差（格子点ごとの
  外れ）を除いた、方式の形の誤差を見る

Δρ = (1/k_格子点 − 1/k_候補) × 1e5 pcm（k∞ は少数群定数の無限媒質の固有値。../compare_burnup_interpolation.py の `k_inf`）で、
燃焼度の点ごと（`burnup` の順）に、ホウ素と減速材のすべての組の中で |Δρ| が最大の値（`worst_pcm`）とその状態（`worst`。
`surrogates` は値だけ）を書く。方式の値はどれも格子点の値の線形結合で、少数群定数の
項目ごとに同じ重みで足す（当てはめの重みは最小 2 乗の射影の行）。`inputs` は読んだテーブルのパスと md5。
"""

import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from compare_branch_cross_terms import table_grid
from compare_burnup_interpolation import item_names, k_inf, unvector, vector
from compare_grid_interpolation import METHODS, MIN_POINTS

AXES = ["burnup", "boron", "moderator", "fuel"]

# 指数と残差を見る項目（名前、少数群定数の `vector` の中の項目と添字）。散乱は [元の群][先の群] で、1→2 は減速の除去
ITEMS = [
    ("高速群の吸収", "absorption", (0,)),
    ("熱群の吸収", "absorption", (1,)),
    ("高速群の νΣf", "nu_fission", (0,)),
    ("熱群の νΣf", "nu_fission", (1,)),
    ("1→2 の散乱", "scatter", (0, 1)),
    ("高速群の拡散係数", "diffusion", (0,)),
    ("熱群の拡散係数", "diffusion", (1,)),
]

EXPONENTS = [round(-1.5 + 0.005 * n, 3) for n in range(701)]

# 基底の関数（定数の項は足す）。物理のモデルからの意味は #0274 の「燃料の温度への依存のモデルと当てはめの式」
BASES = {
    "√T": [math.sqrt],
    "√T、T": [math.sqrt, lambda t: t],
    "√T、ln T": [math.sqrt, math.log],
    "√T、T、1/√T": [math.sqrt, lambda t: t, lambda t: t**-0.5],
    "√T、T、ln T": [math.sqrt, lambda t: t, math.log],
    "√T、T、T^1.5": [math.sqrt, lambda t: t, lambda t: t**1.5],
}

VARIABLES = {"√T": math.sqrt, "T": lambda t: t, "ln T": math.log}

# (種類, 変数か基底, 方式)。内挿は ../compare_grid_interpolation.py の `METHODS` と同じ重み
SCHEMES = [("interpolate", v, m) for v in VARIABLES for m in METHODS] + [
    ("fit", b, "最小 2 乗") for b in BASES
]

# 本番の間隔の区間の中点の見込みに使う滑らかな曲線（`BASES` の名前）。√T、T、T^1.5 は √T の 3 次式で、√T の 4 点のラグランジュが
# 厳密に再現して形の誤差が 0 になるので使わない。どれも内挿の変数（`VARIABLES`）の多項式でない
SURROGATES = ["√T、T、1/√T", "√T、T、ln T"]

# 格子点の選び方（燃料の温度の軸の添字。軸は 7 点で、端の 2 点をどれも含む）。None は内側の点を 1 つずつ除く
FUEL_POINTS = 7
DESIGNS = {"loo": None, "4 点": [0, 2, 4, 6], "3 点": [0, 3, 6]}


def design_matrix(
    basis: list, xs: list[float], origin: float | None = None
) -> np.ndarray:
    """行が点、列が定数の項と基底の関数（`origin` があれば、そこでの値を引く）の行列"""
    shift = [0.0] * len(basis) if origin is None else [f(origin) for f in basis]
    return np.array(
        [[1.0] + [f(x) - c for f, c in zip(basis, shift, strict=True)] for x in xs]
    )


def projection(basis: list, xs: list[float]) -> np.ndarray:
    """xs の値を最小 2 乗の当てはめの xs の値に移す行列"""
    a = design_matrix(basis, xs)
    return a @ np.linalg.pinv(a)


def fit_weights(basis: list, xs: list[float], q: float) -> dict[int, float]:
    """xs の値への最小 2 乗の当てはめを q で引く重み（射影の行）"""
    if len(xs) < len(basis) + 1:
        raise ValueError(f"{len(basis) + 1} terms need as many points, not {len(xs)}")
    row = design_matrix(basis, [q]) @ np.linalg.pinv(design_matrix(basis, xs))
    return dict(enumerate(row[0].tolist()))


def scheme_points(scheme: tuple) -> int:
    kind, name, method = scheme
    return len(BASES[name]) + 1 if kind == "fit" else MIN_POINTS[method]


def scheme_weights(scheme: tuple, xs: list[float], q: float) -> dict[int, float]:
    kind, name, method = scheme
    if kind == "fit":
        return fit_weights(BASES[name], xs, q)
    g = VARIABLES[name]
    return METHODS[method]([g(x) for x in xs], g(q))


def load(table: dict):
    names = [a["name"] for a in table["axes"]]
    if names != AXES:
        raise SystemExit(f"the table needs the axes {AXES}, not {names}")
    if any("nominal" not in a for a in table["axes"][1:]):
        raise SystemExit("the table needs the nominal value of every state axis")
    # 直積の順、集合体を分けないこと、状態の数（compare_grid_interpolation.py と同じ検査）
    table_grid("table", table, {a["name"]: a["nominal"] for a in table["axes"][1:]})
    points = [list(a["points"]) for a in table["axes"]]
    shape = [len(p) for p in points]
    if len(points[3]) != FUEL_POINTS:
        raise SystemExit(
            f"the designs pick {FUEL_POINTS} fuel points, the table has {len(points[3])}"
        )
    # 係数の粗さは、ホウ素と減速材の線ごとの 2 次式からの残差で、3 点以下では 0 になる
    if min(len(points[1]), len(points[2])) < 4:
        raise SystemExit("the roughness needs 4 boron and 4 moderator points")
    if min(points[3]) <= 0:
        raise SystemExit("the fuel temperatures must be above 0 K")
    if any(b <= a for a, b in itertools.pairwise(points[3])):
        raise SystemExit("the fuel temperatures must increase")
    vectors = np.array([vector(s) for s in table["states"]]).reshape(*shape, -1)
    return points, vectors, table["states"][0]


def item_column(like: dict, item: str, index: tuple) -> int:
    names = item_names(like)
    start = names.index(item)
    return start + int(np.ravel_multi_index(index, np.asarray(like[item]).shape))


def best_exponents(t: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """ys の行ごとに、a + b T^p の残差の 2 乗和が最小の p（`EXPONENTS` の格子の上）。温度で変わらない行（変化が値の 1e-12 倍
    以下）は NaN。最小が格子の端なら、べきが格子の外にあるのでエラー"""
    yc = ys - ys.mean(axis=1, keepdims=True)
    errors = []
    for p in EXPONENTS:
        x = np.log(t) if p == 0 else t**p
        xc = x - x.mean()
        # 1 次の回帰の残差の 2 乗和 = Σ yc² − (Σ yc xc)² / Σ xc²
        errors.append((yc**2).sum(axis=1) - (yc @ xc) ** 2 / (xc @ xc))
    best = np.argmin(np.array(errors), axis=0)
    flat = np.abs(yc).max(axis=1) <= 1e-12 * np.abs(ys).max(axis=1)
    if np.any(~flat & ((best == 0) | (best == len(EXPONENTS) - 1))):
        raise SystemExit(
            f"an exponent lies outside [{EXPONENTS[0]}, {EXPONENTS[-1]}]; widen EXPONENTS"
        )
    return np.where(flat, np.nan, np.array(EXPONENTS)[best])


def quadratic_roughness(x: np.ndarray, y: np.ndarray) -> float:
    """y を x の 2 次式に当てはめた残差の絶対値の最大"""
    c = np.polyfit(x, y, 2)
    return float(np.max(np.abs(y - np.polyval(c, x))))


def summarize(table: dict) -> dict:
    points, vectors, like = load(table)
    burnups, borons, moderators, fuels = points
    t = np.array(fuels)
    # load が状態の軸のすべてに公称の値があることを確かめる
    nominal = table["axes"][3]["nominal"]
    rest = list(itertools.product(range(len(borons)), range(len(moderators))))

    def rebuilt_k(u, b, m, weights):
        v = sum(c * vectors[u, b, m, k] for k, c in weights.items())
        return k_inf(unvector(v, like))

    def truth_of(weights):
        """燃料の温度の添字の重みで作った状態の 1e5/k∞（燃焼度、ホウ素、減速材）"""
        out = np.empty(inverse_k.shape[:3])
        for u in range(len(burnups)):
            for b, m in rest:
                out[u, b, m] = 1e5 / rebuilt_k(u, b, m, weights)
        return out

    def errors(truth, weights, state=True):
        """1e5/k∞ の真値 `truth`（燃焼度、ホウ素、減速材）と、重み `weights`（燃料の温度の添字）で作った状態の Δρ の、
        燃焼度の点ごとの |Δρ| の最大の値（`state` なら、その状態も）"""
        out = []
        for u in range(len(burnups)):
            worst = None
            for b, m in rest:
                rho = truth[u, b, m] - 1e5 / rebuilt_k(u, b, m, weights)
                if worst is None or abs(rho) > abs(worst[0]):
                    worst = (rho, b, m)
            entry = {"worst_pcm": worst[0]}
            if state:
                entry["worst"] = {
                    "boron": borons[worst[1]],
                    "moderator": moderators[worst[2]],
                }
            out.append(entry)
        return out

    columns = [(label, item_column(like, item, index)) for label, item, index in ITEMS]
    inverse_k = 1e5 / np.array(
        [
            [
                [
                    [
                        k_inf(unvector(vectors[u, b, m, f], like))
                        for f in range(len(fuels))
                    ]
                    for m in range(len(moderators))
                ]
                for b in range(len(borons))
            ]
            for u in range(len(burnups))
        ]
    )
    exponents = []
    for label, column in [*columns, ("1/k∞", None)]:
        rows = []
        for u in range(len(burnups)):
            ys = [
                vectors[u, b, m, :, column]
                if column is not None
                else inverse_k[u, b, m]
                for b, m in rest
            ]
            ys = np.array(ys)
            ps = best_exponents(t, ys)
            change = (ys[:, -1] - ys[:, 0]) / ys[:, 0]
            changing = ps[~np.isnan(ps)]
            rows.append(
                {
                    "burnup": burnups[u],
                    # 温度で変わらない状態の数（べきを持たない）
                    "flat": int(np.isnan(ps).sum()),
                    "min": float(changing.min()) if changing.size else None,
                    "median": float(np.median(changing)) if changing.size else None,
                    "max": float(changing.max()) if changing.size else None,
                    "relative_change_median": float(np.median(change)),
                }
            )
        exponents.append({"item": label, "burnups": rows})

    fits = []
    for name, basis in BASES.items():
        onto = projection(basis, fuels)
        centered_inverse = np.linalg.pinv(design_matrix(basis, fuels, nominal))
        items = []
        for label, column in columns:
            # (燃焼度, ホウ素, 減速材, 燃料の温度)
            y = vectors[..., column]
            r = (y - y @ onto.T) / y.mean(axis=-1, keepdims=True) * 1e5
            c = y @ centered_inverse.T
            coefficients = c[..., 1:] / c[..., :1]
            r = r.reshape(-1, len(fuels))
            roughness = []
            for u, c in enumerate(coefficients):
                row = []
                for j in range(c.shape[2]):
                    cj = c[:, :, j]
                    width = float(cj.max() - cj.min())
                    lines = [
                        quadratic_roughness(np.array(moderators), cj[b])
                        for b in range(len(borons))
                    ]
                    lines += [
                        quadratic_roughness(np.array(borons), cj[:, m])
                        for m in range(len(moderators))
                    ]
                    row.append(
                        {
                            "min": float(cj.min()),
                            "max": float(cj.max()),
                            # 状態によらず同じ係数（幅が丸めの桁）は跳ばない
                            "roughness": max(lines) / width
                            if width > 1e-12 * float(np.abs(cj).max())
                            else 0.0,
                        }
                    )
                roughness.append({"burnup": burnups[u], "coefficients": row})
            items.append(
                {
                    "item": label,
                    "residual_min": r.min(0).tolist(),
                    "residual_max": r.max(0).tolist(),
                    "coefficients": roughness,
                }
            )
        k_rows = [
            {
                "point": fuels[i],
                "points": errors(inverse_k[..., i], dict(enumerate(onto[i].tolist()))),
            }
            for i in range(len(fuels))
        ]
        fits.append({"basis": name, "items": items, "k": k_rows})

    designs = []
    for design, chosen in DESIGNS.items():
        rows = []
        skipped = []
        cases = (
            [
                (i, [k for k in range(len(fuels)) if k != i])
                for i in range(1, len(fuels) - 1)
            ]
            if chosen is None
            else [
                (i, chosen)
                for i in range(len(fuels))
                if i not in chosen and chosen[0] < i < chosen[-1]
            ]
        )
        for scheme in SCHEMES:
            if any(len(c) < scheme_points(scheme) for _, c in cases):
                skipped.append(list(scheme))
                continue
            centers = []
            for i, use in cases:
                weights = scheme_weights(scheme, [fuels[k] for k in use], fuels[i])
                centers.append(
                    {
                        "point": fuels[i],
                        "points": errors(
                            inverse_k[..., i], {use[m]: c for m, c in weights.items()}
                        ),
                    }
                )
            rows.append({"scheme": list(scheme), "centers": centers})
        designs.append(
            {
                "design": design,
                "points": None if chosen is None else [fuels[k] for k in chosen],
                "schemes": rows,
                "skipped": skipped,
            }
        )

    surrogates = []
    for name in SURROGATES:
        onto = projection(BASES[name], fuels)
        for design, chosen in DESIGNS.items():
            use = list(range(len(fuels))) if chosen is None else chosen
            # 区間の中点と、曲線の中点の状態の 1e5/k∞（方式によらない）
            midpoints = [(fuels[a] + fuels[b]) / 2 for a, b in itertools.pairwise(use)]
            truths = [truth_of(fit_weights(BASES[name], fuels, q)) for q in midpoints]
            rows = []
            skipped = []
            for scheme in SCHEMES:
                # 曲線と同じ基底の当てはめは、曲線そのものを返す
                if len(use) < scheme_points(scheme) or scheme == (
                    "fit",
                    name,
                    "最小 2 乗",
                ):
                    skipped.append(list(scheme))
                    continue
                centers = []
                for q, truth in zip(midpoints, truths, strict=True):
                    # 方式は、選んだ格子点の、滑らかな曲線の値（射影）から引く
                    dense = np.zeros(len(fuels))
                    weights = scheme_weights(scheme, [fuels[k] for k in use], q)
                    for m, c in weights.items():
                        dense[use[m]] += c
                    row = dense @ onto
                    centers.append(
                        {
                            "point": q,
                            "points": errors(
                                truth, dict(enumerate(row.tolist())), state=False
                            ),
                        }
                    )
                rows.append({"scheme": list(scheme), "centers": centers})
            surrogates.append(
                {
                    "basis": name,
                    "design": design,
                    "points": [fuels[k] for k in use],
                    "schemes": rows,
                    "skipped": skipped,
                }
            )

    return {
        "format": "xlc-fuel-temperature-models",
        "burnup": burnups,
        "fuel": fuels,
        "nominal": nominal,
        "exponents": exponents,
        "fits": fits,
        "designs": designs,
        "surrogates": surrogates,
    }


def scheme_label(scheme: list) -> str:
    kind, name, method = scheme
    return f"当てはめ {name}" if kind == "fit" else f"内挿 {name} {method}"


def scheme_lines(group: dict, skipped: str) -> list[str]:
    """方式ごとに、|Δρ| の最大と、点ごとに燃焼度の点の中で |Δρ| が最大の値"""
    lines = []
    for s in group["schemes"]:
        cells = "、".join(
            f"{c['point']:g}: "
            + f"{max(c['points'], key=lambda p: abs(p['worst_pcm']))['worst_pcm']:+.2f}"
            for c in s["centers"]
        )
        worst = max(abs(p["worst_pcm"]) for c in s["centers"] for p in c["points"])
        lines.append(f"  {scheme_label(s['scheme'])}: 最大 {worst:.2f}（{cells}）")
    for s in group["skipped"]:
        lines.append(f"  {scheme_label(s)}: {skipped}")
    return lines


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    fuels = summary["fuel"]
    burnups = summary["burnup"]
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        f"a + b T^p の p（燃焼度の点ごとに 最小/中央値/最大、{fuels[0]:g}→{fuels[-1]:g} K の相対変化の中央値）:"
    )
    for e in summary["exponents"]:
        cells = "、".join(
            f"{r['burnup']:g}: "
            + (
                f"{r['min']:.3f}/{r['median']:.3f}/{r['max']:.3f}"
                if r["min"] is not None
                else "なし"
            )
            + (f"（温度で変わらない状態 {r['flat']}）" if r["flat"] else "")
            + f"（{r['relative_change_median']:+.2e}）"
            for r in e["burnups"]
        )
        lines.append(f"  {e['item']}: {cells}")
    for f in summary["fits"]:
        lines.append(f"当てはめ {f['basis']}（軸のすべての点）:")
        for i in f["items"]:
            cells = " ".join(
                f"{lo:+.1f}〜{hi:+.1f}"
                for lo, hi in zip(i["residual_min"], i["residual_max"], strict=True)
            )
            lines.append(f"  {i['item']} の残差（1e-5）: {cells}")
            for c in i["coefficients"]:
                coefficients = "、".join(
                    f"{x['min']:+.3e}〜{x['max']:+.3e}（{x['roughness']:.3f}）"
                    for x in c["coefficients"]
                )
                lines.append(
                    f"    {c['burnup']:g} MWd/kg の係数（粗さ）: {coefficients}"
                )
        cells = " ".join(
            f"{max(abs(p['worst_pcm']) for p in r['points']):.2f}" for r in f["k"]
        )
        lines.append(
            f"  Δρ の最大（pcm、{'、'.join(f'{x:g}' for x in fuels)} K）: {cells}"
        )
        for r in f["k"]:
            values = "、".join(f"{p['worst_pcm']:+.2f}" for p in r["points"])
            lines.append(
                f"    {r['point']:g} K（燃焼度 {'、'.join(f'{b:g}' for b in burnups)}）: {values}"
            )
    for d in summary["designs"]:
        chosen = (
            "内側の点を 1 つずつ除く"
            if d["points"] is None
            else "、".join(f"{x:g}" for x in d["points"])
        )
        lines.append(
            f"選び方 {d['design']}（{chosen}）の Δρ（pcm。点ごとに燃焼度の点の中で |Δρ| が最大の値）:"
        )
        lines += scheme_lines(d, "点が足りない")
    for g in summary["surrogates"]:
        chosen = "、".join(f"{x:g}" for x in g["points"])
        lines.append(
            f"当てはめ {g['basis']} の曲線を真値とみなし、{chosen} K から引いた区間の中点の Δρ（見込み、pcm）:"
        )
        lines += scheme_lines(g, "点が足りないか、曲線と同じ基底")
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
            json.dumps(summary, ensure_ascii=False, separators=(",", ":")) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
