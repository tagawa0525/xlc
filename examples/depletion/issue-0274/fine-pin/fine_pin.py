"""燃料の温度を細かく取ったピンセルのテーブルで、物理の基底の当てはめの残差が、テーブルの格子点だけの片寄りか、滑らかな曲がりかを
見る（#0274 の「燃料の温度を細かく取ったピンセル」）。

使い方: python3 fine_pin.py --json <要約の出力先> --table <直積のテーブル> [--reference <名前>=<7 点のピンセル>]
            [--library <ライブラリ> <比べるライブラリ>] <名前>=<ピンセルのテーブル>:<燃焼度>...
        python3 fine_pin.py --report <要約の JSON>

ピンセルのテーブルは燃料の温度の軸だけが複数の点を持ち（ほかの軸は 1 点）、`GRID`（直積のテーブルの燃料の温度の 7 点）を
含むこと。`GRID` の点を「格子点」、ほかを「足した点」と呼ぶ。テーブルごとに次を書く:

- `fits`: 基底（../fuel_models.py の `BASES`）ごとに、すべての点か格子点だけに最小 2 乗で当てはめた残差。項目（`ITEMS`）は点の値の
  平均で割った値（1e-5）、k∞ は Δρ（pcm）。すべての点の当てはめでは、内側の点の残差から、両隣の点の残差を温度で線形に内挿した値を
  引いた「跳び」も書く（残差が滑らかな曲がりなら、跳びは残差より桁で小さく、格子点の片寄りなら格子点の跳びが残差と同じ桁になる）
- `predictions`: 格子点だけから方式（../fuel_models.py の `SCHEMES`）で引いた足した点の Δρ（pcm）。本番の間隔の、方式の実測の誤差
- `product`: 直積のテーブルの公称の状態（ホウ素と減速材の軸の `nominal`）の、その燃焼度の点の、格子点の当てはめの残差（ピンセルの
  格子点の当てはめと形を比べる）

`--reference` は、名前のピンセルと同じピンセルを格子点だけで、ほかのライブラリで解いたテーブルで、格子点の状態の少数群定数の
最大の相対差を `reference` に書く。`--library` は 2 つのライブラリの核種と温度の組のうち、両方にあるもののデータが同じバイトか（`library`。
同じでない組を並べる）。Δρ = (1/k_点 − 1/k_候補) × 1e5、k∞ は少数群定数の無限媒質の固有値（../../compare_burnup_interpolation.py の
`k_inf`）。`inputs` は読んだファイルのパスと md5。
"""

import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path

import msgpack
import numpy as np

# ../fuel_models.py と ../../compare_burnup_interpolation.py
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from compare_burnup_interpolation import k_inf, unvector, vector
from fuel_models import BASES, ITEMS, SCHEMES, item_column, projection, scheme_weights
from fuel_models import load as load_product

GRID = [560.0, 600.0, 750.0, 900.0, 1050.0, 1200.0, 1350.0]

# 跳びを書く項目（`ITEMS` の名前）と k∞
JUMP_ITEMS = ["高速群の吸収", "熱群の吸収", "高速群の νΣf"]


def load_pin(table: dict):
    """燃料の温度の点、点ごとの少数群定数のベクトル、形の見本の状態"""
    axes = {a["name"]: a["points"] for a in table["axes"]}
    if set(axes) - {"boron", "moderator", "fuel"}:
        raise SystemExit(f"the pin table has unexpected axes {sorted(axes)}")
    if any(len(p) != 1 for name, p in axes.items() if name != "fuel"):
        raise SystemExit("every axis of the pin table but fuel must hold one point")
    fuels = [float(t) for t in axes["fuel"]]
    if len(table["states"]) != len(fuels):
        raise SystemExit("the pin table must hold one state per fuel point")
    if any(b <= a for a, b in itertools.pairwise(fuels)):
        raise SystemExit("the fuel temperatures must increase")
    if not set(GRID) <= set(fuels):
        raise SystemExit(f"the pin table needs the grid points {GRID} K")
    return fuels, np.array([vector(s) for s in table["states"]]), table["states"][0]


def inverse_k(v: np.ndarray, like: dict) -> float:
    return 1e5 / k_inf(unvector(v, like))


def residuals(
    basis: list, xs: list[float], vectors: np.ndarray, like: dict, columns
) -> dict:
    """xs の点の vectors を基底 basis に当てはめた残差（項目は平均で割った 1e-5、k∞ は Δρ の pcm）"""
    fitted = projection(basis, xs) @ vectors
    out = {
        label: ((vectors[:, c] - fitted[:, c]) / vectors[:, c].mean() * 1e5).tolist()
        for label, c in columns
    }
    out["k"] = [
        inverse_k(v, like) - inverse_k(f, like)
        for v, f in zip(vectors, fitted, strict=True)
    ]
    return out


def jumps(xs: list[float], r: list[float]) -> list[float]:
    """内側の点ごとに、残差から両隣の残差を xs で線形に内挿した値を引いた値"""
    out = []
    for i in range(1, len(xs) - 1):
        w = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1])
        out.append(r[i] - ((1 - w) * r[i - 1] + w * r[i + 1]))
    return out


def summarize_pin(table: dict) -> dict:
    fuels, vectors, like = load_pin(table)
    columns = [(label, item_column(like, item, index)) for label, item, index in ITEMS]
    grid = [fuels.index(t) for t in GRID]
    added = [i for i in range(len(fuels)) if i not in grid]
    fits = []
    for name, basis in BASES.items():
        whole = residuals(basis, fuels, vectors, like, columns)
        fits.append(
            {
                "basis": name,
                "all": whole,
                "jumps": {
                    label: jumps(fuels, whole[label]) for label in [*JUMP_ITEMS, "k"]
                },
                "grid": residuals(basis, GRID, vectors[grid], like, columns),
            }
        )
    predictions = []
    for scheme in SCHEMES:
        rows = []
        for i in added:
            weights = scheme_weights(scheme, GRID, fuels[i])
            v = sum(c * vectors[grid[m]] for m, c in weights.items())
            rows.append(inverse_k(vectors[i], like) - inverse_k(v, like))
        predictions.append({"scheme": list(scheme), "pcm": rows})
    return {
        "fuel": fuels,
        "grid": GRID,
        "added": [fuels[i] for i in added],
        "k_inf": [k_inf(unvector(v, like)) for v in vectors],
        "fits": fits,
        "predictions": predictions,
    }


def summarize_product(table: dict, burnup: float) -> dict:
    """直積のテーブルの公称の状態の、燃焼度 burnup の点の、格子点の当てはめの残差"""
    points, vectors, like = load_product(table)
    burnups, borons, moderators, fuels = points
    if fuels != GRID:
        raise SystemExit(f"the product table needs the fuel points {GRID} K")
    if burnup not in burnups:
        raise SystemExit(f"the product table has no burnup point {burnup}")
    nominal = {a["name"]: a["nominal"] for a in table["axes"][1:]}
    u = burnups.index(burnup)
    b = borons.index(nominal["boron"])
    m = moderators.index(nominal["moderator"])
    columns = [(label, item_column(like, item, index)) for label, item, index in ITEMS]
    return {
        "burnup": burnup,
        "boron": nominal["boron"],
        "moderator": nominal["moderator"],
        "fits": [
            {
                "basis": name,
                "grid": residuals(basis, GRID, vectors[u, b, m], like, columns),
            }
            for name, basis in BASES.items()
        ],
    }


def reference_difference(table: dict, reference: dict) -> float:
    """reference の燃料の温度の点の状態と、table の同じ点の状態の少数群定数の相対差の最大"""
    fuels, vectors, _ = load_pin(table)
    axes = {a["name"]: a["points"] for a in reference["axes"]}
    if len(reference["states"]) != len(axes["fuel"]):
        raise SystemExit("the reference must hold one state per fuel point")
    worst = 0.0
    for t, s in zip(axes["fuel"], reference["states"], strict=True):
        if float(t) not in fuels:
            raise SystemExit(f"the reference point {t} K is not in the pin table")
        a, r = vectors[fuels.index(float(t))], vector(s)
        scale = np.where(np.abs(r) > 0, np.abs(r), 1.0)
        worst = max(worst, float(np.max(np.abs(a - r) / scale)))
    return worst


def library_difference(library: dict, other: dict) -> dict:
    """両方のライブラリにある核種と温度の組の数と、データが同じバイトでない組"""
    theirs = {
        (n["name"], t["temperature"]): msgpack.packb(t)
        for n in other["nuclides"]
        for t in n["temperatures"]
    }
    shared, differ = 0, []
    for n in library["nuclides"]:
        for t in n["temperatures"]:
            key = (n["name"], t["temperature"])
            if key in theirs:
                shared += 1
                if msgpack.packb(t) != theirs[key]:
                    differ.append(list(key))
    return {"shared": shared, "differ": differ}


def read(path: str, inputs: list) -> dict:
    data = Path(path).read_bytes()
    inputs.append({"path": path, "md5": hashlib.md5(data).hexdigest()})
    return msgpack.unpackb(data, strict_map_key=False)


def summarize(args) -> dict:
    inputs = []
    product = read(args.table, inputs)
    pins, tables = [], {}
    for spec in args.pins:
        name, rest = spec.split("=", 1)
        path, burnup = rest.rsplit(":", 1)
        tables[name] = read(path, inputs)
        entry = {"name": name, "burnup": float(burnup), **summarize_pin(tables[name])}
        entry["product"] = summarize_product(product, float(burnup))
        pins.append(entry)
    summary = {"format": "xlc-fine-pin", "pins": pins}
    if args.reference:
        name, path = args.reference.split("=", 1)
        if name not in tables:
            raise SystemExit(f"--reference names {name}, not one of {sorted(tables)}")
        summary["reference"] = {
            name: reference_difference(tables[name], read(path, inputs))
        }
    if args.library:
        mine, theirs = (read(p, inputs) for p in args.library)
        summary["library"] = library_difference(mine, theirs)
    summary["inputs"] = inputs
    return summary


def scheme_label(scheme: list) -> str:
    kind, name, method = scheme
    return f"当てはめ {name}" if kind == "fit" else f"内挿 {name} {method}"


def cells(values: list[float], digits: int = 2) -> str:
    return " ".join(f"{x:+.{digits}f}" for x in values)


def report(summary: dict) -> list[str]:
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    if "library" in summary:
        lib = summary["library"]
        lines.append(
            f"ライブラリ: 両方にある核種と温度の組 {lib['shared']}、同じバイトでない組 {len(lib['differ'])}"
            + (f"（{lib['differ']}）" if lib["differ"] else "")
        )
    for name, diff in summary.get("reference", {}).items():
        lines.append(f"{name}: 7 点のピンセルとの少数群定数の相対差の最大 {diff:.3e}")
    for p in summary["pins"]:
        fuels = p["fuel"]
        lines.append(
            f"== {p['name']}（燃焼度 {p['burnup']:g}）: 燃料の温度 {'、'.join(f'{t:g}' for t in fuels)} K"
        )
        lines.append("  k∞: " + " ".join(f"{k:.6f}" for k in p["k_inf"]))
        marks = " ".join("*" if t in p["grid"] else "." for t in fuels)
        lines.append(f"  格子点（*）: {marks}")
        for f in p["fits"]:
            lines.append(f"  基底 {f['basis']}、すべての点の残差:")
            for label in [*JUMP_ITEMS, "熱群の νΣf", "k"]:
                lines.append(f"    {label}: {cells(f['all'][label])}")
            for label, values in f["jumps"].items():
                inner = fuels[1:-1]
                on = [
                    abs(v) for t, v in zip(inner, values, strict=True) if t in p["grid"]
                ]
                off = [
                    abs(v)
                    for t, v in zip(inner, values, strict=True)
                    if t not in p["grid"]
                ]
                lines.append(
                    f"    跳び {label}: {cells(values)}（|跳び| の最大、格子点 {max(on):.2f} / 足した点 {max(off):.2f}）"
                )
            lines.append(f"  基底 {f['basis']}、格子点だけの残差:")
            product = next(x for x in p["product"]["fits"] if x["basis"] == f["basis"])
            for label in [*JUMP_ITEMS, "k"]:
                lines.append(
                    f"    {label}: ピンセル {cells(f['grid'][label])} / 直積の公称 {cells(product['grid'][label])}"
                )
        lines.append(
            f"  格子点から引いた足した点（{'、'.join(f'{t:g}' for t in p['added'])} K）の Δρ（pcm）:"
        )
        for s in p["predictions"]:
            worst = max(abs(x) for x in s["pcm"])
            lines.append(
                f"    {scheme_label(s['scheme'])}: 最大 {worst:.2f}（{cells(s['pcm'])}）"
            )
    return lines


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--report":
        print("\n".join(report(json.loads(Path(sys.argv[2]).read_text()))))
        return
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--json", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--reference")
    parser.add_argument("--library", nargs=2)
    parser.add_argument("pins", nargs="+")
    args = parser.parse_args()
    summary = summarize(args)
    Path(args.json).write_text(
        json.dumps(summary, ensure_ascii=False, separators=(",", ":")) + "\n"
    )
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
