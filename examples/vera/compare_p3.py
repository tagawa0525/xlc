"""VERA Problem 3A、3B の xlc の炉心の解を、仕様書の KENO-VI の値と、OpenMC の参照解（あれば）と比べる。

使い方: python3 compare_p3.py <問題（3A、3B）> <spec.txt> <core の標準出力> <core の pins.msgpack> <core の layers.msgpack>
        [<OpenMC の結果のディレクトリ>]

- k: core の標準出力の `k-eff` と、表 P3-4 の KENO-VI の k（差は (k_xlc − k_ref) × 1e5 pcm）、OpenMC の k.json
- 軸方向の出力: 付録 D の 49 層の KENO-VI の値（燃料の高さの平均で 1 の線密度）と、OpenMC の layers.msgpack の `axial`。
  xlc の層（`xlc-axial-pin-powers`）は参照の層を細かく分けたものでよく、参照の層ごとに厚さで重みを付けて平均する。
  層の境界は、どちらも燃料の下端からの高さで合わせる（xlc は炉心の下端から、参照は下部炉心板の上面から）。
  差は層ごとの (xlc − 参照)（出力の単位のまま。軸方向の出力は 1 程度）と、軸方向のオフセット
  AO = (P_上半分 − P_下半分) / (P_上半分 + P_下半分)（燃料の高さの中央で分ける）
- 全長で積分したピン出力: KENO-VI の 1/8 の表（compare_keno.py と同じ読み方、対称な 8 本の平均）と、OpenMC の
  pins.msgpack（ピンごと）。差は燃料のピンの相対差の RMS と最大
- 層ごとのピン出力（OpenMC があるとき）: 参照の層ごとの、燃料のピンの (xlc / OpenMC − 1) の RMS と最大と、全体の最大
"""

import json
import math
import re
import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import itertools

from compare_keno import keno_octant

CAPTIONS = {
    "3A": "Figure P3-3: Problem 3A CE KENO-VI Radial Power Distribution Results",
    "3B": "Figure P3-4: Problem 3B CE KENO-VI Radial Power Distribution Results",
}
# 燃料の下端（仕様書の高さ、cm）
FUEL_BOTTOM = 11.951
TOLERANCE = 1e-6


def keno_k(spec: list[str], problem: str) -> tuple[float, float]:
    """付録 D の `<問題> <k> <σ>` の行"""
    for line in spec:
        m = re.match(rf"^{problem}\s+(\d\.\d+)\s+(\d\.\d+)\s*$", line)
        if m:
            return float(m.group(1)), float(m.group(2))
    raise SystemExit(f"KENO-VI k of {problem} not found")


def keno_axial(spec: list[str], problem: str) -> tuple[list[float], list[float]]:
    """付録 D の「Problem 3 Axial Power Results」の 49 層: 燃料の下端からの境界（下から）と出力"""
    start = next(
        (i for i, line in enumerate(spec) if "Problem 3 Axial Power Results" in line),
        None,
    )
    if start is None:
        raise SystemExit("the axial power table is not found")
    column = {"3A": 3, "3B": 5}[problem]
    rows = []
    for line in spec[start : start + 60]:
        m = re.match(
            r"\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+) \+/- ([\d.]+)% ([\d.]+) \+/- ([\d.]+)%",
            line,
        )
        if m:
            rows.append(
                (float(m.group(2)), float(m.group(3)), float(m.group(column + 1)))
            )
    if len(rows) != 49:
        raise SystemExit(f"expected 49 axial levels, found {len(rows)}")
    rows.sort()
    bounds = [rows[0][0] - rows[0][1] - FUEL_BOTTOM]
    for top, thickness, _ in rows:
        if abs(top - thickness - FUEL_BOTTOM - bounds[-1]) > 1e-3:
            raise SystemExit(
                f"axial level ending at {top} does not follow the one below"
            )
        bounds.append(top - FUEL_BOTTOM)
    return bounds, [r[2] for r in rows]


def read_layers(path: Path) -> dict:
    data = msgpack.unpackb(path.read_bytes())
    if data.get("format") != "xlc-axial-pin-powers":
        raise SystemExit(f"{path}: not xlc-axial-pin-powers")
    return data


def coarsen(bounds: list[float], values: list, reference: list[float]) -> list:
    """`bounds` の層の値（数か、数の入れ子の表）を、`reference` の層ごとに厚さで重みを付けて平均する。境界は燃料の
    下端から。`reference` の境界は、どれも `bounds` の境界で、両端が `bounds` の両端と同じこと（燃料の高さが同じ。
    層を捨てると、燃料の高さの平均で 1 にした規格化がそろわない）"""
    origin = bounds[0]
    z = [b - origin for b in bounds]
    index = []
    for r in reference:
        k = min(range(len(z)), key=lambda i: abs(z[i] - r))
        if abs(z[k] - r) > TOLERANCE * max(1.0, z[-1]):
            raise SystemExit(
                f"the reference bound {r} is not a bound of the xlc layers"
            )
        index.append(k)
    if index[0] != 0 or index[-1] != len(z) - 1:
        raise SystemExit(
            f"the reference layers span {reference[0]}-{reference[-1]} cm, the xlc layers {z[0]}-{z[-1]} cm"
        )

    def weighted(a: int, b: int):
        total = z[b] - z[a]

        def add(items):
            first = items[0][1]
            if isinstance(first, list):
                return [add([(w, v[i]) for w, v in items]) for i in range(len(first))]
            return sum(w * v for w, v in items) / total

        return add([(z[k + 1] - z[k], values[k]) for k in range(a, b)])

    return [weighted(a, b) for a, b in itertools.pairwise(index)]


def offset(bounds: list[float], axial: list[float]) -> float:
    middle = (bounds[0] + bounds[-1]) / 2
    lower = upper = 0.0
    for lo, hi, p in zip(bounds, bounds[1:], axial):
        lower += p * max(0.0, min(hi, middle) - lo)
        upper += p * max(0.0, hi - max(lo, middle))
    return (upper - lower) / (upper + lower)


def axial_report(
    name: str, bounds: list[float], ours: list[float], ref: list[float]
) -> None:
    if not (len(ours) == len(ref) == len(bounds) - 1):
        raise SystemExit(
            f"{len(ours)} xlc and {len(ref)} {name} levels for {len(bounds) - 1} layers"
        )
    diffs = [a - b for a, b in zip(ours, ref)]
    rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    # 層ごとの相対差 (xlc / 参照 − 1) の RMS
    relative = math.sqrt(sum((a / b - 1) ** 2 for a, b in zip(ours, ref)) / len(ref))
    worst = max(range(len(diffs)), key=lambda i: abs(diffs[i]))
    print(
        f"axial vs {name}: rms {rms:.4f} (relative {100 * relative:.2f}%), max {diffs[worst]:+.4f} at level {worst + 1} "
        f"({bounds[worst]:.3f}-{bounds[worst + 1]:.3f} cm above the fuel bottom), "
        f"AO xlc {100 * offset(bounds, ours):+.3f}% {name} {100 * offset(bounds, ref):+.3f}%"
    )
    for k, (a, b) in enumerate(zip(ours, ref)):
        print(
            f"  level {k + 1:2} {bounds[k]:8.3f}-{bounds[k + 1]:8.3f}: xlc {a:.5f} {name} {b:.5f} diff {a - b:+.5f}"
        )


def pin_stats(
    ours: list[list[float]], ref: list[list[float]]
) -> tuple[float, float, tuple]:
    """燃料のピン（どちらも正）の (xlc / 参照 − 1) の RMS、最大、その位置（行, 列）。2つの表は同じ形であること"""
    if len(ours) != len(ref) or any(len(a) != len(b) for a, b in zip(ours, ref)):
        raise SystemExit("the two pin maps differ in shape")
    diffs = [
        (a / b - 1, (i, j))
        for i, (ra, rb) in enumerate(zip(ours, ref))
        for j, (a, b) in enumerate(zip(ra, rb))
        if a > 0 and b > 0
    ]
    if any((a > 0) != (b > 0) for ra, rb in zip(ours, ref) for a, b in zip(ra, rb)):
        raise SystemExit("the two maps have fission in different cells")
    if not diffs:
        raise SystemExit("the two maps have no cell with fission")
    rms = math.sqrt(sum(d * d for d, _ in diffs) / len(diffs))
    worst = max(diffs, key=lambda d: abs(d[0]))
    return rms, worst[0], worst[1]


def keno_pins(spec: list[str], problem: str, rows: list[list[float]]) -> None:
    """compare_keno.py と同じ比べ方（燃料のピンの平均で 1、対称な 8 本の平均）"""
    if len(rows) != 17 or any(len(r) != 17 for r in rows):
        raise SystemExit("expected a 17 x 17 map")
    fuel = [x for r in rows for x in r if x > 0]
    mean = sum(fuel) / len(fuel)
    diffs = []
    for (i, j), keno in keno_octant(spec, CAPTIONS[problem]).items():
        images = {
            (8 + a, 8 + b)
            for a, b in [
                (i, j),
                (j, i),
                (-i, j),
                (j, -i),
                (i, -j),
                (-j, i),
                (-i, -j),
                (-j, -i),
            ]
        }
        ours = sum(rows[r][c] for r, c in images) / len(images) / mean
        diffs.append((ours / keno - 1, (i, j)))
    rms = math.sqrt(sum(d * d for d, _ in diffs) / len(diffs))
    worst, at = max(diffs, key=lambda d: abs(d[0]))
    print(
        f"pins vs KENO-VI: {len(diffs)} positions, rms {100 * rms:.3f}%, max {100 * worst:+.3f}% at {at}"
    )


def main() -> None:
    if len(sys.argv) not in (6, 7):
        raise SystemExit(__doc__)
    problem = sys.argv[1]
    if problem not in CAPTIONS:
        raise SystemExit(f"unknown problem {problem}")
    spec = Path(sys.argv[2]).read_text().split("\n")
    lines = [
        l for l in Path(sys.argv[3]).read_text().splitlines() if l.startswith("k-eff ")
    ]
    if len(lines) != 1:
        raise SystemExit(f"{sys.argv[3]}: expected one k-eff line")
    k = float(lines[0].split()[1])
    pins = msgpack.unpackb(Path(sys.argv[4]).read_bytes())["rows"]
    layers = read_layers(Path(sys.argv[5]))
    openmc = Path(sys.argv[6]) if len(sys.argv) == 7 else None

    k_keno, s_keno = keno_k(spec, problem)
    print(
        f"{problem} k xlc {k:.6f} KENO-VI {k_keno:.6f} +- {s_keno:.6f} diff {(k - k_keno) * 1e5:+.0f} pcm"
    )
    if openmc:
        ref = json.loads((openmc / "k.json").read_text())
        print(
            f"{problem} k xlc {k:.6f} OpenMC {ref['k']:.6f} +- {ref['k-std']:.6f} "
            f"diff {(k - ref['k']) * 1e5:+.0f} +- {ref['k-std'] * 1e5:.0f} pcm"
        )

    keno_bounds, keno_power = keno_axial(spec, problem)
    ours_axial = coarsen(layers["bounds"], layers["axial"], keno_bounds)
    axial_report("KENO-VI", keno_bounds, ours_axial, keno_power)
    keno_pins(spec, problem, pins)

    if openmc:
        reference = read_layers(openmc / "layers.msgpack")
        ref_bounds = [b - reference["bounds"][0] for b in reference["bounds"]]
        axial_report(
            "OpenMC",
            ref_bounds,
            coarsen(layers["bounds"], layers["axial"], ref_bounds),
            reference["axial"],
        )
        ref_pins = msgpack.unpackb((openmc / "pins.msgpack").read_bytes())["rows"]
        rms, worst, at = pin_stats(pins, ref_pins)
        print(f"pins vs OpenMC: rms {100 * rms:.3f}%, max {100 * worst:+.3f}% at {at}")
        ours_layers = coarsen(layers["bounds"], layers["layers"], ref_bounds)
        if len(reference["layers"]) != len(ours_layers):
            raise SystemExit(
                f"{len(reference['layers'])} OpenMC layer maps for {len(ours_layers)} layers"
            )
        overall = None
        for k_level, (a, b) in enumerate(zip(ours_layers, reference["layers"])):
            rms, worst, at = pin_stats(a, b)
            print(
                f"  level {k_level + 1:2}: pins rms {100 * rms:.3f}%, max {100 * worst:+.3f}% at {at}"
            )
            if overall is None or abs(worst) > abs(overall[0]):
                overall = (worst, k_level + 1, at)
        if overall is None:
            raise SystemExit("the OpenMC results have no layers")
        print(
            f"layer pins vs OpenMC: max {100 * overall[0]:+.3f}% at level {overall[1]} {overall[2]}"
        )


if __name__ == "__main__":
    main()
