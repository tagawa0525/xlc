"""VERA Problem 5-2D の xlc の炉心の解を、参照（KENO-VI、OpenMC、最小の入力の MOC）と比べ、炉心の解の対称性を調べる。

使い方:
  python3 compare_p5.py keno <問題（5A、5B、5C）> <spec.txt> <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py openmc <core の標準出力> <core の pins.msgpack> <OpenMC の結果のディレクトリ>
  python3 compare_p5.py moc <lattice の標準出力> <lattice の pins.msgpack> <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py symmetry <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py k <lattice の標準出力> <core の標準出力>

- k: core の標準出力の `k-eff`。KENO-VI は付録 J の表（差は (k_xlc / k_ref − 1) × 1e5 pcm）、OpenMC は結果の k.json、
  MOC は lattice の標準出力の `k`（外周が真空の格子では実効増倍率）
- 集合体の出力: core の標準出力の集合体ごとの核分裂率（集合体の平均で 1）。KENO-VI は付録 J の 1/8 の表（5 桁）を、H 列と
  8 行で鏡映して全炉心に広げる。OpenMC は assemblies.msgpack（15 × 15）。MOC はピンの核分裂率の 17 × 17 の和を集合体の
  平均で割る。差は (xlc − 参照)（出力の単位のまま）と相対差の RMS と最大
- ピン出力: どちらも燃料のピンの平均で 1 に規格化し直し、燃料のない行と列（集合体の間の余白、反射体）を除いて重ねる。
  相対差 (xlc / 参照 − 1) の RMS と最大、炉心の外周の集合体（反射体に面か角の点で接する集合体）のピンと、出力 0.8 未満の
  ピンに分けた RMS と最大。KENO-VI はピン出力を公開していないので、最大のピン出力（表 M3-6）だけを並べる
- 対称性: 集合体の出力とピン出力を、全炉心の中心のまわりの 90° の回転と、x と y の鏡映で移したものとの差の最大
"""

import json
import re
import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_k

# 付録 J の見出し（Problem 5-2D ENDF/B-VII.0 Assembly Power Results）
ASSEMBLY_TABLE = "{problem}-2D Assembly Powers"
# 表 M3-6 の最大のピン出力
KENO_MAX_PIN = {"5A": 1.44402, "5B": 1.60240, "5C": 1.62997}
QUARTER = 8


def core_output(path: str) -> tuple[float, np.ndarray]:
    """core の標準出力の k と、集合体の位置ごとの核分裂率（行は上から。反射体と `-` の位置は NaN）"""
    k = None
    rows = []
    reading = False
    for line in Path(path).read_text().splitlines():
        if line.startswith("k-eff"):
            if k is not None:
                raise SystemExit(f"{path}: more than one k-eff")
            k = float(line.split()[1])
        elif line.startswith("(`.` is a reflector site"):
            reading = True
        elif reading and line.strip():
            try:
                rows.append(
                    [float(x) if x not in (".", "-") else np.nan for x in line.split()]
                )
            except ValueError:
                raise SystemExit(
                    f"{path}: not a row of assembly fission rates: {line}"
                ) from None
    if k is None or not rows:
        raise SystemExit(f"{path}: no k-eff or assembly fission rates")
    if len({len(r) for r in rows}) != 1:
        raise SystemExit(f"{path}: assembly rows differ in length")
    return k, np.array(rows)


def fuel_assemblies(rates: np.ndarray) -> np.ndarray:
    """反射体の輪を除いた集合体の位置（燃料のある行と列）"""
    keep_r = ~np.all(np.isnan(rates), axis=1)
    keep_c = ~np.all(np.isnan(rates), axis=0)
    return rates[keep_r][:, keep_c]


def read_pins(path: str) -> np.ndarray:
    rows = np.array(
        msgpack.unpackb(Path(path).read_bytes(), raw=False)["rows"], dtype=float
    )
    rows = rows[rows.any(axis=1)][:, rows.any(axis=0)]
    return rows / rows[rows > 0].mean()


def keno_k(spec: list[str], problem: str) -> tuple[float, float]:
    for line in spec:
        m = re.match(rf"^{problem}-2D\s+(\d\.\d+)\s+(\d\.\d+)\s*$", line)
        if m:
            return float(m.group(1)), float(m.group(2))
    raise SystemExit(f"KENO-VI k of {problem}-2D not found")


def keno_assemblies(spec: list[str], problem: str) -> np.ndarray:
    """付録 J の 1/8 の表（行 8〜15、H 列から）を全炉心の 15 × 15 に広げる"""
    title = ASSEMBLY_TABLE.format(problem=problem)
    start = next(
        (i for i, line in enumerate(spec) if line.strip().startswith(title)), None
    )
    if start is None:
        raise SystemExit(f"{title} not found")
    octant = []
    for line in spec[start + 2 : start + 2 + QUARTER]:
        values = [float(x) for x in re.findall(r"(?<![\d.])(\d\.\d{5})(?![\d%])", line)]
        octant.append(values)
    quarter = np.full((QUARTER, QUARTER), np.nan)
    for r, values in enumerate(octant):
        for c, v in enumerate(values):
            quarter[r, c] = v
    # 1/8 の表の欠けた位置（下三角の外）を、対角の鏡映で埋める
    for r in range(QUARTER):
        for c in range(QUARTER):
            if np.isnan(quarter[r, c]) and not np.isnan(quarter[c, r]):
                quarter[r, c] = quarter[c, r]
    if np.count_nonzero(~np.isnan(quarter)) != 56:
        raise SystemExit(
            f"{title}: expected 56 assemblies in the quarter, read {np.count_nonzero(~np.isnan(quarter))}"
        )
    return unfold(quarter)


def unfold(quarter: np.ndarray) -> np.ndarray:
    """右下の 1/4（先頭の行と列が対称面の上）を鏡面反射で全体に広げる（openmc_core_2d.py と同じ）"""
    rows = np.concatenate([quarter[:0:-1], quarter], axis=0)
    return np.concatenate([rows[:, :0:-1], rows], axis=1)


def summary(name: str, d: np.ndarray, where: np.ndarray | None = None) -> None:
    if d.size == 0:
        print(f"{name}: none")
        return
    i = np.abs(d).argmax()
    at = f" at {tuple(int(x) for x in where[i])}" if where is not None else ""
    print(
        f"{name}: {d.size}, rms {100 * np.sqrt((d * d).mean()):.3f}%, max {100 * d[i]:+.3f}%{at}"
    )


def compare_assemblies(
    ours: np.ndarray, reference: np.ndarray, std: np.ndarray | None = None
) -> None:
    ours = fuel_assemblies(ours)
    if ours.shape != reference.shape:
        raise SystemExit(
            f"assembly maps differ in shape: {ours.shape} and {reference.shape}"
        )
    fuel = ~np.isnan(reference) & (reference > 0)
    if not np.array_equal(fuel, ~np.isnan(ours)):
        raise SystemExit("assembly maps differ in fuel positions")
    ours = ours / ours[fuel].mean()
    reference = reference / reference[fuel].mean()
    rel = ours[fuel] / reference[fuel] - 1
    absolute = ours[fuel] - reference[fuel]
    where = np.argwhere(fuel)
    summary("assemblies (relative)", rel, where)
    i = np.abs(absolute).argmax()
    print(
        f"assemblies (absolute): rms {np.sqrt((absolute**2).mean()):.4f}, max {absolute[i]:+.4f} "
        f"at {tuple(int(x) for x in where[i])} (reference {reference[fuel][i]:.4f})"
    )
    edge = edge_assemblies(fuel)
    summary(
        "peripheral assemblies (relative)", rel[edge[fuel]], np.argwhere(fuel & edge)
    )
    if std is not None:
        s = std[fuel] / reference[fuel]
        print(
            f"reference relative std: median {100 * np.median(s):.3f}%, max {100 * s.max():.3f}%"
        )
    print("quarter (rows 8-15 from the top, columns H-A), xlc / reference:")
    n = fuel.shape[0] // 2
    for r in range(n, fuel.shape[0]):
        print(
            " ".join(
                f"{ours[r, c]:.4f}/{reference[r, c]:.4f}" if fuel[r, c] else " " * 13
                for c in range(n, fuel.shape[1])
            )
        )


def edge_assemblies(fuel: np.ndarray) -> np.ndarray:
    """燃料でない位置か図の外に、面か角の点で接する集合体"""
    padded = np.pad(fuel, 1, constant_values=False)
    out = np.zeros_like(fuel)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            out |= ~padded[
                1 + di : 1 + di + fuel.shape[0], 1 + dj : 1 + dj + fuel.shape[1]
            ]
    return out & fuel


def compare_pins(
    ours: np.ndarray, reference: np.ndarray, std: np.ndarray | None = None
) -> None:
    if ours.shape != reference.shape:
        raise SystemExit(
            f"pin maps differ in shape: {ours.shape} and {reference.shape}"
        )
    fuel = reference > 0
    if not np.array_equal(fuel, ours > 0):
        raise SystemExit("the two pin maps have different fuel cells")
    d = ours[fuel] / reference[fuel] - 1
    where = np.argwhere(fuel)
    summary("pins", d, where)
    n = 17
    blocks = fuel.reshape(fuel.shape[0] // n, n, fuel.shape[1] // n, n).any(axis=(1, 3))
    edge = np.kron(edge_assemblies(blocks), np.ones((n, n), dtype=bool)) & fuel
    summary(
        "pins of peripheral assemblies",
        ours[edge] / reference[edge] - 1,
        np.argwhere(edge),
    )
    low = fuel & (reference < 0.8)
    summary("pins with power < 0.8", ours[low] / reference[low] - 1, np.argwhere(low))
    print(f"max pin power: xlc {ours.max():.5f}, reference {reference.max():.5f}")
    if std is not None:
        s = std[fuel] / reference[fuel]
        print(
            f"reference relative std: median {100 * np.median(s):.3f}%, max {100 * s.max():.3f}%"
        )


def keno(problem: str, spec_path: str, core_out: str, core_pins: str) -> None:
    spec = Path(spec_path).read_text().splitlines()
    k, rates = core_output(core_out)
    k_ref, sigma = keno_k(spec, problem)
    print(
        f"k: xlc {k:.6f}, KENO-VI {k_ref:.6f} +- {sigma:.6f}, {1e5 * (k / k_ref - 1):+.1f} pcm"
    )
    compare_assemblies(rates, keno_assemblies(spec, problem))
    pins = read_pins(core_pins)
    print(f"max pin power: xlc {pins.max():.5f}, KENO-VI {KENO_MAX_PIN[problem]:.5f}")


def openmc(core_out: str, core_pins: str, results: str) -> None:
    k, rates = core_output(core_out)
    ref = json.loads((Path(results) / "k.json").read_text())
    print(
        f"k: xlc {k:.6f}, OpenMC {ref['k']:.6f} +- {ref['k-std']:.6f}, {1e5 * (k / ref['k'] - 1):+.1f} pcm"
    )
    assemblies = np.array(
        msgpack.unpackb((Path(results) / "assemblies.msgpack").read_bytes())["rows"]
    )
    std = np.array(
        msgpack.unpackb((Path(results) / "assemblies-std.msgpack").read_bytes())["rows"]
    )
    compare_assemblies(rates, np.where(assemblies > 0, assemblies, np.nan), std)
    pins = np.array(
        msgpack.unpackb((Path(results) / "pins.msgpack").read_bytes())["rows"]
    )
    pins_std = np.array(
        msgpack.unpackb((Path(results) / "pins-std.msgpack").read_bytes())["rows"]
    )
    mean = pins[pins > 0].mean()
    compare_pins(read_pins(core_pins), pins / mean, pins_std / mean)


def lattice_k(path: str) -> float:
    return compare_k.lattice_k(Path(path))


def moc(lattice_out: str, lattice_pins: str, core_out: str, core_pins: str) -> None:
    k_ref = lattice_k(lattice_out)
    k, rates = core_output(core_out)
    print(f"k: core {k:.6f}, MOC {k_ref:.6f}, {1e5 * (k / k_ref - 1):+.1f} pcm")
    reference = read_pins(lattice_pins)
    n = 17
    blocks = reference.reshape(
        reference.shape[0] // n, n, reference.shape[1] // n, n
    ).sum(axis=(1, 3))
    compare_assemblies(rates, blocks / blocks.mean())
    compare_pins(read_pins(core_pins), reference)


def symmetry(core_out: str, core_pins: str) -> None:
    k, rates = core_output(core_out)
    print(f"k: {k:.6f}")
    a = np.nan_to_num(fuel_assemblies(rates))
    p = read_pins(core_pins)
    for name, m in (("assemblies", a), ("pins", p)):
        scale = np.abs(m).max()
        for op, image in (
            ("rotation by 90", np.rot90(m)),
            ("rotation by 180", np.rot90(m, 2)),
            ("mirror in x", m[:, ::-1]),
            ("mirror in y", m[::-1, :]),
            ("mirror in the diagonal", m.T),
        ):
            print(
                f"{name}, {op}: max |difference| / max {np.abs(image - m).max() / scale:.2e}"
            )


def main() -> None:
    match sys.argv[1:]:
        case ["keno", problem, spec, out, pins]:
            keno(problem, spec, out, pins)
        case ["openmc", out, pins, results]:
            openmc(out, pins, results)
        case ["moc", lattice_out, lattice_pins, out, pins]:
            moc(lattice_out, lattice_pins, out, pins)
        case ["symmetry", out, pins]:
            symmetry(out, pins)
        case ["k", lattice_out, out]:
            k_ref = lattice_k(lattice_out)
            k, _ = core_output(out)
            print(f"k: core {k:.6f}, MOC {k_ref:.6f}, {1e5 * (k / k_ref - 1):+.1f} pcm")
        case _:
            raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
