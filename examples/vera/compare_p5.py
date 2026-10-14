"""VERA Problem 5-2D の xlc の炉心の解を、参照（KENO-VI、OpenMC、最小の入力の MOC）と比べ、炉心の解の対称性を調べる。

使い方:
  python3 compare_p5.py --json <要約の出力先> keno <問題（5A、5B、5C）> <spec.txt> <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py --json <要約の出力先> openmc <問題> <core の標準出力> <core の pins.msgpack> <OpenMC の結果のディレクトリ>
  python3 compare_p5.py --report <要約の JSON>
  python3 compare_p5.py moc <lattice の標準出力> <lattice の pins.msgpack> <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py symmetry <core の標準出力> <core の pins.msgpack>
  python3 compare_p5.py k <lattice の標準出力> <core の標準出力>

- `keno` と `openmc` は比べた結果を JSON だけに書く（標準出力には何も出さない）。JSON は検証の記録の表を描くツール
  （examples/verification_tables.py）も読む（ADR-0176）。割合は小数（0.01 = 1%）、k の差は pcm、位置は [行, 列]（集合体は 15 × 15 の
  地図の添字で、名前は列を西から R P N M L K J H G F E D C B A、行を北から 01 とする）。`--report` は JSON から文章の報告を標準出力に描く。
  1 行目は読んだファイルのパスと md5 の先頭 8 桁（`inputs`）で、JSON が持つ数だけを使う。OpenMC との比較では、差が OpenMC の統計誤差
  （`assemblies-std.msgpack`、`pins-std.msgpack`、`k.json` の `k-std`）の何倍かも書く。集合体とピンは、差の最大の倍率とその位置、2σ を超える数
- k: core の標準出力の `k-eff`。KENO-VI は付録 J の表（差は (k_xlc / k_ref − 1) × 1e5 pcm）、OpenMC は結果の k.json、
  MOC は lattice の標準出力の `k`（外周が真空の格子では実効増倍率）
- 集合体の出力: core の標準出力の集合体ごとの核分裂率（集合体の平均で 1）。KENO-VI は付録 J の 1/8 の表（5 桁）を、H 列と
  8 行で鏡映して全炉心に広げる。OpenMC は assemblies.msgpack（15 × 15）。MOC はピンの核分裂率の 17 × 17 の和を集合体の
  平均で割る。差は (xlc − 参照)（出力の単位のまま）と相対差の RMS と最大
- ピン出力: どちらも燃料のピンの平均で 1 に規格化し直し、燃料のない行と列（集合体の間の余白、反射体）を除いて重ねる。
  相対差 (xlc / 参照 − 1) の RMS と最大、炉心の外周の集合体（反射体に面か角の点で接する集合体）のピンと、出力 0.8 未満の
  ピンに分けた RMS と最大。KENO-VI はピン出力を公開していないので、最大のピン出力（表 M3-6）だけを並べる
- 対称性: 集合体の出力とピン出力を、全炉心の中心のまわりの 90° と 180° の回転と、x、y、対角の鏡映で移したものとの差の最大
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_k

FORMAT = "xlc-compare-p5"
# 全炉心の列の名前（西から。I、O、Q は使わない。p5_inputs.py と同じ）
COLUMNS = "RPNMLKJHGFEDCBA"
# 付録 J の見出し（Problem 5-2D ENDF/B-VII.0 Assembly Power Results）
ASSEMBLY_TABLE = "{problem}-2D Assembly Powers"
# 表 M3-6 の見出しと、問題の列の順
POWER_TABLE = "Table M3-6: Problem 5-2D Reference Solution Power Statistics"
PROBLEMS = ["5A", "5B", "5C"]
QUARTER = 8
# 付録 J の 1/8 の表の行（8〜15 行）ごとの集合体の数（H 列から。図 M3-1 の 1/4 炉心の下三角）
OCTANT_ROWS = [1, 2, 3, 4, 5, 6, 6, 4]


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


def keno_max_pin(spec: list[str], problem: str) -> float:
    """表 M3-6 のピン出力の最大（集合体の最大の行の次の `Maximum` の行）"""
    # 目次にも同じ見出しがあるので、`Maximum` の行が続くものを取る
    start = next(
        (
            i
            for i, line in enumerate(spec)
            if POWER_TABLE in line
            and any(l.split()[:1] == ["Maximum"] for l in spec[i : i + 4])
        ),
        None,
    )
    if start is None:
        raise SystemExit(f"{POWER_TABLE} not found")
    rows = [
        line.split()[1:]
        for line in spec[start : start + 12]
        if line.split()[:1] == ["Maximum"] and line.split()[1] != "Uncertainty"
    ]
    if len(rows) != 2 or len(rows[1]) != len(PROBLEMS):
        raise SystemExit(
            f"{POWER_TABLE}: expected two Maximum rows of {len(PROBLEMS)} values"
        )
    return float(rows[1][PROBLEMS.index(problem)])


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
        if len(values) != OCTANT_ROWS[len(octant)]:
            raise SystemExit(
                f"{title}: row {8 + len(octant)} has {len(values)} values, expected {OCTANT_ROWS[len(octant)]}"
            )
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


def stat(d: np.ndarray, where: np.ndarray | None = None) -> dict:
    """相対差の数、RMS、絶対値が最大のもの（符号つき）とその位置 [行, 列]"""
    if d.size == 0:
        return {"n": 0, "rms": None, "max": None, "at": None}
    i = int(np.abs(d).argmax())
    return {
        "n": int(d.size),
        "rms": float(np.sqrt((d * d).mean())),
        "max": float(d[i]),
        "at": [int(x) for x in where[i]] if where is not None else None,
    }


def stat_line(name: str, d: dict) -> str:
    if d["n"] == 0:
        return f"{name}: none"
    at = f" at {tuple(d['at'])}" if d["at"] is not None else ""
    return (
        f"{name}: {d['n']}, rms {100 * d['rms']:.3f}%, max {100 * d['max']:+.3f}%{at}"
    )


def std_summary(std: np.ndarray, reference: np.ndarray) -> dict:
    """参照の相対標準偏差（標準偏差 / 参照）の中央値と最大"""
    s = std / reference
    return {"median": float(np.median(s)), "max": float(s.max())}


def sigma_summary(
    ours: np.ndarray, reference: np.ndarray, std: np.ndarray, where: np.ndarray
) -> dict:
    """|xlc − 参照| / 参照の標準偏差 の最大とその位置、2σ を超える数（等しいときは超えない）"""
    ratio = np.abs(ours - reference) / std
    i = int(ratio.argmax())
    return {
        "max": float(ratio[i]),
        "at": [int(x) for x in where[i]],
        "beyond_2sigma": int((ratio > 2).sum()),
        "n": int(ratio.size),
    }


def assembly_summary(
    ours: np.ndarray, reference: np.ndarray, std: np.ndarray | None = None
) -> dict:
    """集合体の出力の差。どちらも燃料の集合体の平均で 1 に規格化する（標準偏差も参照の平均で割る）"""
    ours = fuel_assemblies(ours)
    if ours.shape != reference.shape:
        raise SystemExit(
            f"assembly maps differ in shape: {ours.shape} and {reference.shape}"
        )
    fuel = ~np.isnan(reference) & (reference > 0)
    if not np.array_equal(fuel, ~np.isnan(ours)):
        raise SystemExit("assembly maps differ in fuel positions")
    ours = ours / ours[fuel].mean()
    mean = reference[fuel].mean()
    reference = reference / mean
    if std is not None:
        std = std / mean
    rel = ours[fuel] / reference[fuel] - 1
    absolute = ours[fuel] - reference[fuel]
    where = np.argwhere(fuel)
    i = int(np.abs(absolute).argmax())
    edge = edge_assemblies(fuel)
    n = fuel.shape[0] // 2
    result = {
        "relative": stat(rel, where),
        "absolute": {
            "rms": float(np.sqrt((absolute**2).mean())),
            "max": float(absolute[i]),
            "at": [int(x) for x in where[i]],
            "reference": float(reference[fuel][i]),
        },
        "peripheral": stat(rel[edge[fuel]], np.argwhere(fuel & edge)),
        # 中心の行と列から右下の 1/4 の [xlc, 参照]（燃料でない位置は null）
        "quarter": [
            [
                [float(ours[r, c]), float(reference[r, c])] if fuel[r, c] else None
                for c in range(n, fuel.shape[1])
            ]
            for r in range(n, fuel.shape[0])
        ],
    }
    if std is not None:
        result["std"] = std_summary(std[fuel], reference[fuel])
        result["sigma"] = sigma_summary(ours[fuel], reference[fuel], std[fuel], where)
    return result


def assembly_lines(d: dict) -> list[str]:
    a = d["absolute"]
    lines = [
        stat_line("assemblies (relative)", d["relative"]),
        (
            f"assemblies (absolute): rms {a['rms']:.4f}, max {a['max']:+.4f} "
            f"at {tuple(a['at'])} (reference {a['reference']:.4f})"
        ),
        stat_line("peripheral assemblies (relative)", d["peripheral"]),
    ]
    if "std" in d:
        lines.append(
            f"reference relative std: median {100 * d['std']['median']:.3f}%, max {100 * d['std']['max']:.3f}%"
        )
    lines.append(
        "lower-right quarter from the centre row and column (rows 8-15, columns H-A in 5-2D), xlc / reference:"
    )
    for row in d["quarter"]:
        lines.append(" ".join(f"{c[0]:.4f}/{c[1]:.4f}" if c else " " * 13 for c in row))
    if "sigma" in d:
        g = d["sigma"]
        lines.append(
            f"assemblies xlc - reference in sigma: max {g['max']:.1f} at {tuple(g['at'])}, "
            f"{g['beyond_2sigma']} of {g['n']} beyond 2 sigma"
        )
    return lines


def assembly_name(at: list[int]) -> str:
    """15 × 15 の集合体の位置 [行, 列] の名前（列は西から R P N M L K J H G F E D C B A、行は北から 01）"""
    return f"{COLUMNS[at[1]]}{at[0] + 1:02d}"


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


def pin_summary(
    ours: np.ndarray, reference: np.ndarray, std: np.ndarray | None = None
) -> dict:
    if ours.shape != reference.shape:
        raise SystemExit(
            f"pin maps differ in shape: {ours.shape} and {reference.shape}"
        )
    fuel = reference > 0
    if not np.array_equal(fuel, ours > 0):
        raise SystemExit("the two pin maps have different fuel cells")
    n = 17
    blocks = fuel.reshape(fuel.shape[0] // n, n, fuel.shape[1] // n, n).any(axis=(1, 3))
    edge = np.kron(edge_assemblies(blocks), np.ones((n, n), dtype=bool)) & fuel
    low = fuel & (reference < 0.8)

    def relative(mask: np.ndarray) -> dict:
        return stat(ours[mask] / reference[mask] - 1, np.argwhere(mask))

    result = {
        "all": relative(fuel),
        "peripheral": relative(edge),
        "low": relative(low),
        "max_pin": {"xlc": float(ours.max()), "reference": float(reference.max())},
    }
    if std is not None:
        result["std"] = std_summary(std[fuel], reference[fuel])
        result["sigma"] = sigma_summary(
            ours[fuel], reference[fuel], std[fuel], np.argwhere(fuel)
        )
    return result


def pin_lines(d: dict) -> list[str]:
    lines = [
        stat_line("pins", d["all"]),
        stat_line("pins of peripheral assemblies", d["peripheral"]),
        stat_line("pins with power < 0.8", d["low"]),
        f"max pin power: xlc {d['max_pin']['xlc']:.5f}, reference {d['max_pin']['reference']:.5f}",
    ]
    if "std" in d:
        lines.append(
            f"reference relative std: median {100 * d['std']['median']:.3f}%, max {100 * d['std']['max']:.3f}%"
        )
        g = d["sigma"]
        lines.append(
            f"pins xlc - reference in sigma: max {g['max']:.1f} at {tuple(g['at'])}, "
            f"{g['beyond_2sigma']} of {g['n']} beyond 2 sigma"
        )
    return lines


def md5_inputs(paths: list[Path]) -> list[dict]:
    """読んだファイルのパス（段階の cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）"""
    return [
        {"path": str(p), "md5": hashlib.md5(p.read_bytes()).hexdigest()} for p in paths
    ]


def keno(problem: str, spec_path: str, core_out: str, core_pins: str) -> dict:
    spec = Path(spec_path).read_text().splitlines()
    k, rates = core_output(core_out)
    k_ref, sigma = keno_k(spec, problem)
    pins = read_pins(core_pins)
    return {
        "format": FORMAT,
        "kind": "keno",
        "problem": problem,
        "inputs": md5_inputs([Path(spec_path), Path(core_out), Path(core_pins)]),
        "k": {
            "xlc": k,
            "reference": k_ref,
            "reference_std": sigma,
            "diff_pcm": 1e5 * (k / k_ref - 1),
        },
        "assemblies": assembly_summary(rates, keno_assemblies(spec, problem)),
        "max_pin": {"xlc": float(pins.max()), "reference": keno_max_pin(spec, problem)},
    }


def openmc(problem: str, core_out: str, core_pins: str, results: str) -> dict:
    k, rates = core_output(core_out)
    directory = Path(results)
    files = [
        directory / name
        for name in [
            "k.json",
            "assemblies.msgpack",
            "assemblies-std.msgpack",
            "pins.msgpack",
            "pins-std.msgpack",
        ]
    ]
    ref = json.loads(files[0].read_text())
    assemblies, std, pins, pins_std = (
        np.array(msgpack.unpackb(f.read_bytes())["rows"]) for f in files[1:]
    )
    # xlc の表と同じく、燃料のない行と列を除いてから規格化する（統計誤差も同じ行と列）
    keep_r, keep_c = pins.any(axis=1), pins.any(axis=0)
    pins, pins_std = (m[keep_r][:, keep_c] for m in (pins, pins_std))
    mean = pins[pins > 0].mean()
    diff_pcm = 1e5 * (k / ref["k"] - 1)
    return {
        "format": FORMAT,
        "kind": "openmc",
        "problem": problem,
        "inputs": md5_inputs([Path(core_out), Path(core_pins), *files]),
        "k": {
            "xlc": k,
            "reference": ref["k"],
            "reference_std": ref["k-std"],
            "diff_pcm": diff_pcm,
            "diff_std_pcm": 1e5 * ref["k-std"],
            "sigma": abs(diff_pcm) / (1e5 * ref["k-std"]),
        },
        "assemblies": assembly_summary(
            rates, np.where(assemblies > 0, assemblies, np.nan), std
        ),
        "pins": pin_summary(read_pins(core_pins), pins / mean, pins_std / mean),
    }


def report(summary: dict) -> list[str]:
    """要約の JSON から文章の報告を描く"""
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs (an older run); rerun the stage")
    k = summary["k"]
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    if summary["kind"] == "keno":
        lines.append(
            f"k: xlc {k['xlc']:.6f}, KENO-VI {k['reference']:.6f} +- {k['reference_std']:.6f}, {k['diff_pcm']:+.1f} pcm"
        )
        lines += assembly_lines(summary["assemblies"])
        m = summary["max_pin"]
        lines.append(f"max pin power: xlc {m['xlc']:.5f}, KENO-VI {m['reference']:.5f}")
        return lines
    lines.append(
        f"k: xlc {k['xlc']:.6f}, OpenMC {k['reference']:.6f} +- {k['reference_std']:.6f}, {k['diff_pcm']:+.1f} pcm"
    )
    lines.append(f"k xlc - OpenMC in sigma: {k['sigma']:.1f}")
    lines += assembly_lines(summary["assemblies"])
    lines += pin_lines(summary["pins"])
    return lines


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
    print("\n".join(assembly_lines(assembly_summary(rates, blocks / blocks.mean()))))
    print("\n".join(pin_lines(pin_summary(read_pins(core_pins), reference))))


def symmetry(core_out: str, core_pins: str) -> None:
    k, rates = core_output(core_out)
    print(f"k: {k:.6f}")
    a = np.nan_to_num(fuel_assemblies(rates))
    p = read_pins(core_pins)
    for name, m in (("assemblies", a), ("pins", p)):
        if m.shape[0] != m.shape[1]:
            raise SystemExit(
                f"{name}: the map is not square ({m.shape}), so it has no 90-degree symmetry to check"
            )
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
    args = sys.argv[1:]
    match args:
        case ["--report", path]:
            summary = json.loads(Path(path).read_text())
            if summary.get("format") != FORMAT:
                raise SystemExit(f"{path}: not {FORMAT}")
            print("\n".join(report(summary)))
        case ["--json", out, "keno", problem, spec, core_out, pins]:
            Path(out).write_text(
                json.dumps(keno(problem, spec, core_out, pins), indent=1) + "\n"
            )
        case ["--json", out, "openmc", problem, core_out, pins, results]:
            Path(out).write_text(
                json.dumps(openmc(problem, core_out, pins, results), indent=1) + "\n"
            )
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
