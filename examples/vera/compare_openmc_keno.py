"""ENDF/B-VII.0 の OpenMC（NJOY で処理）を、同じ ENDF/B-VII.0 のライブラリの xlc と、仕様書の KENO-VI と比べる（#0131）。

使い方: python3 compare_openmc_keno.py --json <要約の出力先> <spec.txt> <問題> <xlc の compare-<問題>.json> <OpenMC の結果のディレクトリ>
                                       [<xlc の pins.msgpack>]

xlc の k は endfb70/compare-<問題>.json（xlc-compare-keno-lattice）の `k.xlc`、OpenMC の k と標準偏差は結果のディレクトリの
k.json、KENO-VI の k は付録の行（compare_p3.py の keno_k）。k の差は (k_A − k_B) × 1e5 pcm で、誤差は OpenMC の
標準偏差（xlc は決定論なので誤差なし）と KENO-VI の σ の二乗和の平方根。集合体（2 で始まる問題）は、OpenMC のピン出力
（pins.msgpack）を KENO-VI の 1/8 の表と比べ（compare_keno.py と同じ読み方）、xlc のピン出力（引数があるとき）を OpenMC の
ピン出力と比べる（どちらも燃料ピンの平均で 1 に規格化、相対差 (A / B − 1) の二乗平均と最大）。

要約（JSON。ADR-0176）: `{"format": "xlc-compare-openmc-keno", "problem", "inputs": [{"path", "md5"}], "k": {"xlc", "openmc",
"openmc_std", "keno", "keno_std", "diff_xlc_openmc_pcm", "diff_xlc_openmc_std_pcm", "diff_openmc_keno_pcm",
"diff_openmc_keno_std_pcm", "diff_xlc_keno_pcm"}, "pins_openmc_vs_keno": {...} | null, "pins_xlc_vs_openmc": {...} | null}`。
割合は小数、位置は [行, 列]。標準出力には何も出さない。
"""

import hashlib
import json
import math
import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_keno import compare
from compare_p3 import keno_k


def normalized(path: Path) -> list[list[float]]:
    rows = msgpack.unpackb(path.read_bytes())["rows"]
    fuel = [x for row in rows for x in row if x > 0]
    mean = sum(fuel) / len(fuel)
    return [[x / mean for x in row] for row in rows]


def pin_differences(ours: Path, reference: Path) -> dict:
    """燃料ピン（どちらも核分裂率が正のセル）の相対差 (ours / reference − 1) の二乗平均と、絶対値が最大の差とその位置"""
    a, b = normalized(ours), normalized(reference)
    if [[x > 0 for x in r] for r in a] != [[x > 0 for x in r] for r in b]:
        raise SystemExit("the two maps have fission in different cells")
    diffs = [
        (x / y - 1, [i, j])
        for i, (ra, rb) in enumerate(zip(a, b))
        for j, (x, y) in enumerate(zip(ra, rb))
        if x > 0
    ]
    worst, at = max(diffs, key=lambda d: abs(d[0]))
    rms = math.sqrt(sum(d * d for d, _ in diffs) / len(diffs))
    return {"positions": len(diffs), "rms": rms, "max": worst, "at": at}


def main() -> None:
    args = sys.argv[1:]
    if len(args) not in (6, 7) or args[0] != "--json":
        raise SystemExit(__doc__)
    spec_path, problem = Path(args[2]), args[3]
    xlc_path, results = Path(args[4]), Path(args[5])
    xlc_pins = Path(args[6]) if len(args) == 7 else None
    spec = spec_path.read_text().split("\n")
    xlc = json.loads(xlc_path.read_text())
    omc = json.loads((results / "k.json").read_text())
    k_keno, s_keno = keno_k(spec, problem)
    k_xlc, k_omc, s_omc = xlc["k"]["xlc"], omc["k"], omc["k-std"]
    files = [spec_path, xlc_path, results / "k.json"]
    pins_keno = pins_xlc = None
    if problem.startswith("2"):
        omc_pins = results / "pins.msgpack"
        files.append(omc_pins)
        rows = msgpack.unpackb(omc_pins.read_bytes())["rows"]
        n, rms, worst, at = compare(spec, problem, rows)
        pins_keno = {"positions": n, "rms": rms, "max": worst, "at": list(at)}
        if xlc_pins is not None:
            files.append(xlc_pins)
            pins_xlc = pin_differences(xlc_pins, omc_pins)
    summary = {
        "format": "xlc-compare-openmc-keno",
        "problem": problem,
        # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する
        "inputs": [
            {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
            for f in files
        ],
        "k": {
            "xlc": k_xlc,
            "openmc": k_omc,
            "openmc_std": s_omc,
            "keno": k_keno,
            "keno_std": s_keno,
            "diff_xlc_openmc_pcm": (k_xlc - k_omc) * 1e5,
            "diff_xlc_openmc_std_pcm": s_omc * 1e5,
            "diff_openmc_keno_pcm": (k_omc - k_keno) * 1e5,
            "diff_openmc_keno_std_pcm": math.hypot(s_omc, s_keno) * 1e5,
            "diff_xlc_keno_pcm": (k_xlc - k_keno) * 1e5,
        },
        "pins_openmc_vs_keno": pins_keno,
        "pins_xlc_vs_openmc": pins_xlc,
    }
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
