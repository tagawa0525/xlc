"""xlc の格子計算（ピンセル 1B〜1D、集合体 2A〜2K、2Q）の k とピン出力を、仕様書の KENO-VI（ENDF/B-VII.0）の値と比べる。

使い方: python3 compare_keno_lattice.py --json <要約の出力先> <spec.txt> <問題> <lattice の出力> [<pins.msgpack>]

k は付録の `<問題> <k> <σ>` の行（compare_p3.py の keno_k と同じ）、ピン出力は compare_keno.py と同じ読み方（1/8 の表、
対称な 8 本の平均）。差は (k_xlc − k_KENO-VI) × 1e5 pcm。ピンセルはピン出力を比べない（`pins` は null）。

要約（JSON。ADR-0176）: `{"format": "xlc-compare-keno-lattice", "problem", "k": {"xlc", "keno", "keno_std", "diff_pcm"},
"pins": {"positions", "rms", "max", "at"} | null}`。割合は小数、位置は [行, 列]。標準出力には何も出さない。
"""

import json
import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_k import lattice_k
from compare_keno import compare
from compare_p3 import keno_k


def summary(
    spec: list[str], problem: str, lattice_out: Path, pins: Path | None
) -> dict:
    k = lattice_k(lattice_out)
    k_keno, sigma = keno_k(spec, problem)
    result = {
        "format": "xlc-compare-keno-lattice",
        "problem": problem,
        "k": {
            "xlc": k,
            "keno": k_keno,
            "keno_std": sigma,
            "diff_pcm": (k - k_keno) * 1e5,
        },
        "pins": None,
    }
    if pins is not None:
        rows = msgpack.unpackb(pins.read_bytes())["rows"]
        n, rms, worst, at = compare(spec, problem, rows)
        result["pins"] = {"positions": n, "rms": rms, "max": worst, "at": list(at)}
    return result


def main() -> None:
    args = sys.argv[1:]
    if len(args) not in (5, 6) or args[0] != "--json":
        raise SystemExit(__doc__)
    spec = Path(args[2]).read_text().split("\n")
    pins = Path(args[5]) if len(args) == 6 else None
    result = summary(spec, args[3], Path(args[4]), pins)
    Path(args[1]).write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
