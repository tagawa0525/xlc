"""xlc の格子計算（ピンセル 1B〜1D、集合体 2A〜2K、2Q）の k とピン出力を、仕様書の KENO-VI（ENDF/B-VII.0）の値と比べる。

使い方: python3 compare_keno_lattice.py <spec.txt> <問題> <lattice の出力> [<pins.msgpack>]

k は付録の `<問題> <k> <σ>` の行（compare_p3.py の keno_k と同じ）、ピン出力は compare_keno.py と同じ読み方（1/8 の表、
対称な 8 本の平均）。差は (k_xlc − k_KENO-VI) × 1e5 pcm。ピンセルはピン出力を比べない。
"""

import sys
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_k import lattice_k
from compare_keno import compare
from compare_p3 import keno_k


def main() -> None:
    if len(sys.argv) not in (4, 5):
        raise SystemExit(__doc__)
    spec = Path(sys.argv[1]).read_text().split("\n")
    problem = sys.argv[2]
    k = lattice_k(Path(sys.argv[3]))
    k_keno, sigma = keno_k(spec, problem)
    print(
        f"{problem} k xlc {k:.6f} KENO-VI {k_keno:.6f} +- {sigma:.6f} diff {(k - k_keno) * 1e5:+.0f} pcm"
    )
    if len(sys.argv) == 5:
        rows = msgpack.unpackb(Path(sys.argv[4]).read_bytes())["rows"]
        n, rms, worst, at = compare(spec, problem, rows)
        print(
            f"{problem} pins xlc vs KENO-VI: {n} positions, rms {100 * rms:.3f}%, max {100 * worst:+.3f}% at {at}"
        )


if __name__ == "__main__":
    main()
