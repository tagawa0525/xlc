"""試行: Problem 3B の Pyrex の端栓の端の境界（15.761 cm）を残した解と除いた解の差（k と、付録 D の 49 層の軸方向の出力）。

使い方: python3 compare.py <spec.txt> <境界ありの標準出力> <境界ありの layers.msgpack> <除いた標準出力> <除いた layers.msgpack>
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "examples/vera"))
from compare_p3 import coarsen, keno_axial, offset, read_layers


def k_eff(path: str) -> float:
    lines = [l for l in Path(path).read_text().splitlines() if l.startswith("k-eff ")]
    if len(lines) != 1:
        raise SystemExit(f"{path}: expected one k-eff line")
    return float(lines[0].split()[1])


def main() -> None:
    if len(sys.argv) != 6:
        raise SystemExit(__doc__)
    spec = Path(sys.argv[1]).read_text().split("\n")
    bounds, _ = keno_axial(spec, "3B")
    k = [k_eff(sys.argv[2]), k_eff(sys.argv[4])]
    axial = [
        coarsen(read_layers(Path(p))["bounds"], read_layers(Path(p))["axial"], bounds)
        for p in (sys.argv[3], sys.argv[5])
    ]
    print(
        f"k with the bound {k[0]:.8f}, merged {k[1]:.8f}, diff {(k[1] / k[0] - 1) * 1e5:+.2f} pcm"
    )
    rel = [b / a - 1 for a, b in zip(*axial)]
    worst = max(range(len(rel)), key=lambda i: abs(rel[i]))
    print(
        f"axial (49 levels) relative diff: max {100 * rel[worst]:+.4f}% at level {worst + 1} "
        f"({bounds[worst]:.3f}-{bounds[worst + 1]:.3f} cm above the fuel bottom), "
        f"AO {100 * offset(bounds, axial[0]):+.4f}% -> {100 * offset(bounds, axial[1]):+.4f}%"
    )


if __name__ == "__main__":
    main()
