"""VERA のピンセル（1B〜1D）の燃料・被覆管・減速材を体積比で混ぜた均質な無限媒質の入力を書く（xlc の `homogeneous` の例
（crates/lattice/examples/homogeneous.rs）と `openmc_homogeneous.py` が読む）。

使い方: python3 homogeneous_medium.py <問題> <出力 (TOML)> [残す核種 ...]

個数密度は、ピッチの正方形の中の体積（ギャップは真空）で平均したもの。残す核種（xlc の名前、例 U-235）を
指定すると、それ以外を除く（残す核種の個数密度は変えない）。温度は燃料、被覆管、減速材が同じ問題（1B）だけ。
"""

import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from openmc_model import (
    CLAD_INNER,
    CLAD_OUTER,
    COLD,
    FUEL,
    FUEL_TEMPERATURE,
    MODERATOR,
    PELLET,
    PITCH,
    ZIRCALOY,
)


def xlc_name(name: str) -> str:
    """OpenMC の核種の名前（U235）を xlc の名前（U-235）にする。"""
    match = re.fullmatch(r"([A-Z][a-z]?)(\d+)", name)
    if match is None:
        raise SystemExit(f"unexpected nuclide name {name}")
    return f"{match[1]}-{match[2]}"


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    problem, output = sys.argv[1], Path(sys.argv[2])
    keep = set(sys.argv[3:])
    if FUEL_TEMPERATURE.get(problem) != COLD or not problem.startswith("1"):
        raise SystemExit(f"{problem}: only a pin cell with one temperature (1B)")
    volumes = {
        "fuel": math.pi * PELLET**2,
        "clad": math.pi * (CLAD_OUTER**2 - CLAD_INNER**2),
        "water": PITCH**2 - math.pi * CLAD_OUTER**2,
    }
    total = PITCH**2
    densities: dict[str, float] = {}
    for part, composition in (("fuel", FUEL), ("clad", ZIRCALOY), ("water", MODERATOR)):
        for name, n in composition.items():
            name = xlc_name(name)
            densities[name] = densities.get(name, 0.0) + n * volumes[part] / total
    unknown = keep - set(densities)
    if unknown:
        raise SystemExit(f"not in the medium: {sorted(unknown)}")
    if keep:
        densities = {k: v for k, v in densities.items() if k in keep}
    lines = [
        f"# VERA {problem} の均質な無限媒質"
        + (f"（{', '.join(sorted(keep))} だけ）" if keep else "")
        + "。homogeneous_medium.py で作った",
        f"temperature = {COLD}",
        "[nuclides]",
    ] + [f'"{k}" = {v:.6e}' for k, v in densities.items()]
    output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
