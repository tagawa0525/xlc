"""燃焼の履歴（lattice の `xlc-depletion-history`）の 1 つの燃焼点で、燃料の物質の領域の数密度を体積で平均した組成を、
lattice の入力の `nuclides` の形で書く（#0274 の燃料の温度を細かく取ったピンセルの燃えた組成）。

使い方: python3 burned_composition.py <履歴> <燃焼度 (MWd/kg)> <物質の名前の接頭辞> <核種>...

物質の名前が接頭辞で始まる領域の数密度 (1/(barn cm)) を、領域の体積で重みを付けて平均する。燃焼度は履歴の点にあること、
核種は履歴の輸送する核種（`transportable`）にあること。
"""

import sys
from pathlib import Path

import msgpack


def average(history: dict, burnup: float, prefix: str, names: list[str]) -> dict:
    """燃焼度 `burnup` の点で、名前が `prefix` で始まる物質の領域の、核種 `names` の体積で平均した数密度"""
    if history.get("format") != "xlc-depletion-history":
        raise SystemExit("the file is not a depletion history")
    points = [p for p in history["points"] if p["burnup"] == burnup]
    if len(points) != 1:
        raise SystemExit(f"the history has no single point at {burnup} MWd/kg")
    missing = sorted(set(names) - set(history["transportable"]))
    if missing:
        raise SystemExit(f"the history does not transport {missing}")
    columns = [history["nuclides"].index(n) for n in names]
    regions = [
        (r["volume"], d)
        for r, d in zip(history["regions"], points[0]["densities"], strict=True)
        if r["material"].startswith(prefix)
    ]
    if not regions:
        raise SystemExit(f"no region holds a material starting with {prefix!r}")
    volume = sum(v for v, _ in regions)
    return {
        n: sum(v * d[c] for v, d in regions) / volume
        for n, c in zip(names, columns, strict=True)
    }


def main() -> None:
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    path, burnup, prefix, *names = sys.argv[1:]
    history = msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)
    densities = average(history, float(burnup), prefix, names)
    print(
        "nuclides = { "
        + ", ".join(f'"{n}" = {x:.8e}' for n, x in densities.items())
        + " }"
    )


if __name__ == "__main__":
    main()
