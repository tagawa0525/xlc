"""#0293 の項目 2: 燃焼度の点の数だけが違う 2 本の燃焼の履歴（steps/out/<燃料>-40.msgpack と -80.msgpack）を、粗い方の点で比べる。

燃焼点ごとに、k∞ の差（細かい − 粗い、pcm）と、主な核種の数密度の相対差（細かい / 粗い − 1）を書く。数密度は、物質の名前の
接頭辞（`fuel` か `gd`）ごとに燃焼の領域の体積で平均する。CE/CM（ADR-0173）は 2 次なので、粗い方の誤差は差のおよそ 4/3 倍
（細かい方の誤差は差の 1/3）。

使い方: python3 compare_steps.py <燃料> [--coarse 40] [--fine 80] [--json <出力>]（steps/out/<燃料>-<点の数>.msgpack を読む）
"""

import argparse
import json
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent

NUCLIDES = ["U-235", "Pu-239", "Pu-240", "Pu-241", "Xe-135", "Sm-149"]
GD_NUCLIDES = ["Gd-155", "Gd-157"]


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def averaged(history: dict, point: dict, prefix: str, nuclide: str) -> float:
    """燃焼点 `point` の、名前が `prefix` で始まる物質の燃焼の領域の体積で平均した `nuclide` の数密度"""
    index = history["nuclides"].index(nuclide)
    pairs = [
        (region["volume"], densities[index])
        for region, densities in zip(history["regions"], point["densities"])
        if region["material"].startswith(prefix)
    ]
    if not pairs:
        sys.exit(f"no depletion region of a material starting with `{prefix}`")
    volume = sum(v for v, _ in pairs)
    return sum(v * d for v, d in pairs) / volume


def compare(coarse: dict, fine: dict, groups: dict[str, list[str]]) -> list[dict]:
    """粗い履歴の燃焼点ごとの、k∞ の差 (pcm) と、物質の接頭辞ごと核種ごとの数密度の相対差"""
    by_burnup = {p["burnup"]: p for p in fine["points"]}
    rows = []
    for point in coarse["points"]:
        other = by_burnup.get(point["burnup"])
        if other is None:
            sys.exit(f"the fine history has no point at {point['burnup']} MWd/kg")
        densities = {}
        for prefix, nuclides in groups.items():
            for nuclide in nuclides:
                low = averaged(coarse, point, prefix, nuclide)
                high = averaged(fine, other, prefix, nuclide)
                densities[f"{prefix}:{nuclide}"] = high / low - 1 if low > 0 else None
        rows.append(
            {
                "burnup": point["burnup"],
                "k_coarse": point["k_inf"],
                "dk_pcm": (other["k_inf"] - point["k_inf"]) * 1e5,
                "densities": densities,
            }
        )
    return rows


def groups_of(fuel: str) -> dict[str, list[str]]:
    if fuel == "gd":
        return {"fuel": NUCLIDES, "gd": NUCLIDES + GD_NUCLIDES}
    return {"fuel": NUCLIDES}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fuel", choices=["u", "mox", "gd"])
    parser.add_argument("--coarse", type=int, default=40, help="粗い方の燃焼度の点の数")
    parser.add_argument("--fine", type=int, default=80, help="細かい方の燃焼度の点の数")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    out = HERE / "steps" / "out"
    rows = compare(
        read(out / f"{args.fuel}-{args.coarse}.msgpack"),
        read(out / f"{args.fuel}-{args.fine}.msgpack"),
        groups_of(args.fuel),
    )
    keys = list(rows[0]["densities"])
    print(f"burnup k{args.coarse} dk(pcm) " + " ".join(keys))
    for r in rows:
        values = " ".join(
            "-" if r["densities"][k] is None else f"{r['densities'][k] * 100:+.3f}%"
            for k in keys
        )
        print(f"{r['burnup']:g} {r['k_coarse']:.5f} {r['dk_pcm']:+.1f} {values}")
    if args.json:
        args.json.write_text(
            json.dumps(
                {
                    "fuel": args.fuel,
                    "coarse": args.coarse,
                    "fine": args.fine,
                    "rows": rows,
                },
                indent=1,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
