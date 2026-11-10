"""#0133 の試行: 2A の燃料のピンのピンセルを、ペレットを径方向に等体積の N 個の燃焼の単位（物質）に分けて燃焼させる入力を書き、
単位の数ごとの履歴を比べる。

使い方:
  python3 pellet_zones.py write              # この置き場に下の 5 本の入力と、その最初の 2 点だけの short-*.toml を書く
  python3 pellet_zones.py compare <基準の履歴> <履歴>...   # 基準（1 つの単位）に対する k∞ と数密度の差、リングの Pu-239 の分布

形状、組成、温度、減速材、燃焼度の点、出力密度、UFG の許容差は ../../p9-product/lattice-2A.toml の燃料のピンと公称の状態
（ホウ素 650 ppm、減速材 580 K、燃料 900 K）に合わせる。
- zones-1、zones-4、zones-8: 燃料の中性子束のリングは単位の数によらず等体積の 8 本で、単位の境界は 8 本の境界に重なる（xlc は
  ゾーンの中を等体積に分けるので、MOC の領域が同じになる）。単位の数だけを変えるため、3 本とも飛程を `track-table = false` で
  格子点ごとに求める（8 つの物質で飛程の表を使えるかは確かめていない）
- zones-1-p9、zones-5-p9: 本番の離散化（燃料 5 リング、飛程の表）のまま、1 つの単位と、5 リングを 1 本ずつ 5 つの単位にしたもの。
  本番のテーブルで分けたときの時間の比と k∞ の差を、同じ離散化で見る
- short-<名前>.toml: 上の入力の燃焼度を最初の 2 点（輸送 5 回）にしたもの。全部を流す前の時間の見積もり
"""

import math
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
PELLET = 0.4096
FLUX_RINGS = 8
FUEL = "{ U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02 }"
ZIRCALOY = """\
Cr-50 = 3.30121e-06
Cr-52 = 6.36606e-05
Cr-53 = 7.21860e-06
Cr-54 = 1.79686e-06
Fe-54 = 8.68307e-06
Fe-56 = 1.36306e-04
Fe-57 = 3.14789e-06
Fe-58 = 4.18926e-07
Zr-90 = 2.18865e-02
Zr-91 = 4.77292e-03
Zr-92 = 7.29551e-03
Zr-94 = 7.39335e-03
Zr-96 = 1.19110e-03
Sn-112 = 4.68066e-06
Sn-114 = 3.18478e-06
Sn-115 = 1.64064e-06
Sn-116 = 7.01616e-05
Sn-117 = 3.70592e-05
Sn-118 = 1.16872e-04
Sn-119 = 4.14504e-05
Sn-120 = 1.57212e-04
Sn-122 = 2.23417e-05
Sn-124 = 2.79392e-05
"""
BURNUP = """\
  0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5,
  4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
  14, 15, 16, 18, 20, 22, 24, 26, 28, 30,
  32, 34, 36, 38, 40, 44, 48, 52, 56, 60,
"""
NUCLIDES = ["U-235", "U-238", "Pu-239", "Pu-240", "Pu-241"]


def fuel_names(zones: int) -> list[str]:
    return ["fuel"] if zones == 1 else [f"fuel-{k}" for k in range(zones)]


def input_text(
    zones: int,
    rings_per_zone: int,
    track_table: bool,
    output: str,
    burnup: str = BURNUP,
) -> str:
    names = fuel_names(zones)
    materials = "".join(
        f"[material.{n}]\ntemperature = 900\nnuclides = {FUEL}\n\n" for n in names
    )
    fuel_zones = "".join(
        f'  {{ material = "{n}", radius = {PELLET * math.sqrt((k + 1) / zones)!r}, rings = {rings_per_zone} }},\n'
        for k, n in enumerate(names)
    )
    table = "" if track_table else "track-table = false\n"
    listed = ", ".join(f'"{n}"' for n in names)
    return f"""\
# #0133 の試行（pellet_zones.py が書く。手で直さない）: 2A の燃料のピンのピンセルを、ペレットを等体積の {zones} 個の燃焼の単位に
# 分け、燃料の中性子束のリング {zones * rings_per_zone} 本、飛程の{"表" if track_table else "格子点ごとの計算"}で、p9-product/lattice-2A.toml の
# 公称の状態のまま燃焼させる。
library = "../../p2/depletion.msgpack"
{table}ufg-tolerance = 1e-4
pitch = 1.26
azimuthal = 30
spacing = 0.03
rows = ["pin"]

{materials}[material.zircaloy]
temperature = 580

[material.zircaloy.nuclides]
{ZIRCALOY}
[material.water]
temperature = 580
water = {{ density = 0.711901, boron = 650, b10 = 0.199 }}

[cell.pin]
zones = [
{fuel_zones}  {{ material = "void", radius = 0.418, rings = 1 }},
  {{ material = "zircaloy", radius = 0.475, rings = 1 }},
  {{ material = "water", radius = 0.63, rings = 3 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8

[depletion]
chain = "../../chain.msgpack"
power = 38.41
burnup = [
{burnup}]
materials = [{listed}]
output = "{output}"
"""


CASES = {
    # 名前: (単位の数, 単位ごとの MOC のリングの数, 飛程の表)
    "zones-1": (1, FLUX_RINGS, False),
    "zones-4": (4, FLUX_RINGS // 4, False),
    "zones-8": (8, FLUX_RINGS // 8, False),
    "zones-1-p9": (1, 5, True),
    "zones-5-p9": (5, 1, True),
}
SHORT_BURNUP = "  0.1, 0.25,\n"


def write() -> None:
    for name, (zones, rings, table) in CASES.items():
        text = input_text(zones, rings, table, f"{name}.msgpack")
        (HERE / f"{name}.toml").write_text(text)
        short = input_text(zones, rings, table, f"short-{name}.msgpack", SHORT_BURNUP)
        (HERE / f"short-{name}.toml").write_text(short)


def is_heavy_metal(name: str) -> bool:
    symbol = name.split("-")[0]
    return symbol in {"Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es"}


def summary(path: Path) -> dict:
    h = msgpack.unpackb(path.read_bytes())
    if h["format"] != "xlc-depletion-history":
        raise SystemExit(f"{path}: not an xlc depletion history")
    regions = h["regions"]
    volumes = [r["volume"] for r in regions]
    total = sum(volumes)
    index = {n: i for i, n in enumerate(h["nuclides"])}
    heavy = [i for n, i in index.items() if is_heavy_metal(n)]
    points = []
    for p in h["points"]:
        d = p["densities"]
        mean = {
            n: sum(d[r][index[n]] * volumes[r] for r in range(len(regions))) / total
            for n in NUCLIDES
        }
        hm = sum(d[r][i] * volumes[r] for r in range(len(regions)) for i in heavy)
        pu = [d[r][index["Pu-239"]] for r in range(len(regions))]
        points.append(
            {
                "burnup": p["burnup"],
                "k": p["k_inf"],
                "mean": mean,
                "hm": hm,
                "pu239_outer_inner": pu[-1] / pu[0] if pu[0] > 0 else None,
            }
        )
    return {
        "path": str(path),
        "regions": len(regions),
        "volume": total,
        "points": points,
        "seconds": h["seconds"],
    }


def compare(base: Path, others: list[Path]) -> None:
    b = summary(base)
    rows = [summary(o) for o in others]
    print(f"base {b['path']}: {b['regions']} regions, volume {b['volume']:.9g} cm2")
    for s in rows:
        print(f"{s['path']}: {s['regions']} regions, volume {s['volume']:.9g} cm2")
        if len(s["points"]) != len(b["points"]):
            raise SystemExit(f"{s['path']}: burnup points differ from the base")
    print("seconds:", {Path(s["path"]).name: s["seconds"] for s in [b, *rows]})
    header = ["burnup", "k(base)"]
    for s in rows:
        tag = Path(s["path"]).stem
        header += [f"dk {tag} (pcm)", f"dHM {tag}"]
        header += [f"d{n} {tag} (%)" for n in NUCLIDES]
        header += [f"Pu239 out/in {tag}"]
    print(" | ".join(header))
    for i, p in enumerate(b["points"]):
        cells = [f"{p['burnup']:g}", f"{p['k']:.6f}"]
        for s in rows:
            q = s["points"][i]
            if abs(q["burnup"] - p["burnup"]) > 1e-9 * max(1.0, p["burnup"]):
                raise SystemExit(f"{s['path']}: burnup {q['burnup']} != {p['burnup']}")
            cells.append(f"{(q['k'] - p['k']) * 1e5:+.1f}")
            cells.append(f"{q['hm'] / p['hm'] - 1:+.1e}")
            for n in NUCLIDES:
                base_n = p["mean"][n]
                cells.append(
                    "–" if base_n == 0 else f"{(q['mean'][n] / base_n - 1) * 100:+.3f}"
                )
            r = q["pu239_outer_inner"]
            cells.append("–" if r is None else f"{r:.3f}")
        print(" | ".join(cells))


def main() -> None:
    match sys.argv[1:]:
        case ["write"]:
            write()
        case ["compare", base, *others] if others:
            compare(Path(base), [Path(o) for o in others])
        case _:
            raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
