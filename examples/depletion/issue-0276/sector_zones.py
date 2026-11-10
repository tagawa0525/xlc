"""#0276 の試行: 燃料のピンを方位の扇形の組ごとにも燃焼させる（ADR-0244 の `[depletion]` の `sectors`）効果を測る、2A の燃料の
ピンを並べた最小の格子の入力を書く。

使い方:
  python3 sector_zones.py write   # この置き場に下の 8 本の入力と、その最初の 2 点だけの short-*.toml を書く

格子（どちらも外周は鏡面反射）:
- guide: 中央が 2A の案内管の 3 × 3。案内管の隣のピンは、辺で接するピン（辺のピンの組）と角で接するピン（角のピンの組）の 2 つの
  組になる。外周が鏡面反射なので、案内管が 3 ピッチごとに並ぶ無限の格子（17 × 17 の集合体の案内管と計装管の割合 25/289 より
  密な 1/9）
- gap: 2A の燃料のピンだけの 3 × 3 の周りに、2A の入力と同じ幅 0.04 cm の水の余白を付けたもの。外周が鏡面反射なので、3 × 3 の
  小さな集合体が 0.08 cm の水ギャップを挟んで並ぶ格子で、角のピンは集合体の角（2 辺が水ギャップ）、辺のピンは集合体の辺、中央の
  ピンは水ギャップに接しないピンを表す。案内管と同じ格子に置かないのは、案内管と水ギャップの効きを分けるため

ケース（格子ごとに 4 本）:
- rings1: 本番に向けた 2A の入力（../p9-product/lattice-2A.toml）と同じ燃焼の単位（燃料 1 つ、MOC の 5 リング）。方位に分けない
- rings5: 燃料の 5 リングを 1 本ずつ 5 つの単位にしたもの（#0275 の zones-5-p9 と同じ分け方）。方位に分けない
- rings5-sectors4、rings5-sectors8: rings5 の各単位を、さらに 4 つ（象限）、8 つの扇形の組に分けたもの
MOC の扇形は 4 本とも 8 つ（2A の入力と同じ）で、方位に分けるかだけが違う。径方向の効きは rings5 − rings1、方位の効きは
rings5-sectors* − rings5 で、どちらも同じ格子で比べる。

形状、組成、温度、減速材、燃焼度の点、出力密度、UFG の許容差、トラックは ../p9-product/lattice-2A.toml の燃料のピンと案内管と
公称の状態（ホウ素 650 ppm、減速材 580 K、燃料 900 K）に合わせる。
- short-<名前>.toml: 上の入力の燃焼度を最初の 2 点（輸送 5 回）にしたもの。全部を流す前の時間の見積もり
"""

import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PELLET = 0.4096
RINGS = 5
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
SHORT_BURNUP = "  0.1, 0.25,\n"
LATTICES = {
    # 名前: (行、余白)
    "guide": ('["pin pin pin", "pin gt pin", "pin pin pin"]', ""),
    "gap": (
        '["pin pin pin", "pin pin pin", "pin pin pin"]',
        'gap = { width = 0.04, material = "water" }\n',
    ),
}
CASES = {
    # 名前: (径方向の単位の数、方位の扇形の組の数。None は分けない)
    "rings1": (1, None),
    "rings5": (RINGS, None),
    "rings5-sectors4": (RINGS, 4),
    "rings5-sectors8": (RINGS, 8),
}
GUIDE_TUBE = """\
# 案内管（中は減速材）
[cell.gt]
zones = [
  { material = "water", radius = 0.561, rings = 3 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8

"""


def fuel_names(units: int) -> list[str]:
    return ["fuel"] if units == 1 else [f"fuel-{k}" for k in range(units)]


def input_text(
    lattice: str, units: int, sectors: int | None, output: str, burnup: str = BURNUP
) -> str:
    rows, gap = LATTICES[lattice]
    names = fuel_names(units)
    materials = "".join(
        f"[material.{n}]\ntemperature = 900\nnuclides = {FUEL}\n\n" for n in names
    )
    # 単位は MOC の 5 リングを内から等体積に分けたもの（xlc は帯の中を等体積に分けるので、単位の数によらず MOC の領域が同じ）
    fuel_zones = "".join(
        f'  {{ material = "{n}", radius = {PELLET * math.sqrt((k + 1) / units)!r}, rings = {RINGS // units} }},\n'
        for k, n in enumerate(names)
    )
    listed = ", ".join(f'"{n}"' for n in names)
    split = "" if sectors is None else f"sectors = {sectors}\n"
    guide = GUIDE_TUBE if lattice == "guide" else ""
    radial = "燃料 1 つ" if units == 1 else f"等体積の {units} つの単位"
    azimuthal = (
        "方位に分けない"
        if sectors is None
        else f"各単位を {sectors} つの扇形の組に分ける"
    )
    return f"""\
# #0276 の試行（sector_zones.py が書く。手で直さない）: 格子 {lattice} の 2A の燃料のピンを、{radial}で、{azimuthal}。
# p9-product/lattice-2A.toml の公称の状態のまま燃焼させる。
library = "../p2/depletion.msgpack"
ufg-tolerance = 1e-4
pitch = 1.26
azimuthal = 30
spacing = 0.03
{gap}rows = {rows}

{materials}[material.zircaloy]
temperature = 580

[material.zircaloy.nuclides]
{ZIRCALOY}
[material.water]
temperature = 580
water = {{ density = 0.711901, boron = 650, b10 = 0.199 }}

# 燃料のピン（3.1%）
[cell.pin]
zones = [
{fuel_zones}  {{ material = "void", radius = 0.418, rings = 1 }},
  {{ material = "zircaloy", radius = 0.475, rings = 1 }},
  {{ material = "water", radius = 0.63, rings = 3 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8

{guide}[depletion]
chain = "../chain.msgpack"
power = 38.41
burnup = [
{burnup}]
materials = [{listed}]
{split}output = "{output}"
"""


def write() -> None:
    for lattice in LATTICES:
        for case, (units, sectors) in CASES.items():
            name = f"{lattice}-{case}"
            text = input_text(lattice, units, sectors, f"{name}.msgpack")
            (HERE / f"{name}.toml").write_text(text)
            short = input_text(
                lattice, units, sectors, f"short-{name}.msgpack", SHORT_BURNUP
            )
            (HERE / f"short-{name}.toml").write_text(short)


def main(argv: list[str]) -> None:
    if argv == ["write"]:
        write()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
