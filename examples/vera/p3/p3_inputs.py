"""VERA Problem 3A と 3B（3D の単一集合体、零出力）を xlc で解く入力を書く。

このファイルのある場所に、問題（3A、3B）ごとに次を書く。寸法と組成の出典と導出は、下の定数のコメントにあり、書いた
入力のコメントにもコピーする。

- `lattice-<問題>-<層>.toml`: 燃料の高さの 2D の単一集合体の格子（外周は鏡面反射）。層は支持格子と挿入物（Pyrex の棒の
  毒物、下の端栓、プレナム）の組み合わせごとで、少数群定数（形状関数つき）を `lattice-<問題>/<層>.msgpack` に書く。
  支持格子のない層は、軸方向の反射体の平板の燃料にする媒質（`medium-<問題>-<層>.msgpack`）も書く（ADR-0128）
- `axial-<問題>.toml`: 集合体の高さ全体（下の水の外端から上の水の外端まで）を格子の x（左が炉心の −z）に並べた1行の
  平板（ADR-0128）。燃料は上の媒質、層の物質は層の断面の中で体積で混ぜる。燃料のノードは炉心の燃料のノードと同じ
  高さに分ける（f* が燃料の端のノードの高さに依存するため。ADR-0128 の確認）
- `core-<問題>.toml` ほか（ADR-0105 の4つ）: 燃料の層は 2D の格子の定数、軸方向の反射体の層は平板の定数（f* の方法）
- `core-<問題>-slab.toml` ほか: すべての層を平板の定数にした、定数を作った体系の等価性の確認

炉心の z は下の水の外端から（仕様書の高さ（下部炉心板の上面から）+ 55 cm）。燃料の中の炉心のノードは、仕様書の
表 P3-5 の 49 層の境界に、挿入物の端（3B の Pyrex の棒の下の端栓、毒物の下端と上端）を足したもの。

使い方: python3 p3_inputs.py
"""

import itertools
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent

# 仕様書の高さ（下部炉心板の上面から、cm）を炉心の z にずらす量: 下の水 50 cm と下部炉心板 5 cm
OFFSET = Fraction("55")

# 表 P3-5 の 49 層の上端（仕様書の高さ、cm）。下端は燃料の下端 11.951 cm
LEVEL_TOPS = [
    "15.817", "24.0281", "32.2393", "40.4504", "48.6616", "56.8727", "65.0839", "73.295", "77.105",
    "85.17", "93.235", "101.3", "109.365", "117.43", "125.495", "129.305", "137.37", "145.435",
    "153.5", "161.565", "169.63", "177.695", "181.505", "189.57", "197.635", "205.7", "213.765",
    "221.83", "229.895", "233.705", "241.77", "249.835", "257.9", "265.965", "274.03", "282.095",
    "285.905", "293.97", "302.035", "310.1", "318.165", "326.23", "334.295", "338.105", "346.0262",
    "353.9474", "361.8686", "369.7898", "377.711",
]  # fmt: skip
FUEL_BOTTOM = Fraction("11.951")
FUEL_TOP = Fraction("377.711")

# 支持格子（表 3）。1 つの質量を 289 のセルに等しく分け、セルの外周の正方形の帯の枠にする（2Q と同じ）。厚さ t は面積
# A = 質量 / 密度 / 高さ / 289 を 4pt − 4t² = A で解いたもの（openmc_model_3d.py の frame_thickness と同じ）。
# 端の 2 つは Inconel-718（1017 g、8.19 g/cm³、3.866 cm、中心 13.884 と 388.2 cm）、中間の 6 つは Zircaloy-4（875 g、
# 6.56 g/cm³、3.810 cm、中心 75.2〜336.2 cm）。燃料の中の支持格子は表 P3-5 の層とちょうど重なる
END_GRID = ("inconel", 1017.0, 8.19, Fraction("3.866"), ["13.884"])
MID_GRID = (
    "zircaloy",
    875.0,
    6.56,
    Fraction("3.81"),
    ["75.2", "127.4", "179.6", "231.8", "284.0", "336.2"],
)
TOP_GRID_CENTRE = Fraction("388.2")
PITCH = 1.26


def frame_thickness(mass: float, density: float, height: Fraction) -> float:
    area = mass / density / float(height) / 289
    return (PITCH - math.sqrt(PITCH * PITCH - area)) / 2


# 3B の Pyrex の棒（表 6）と、シンブルプラグ（表 9）。解釈は openmc_model_3d.py と同じ: 端栓の高さ 2.54 cm は下の端栓
# （SS304 の中実の棒、半径は被覆管の外側）、プレナム（被覆管の中は真空。内管は寸法がないので置かない）は上部ノズルの
# 下端で切る。シンブルプラグは SS304 の中実の棒（半径 0.538 cm）で、383.31 cm を下端とする
POISON_BOTTOM = Fraction("15.761")
POISON_TOP = POISON_BOTTOM + Fraction("360.68")
PLUG_BOTTOM = POISON_BOTTOM - Fraction("2.54")
THIMBLE_BOTTOM = Fraction("383.31")
THIMBLE_TOP = THIMBLE_BOTTOM + Fraction("11.0")

# 軸方向の構造（仕様書の高さ、cm。表 1、4、5）
NOZZLE_BOTTOM_TOP = Fraction("6.053")
ROD_BOTTOM = Fraction("10.281")
PLENUM_TOP = FUEL_TOP + Fraction("16.0")
ROD_TOP = PLENUM_TOP + Fraction("1.67")
NOZZLE_TOP_BOTTOM = Fraction("397.51")
NOZZLE_TOP_TOP = Fraction("406.337")
PLATE_TOP = NOZZLE_TOP_TOP + Fraction("7.6")

# 列の幅の目安 (cm)。格子計算の CMFD の収束で 1 cm 程度に制限する（#0069。examples/reflector/axial_slab.py と同じ）
WIDTH = 1.0
# 平板の格子の行の高さ (cm)。等価性の炉心はピッチを同じ値にして、ノード平均の中性子束の規格化を平板とそろえる
ROW = "1.0"

# 問題ごとの温度 (K)、燃料と冷却材の組成（表 P3-1、P3-3）
PROBLEMS = {
    "3A": {
        "temperature": 600,
        "fuel": "U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02",
        "water": "H-1 = 4.96224e-02, O-16 = 2.48112e-02, B-10 = 1.07070e-05, B-11 = 4.30971e-05",
        "description": "濃縮度 3.10%、ホウ素 1300 ppm、600 K",
    },
    "3B": {
        "temperature": 565,
        "fuel": "U-234 = 5.09503e-06, U-235 = 6.06733e-04, U-236 = 2.76809e-06, U-238 = 2.22663e-02, O-16 = 4.57617e-02",
        "water": "H-1 = 4.96340e-02, O-16 = 2.48170e-02, B-10 = 8.77976e-06, B-11 = 3.53397e-05",
        "description": "濃縮度 2.619%、ホウ素 1066 ppm、565 K、Pyrex の棒 16 本とシンブルプラグ 8 本",
    },
}

ZIRCALOY = """Cr-50 = 3.30121e-06
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
Sn-124 = 2.79392e-05"""

SS304 = """C-12 = 3.17461e-04
C-13 = 3.43358e-06
Si-28 = 1.58197e-03
Si-29 = 8.03653e-05
Si-30 = 5.30394e-05
P-31 = 6.99938e-05
Cr-50 = 7.64915e-04
Cr-52 = 1.47506e-02
Cr-53 = 1.67260e-03
Cr-54 = 4.16346e-04
Mn-55 = 1.75387e-03
Fe-54 = 3.44776e-03
Fe-56 = 5.41225e-02
Fe-57 = 1.24992e-03
Fe-58 = 1.66342e-04
Ni-58 = 5.30854e-03
Ni-60 = 2.04484e-03
Ni-61 = 8.88879e-05
Ni-62 = 2.83413e-04
Ni-64 = 7.21770e-05"""

INCONEL = """Si-28 = 4.04885e-03
Si-29 = 2.05685e-04
Si-30 = 1.35748e-04
Ti-46 = 2.12518e-04
Ti-47 = 1.91652e-04
Ti-48 = 1.89901e-03
Ti-49 = 1.39360e-04
Ti-50 = 1.33435e-04
Cr-50 = 6.18222e-04
Cr-52 = 1.19218e-02
Cr-53 = 1.35184e-03
Cr-54 = 3.36501e-04
Fe-54 = 3.61353e-04
Fe-56 = 5.67247e-03
Fe-57 = 1.31002e-04
Fe-58 = 1.74340e-05
Ni-58 = 4.17608e-02
Ni-60 = 1.60862e-02
Ni-61 = 6.99255e-04
Ni-62 = 2.22953e-03
Ni-64 = 5.67796e-04"""

PYREX = "B-10 = 9.63266e-04, B-11 = 3.90172e-03, O-16 = 4.67761e-02, Si-28 = 1.81980e-02, Si-29 = 9.24474e-04, Si-30 = 6.10133e-04"


def materials(problem: str, names: list[str]) -> str:
    """物質の TOML（`names` の順。燃料、冷却材、Zircaloy-4、SS304、Inconel-718、Pyrex。Hf は除く。ADR-0054）"""
    p = PROBLEMS[problem]
    t = p["temperature"]
    blocks = {
        "fuel": f"[material.fuel]\ntemperature = {t}\nnuclides = {{ {p['fuel']} }}",
        "water": f"[material.water]\ntemperature = {t}\nnuclides = {{ {p['water']} }}",
        "zircaloy": f"[material.zircaloy]\ntemperature = {t}\n\n[material.zircaloy.nuclides]\n{ZIRCALOY}",
        "ss304": f"[material.ss304]\ntemperature = {t}\n\n[material.ss304.nuclides]\n{SS304}",
        "inconel": f"[material.inconel]\ntemperature = {t}\n\n[material.inconel.nuclides]\n{INCONEL}",
        "pyrex": f"[material.pyrex]\ntemperature = {t}\nnuclides = {{ {PYREX} }}",
    }
    return "\n\n".join(blocks[n] for n in names)


# 17×17 の配置（行は上から）。案内管 24 本の位置（仕様書の図 1。openmc_model.py の GUIDE_TUBES）のうち、3B で Pyrex の棒を
# 入れる 16 本（図 5 の 1/8 の中の、中央からの (3, 0)、(5, 5)、(6, 3)）を `P`、シンブルプラグを入れる 8 本（(3, 3)、
# (6, 0)）を `T`、中央の計装管を `I` とする（openmc_model.py の PYREX_OCTANTS と octant_images）
GUIDE_TUBES = [
    (2, 5), (2, 8), (2, 11), (3, 3), (3, 13), (5, 2), (5, 5), (5, 8), (5, 11), (5, 14), (8, 2), (8, 5),
    (8, 11), (8, 14), (11, 2), (11, 5), (11, 8), (11, 11), (11, 14), (13, 3), (13, 13), (14, 5), (14, 8), (14, 11),
]  # fmt: skip
PYREX_OCTANT = {(3, 0), (5, 5), (6, 3)}


def octant_images(octant: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """1/8 の中の位置（中央からの行と列の数）を、8 つの対称操作で移した位置（行, 列。0 始まり）"""
    return {
        (8 + si * x, 8 + sj * y)
        for a, b in octant
        for x, y in ((a, b), (b, a))
        for si in (1, -1)
        for sj in (1, -1)
    }


def layout() -> list[str]:
    pyrex = octant_images(PYREX_OCTANT)
    if not pyrex <= set(GUIDE_TUBES) or len(pyrex) != 16:
        raise SystemExit("the Pyrex rods are not 16 guide tubes")
    rows = []
    for i in range(17):
        row = ""
        for j in range(17):
            if (i, j) == (8, 8):
                row += "I"
            elif (i, j) in pyrex:
                row += "P"
            elif (i, j) in GUIDE_TUBES:
                row += "T"
            else:
                row += "."
        rows.append(row)
    return rows


LAYOUT = layout()

# 2D の格子のセル（帯は内側から。examples/vera/p2 の 2A、2F と同じ分け方）
CELLS = {
    "pin": """zones = [
  { material = "fuel", radius = 0.4096, rings = 5 },
  { material = "void", radius = 0.418, rings = 1 },
  { material = "zircaloy", radius = 0.475, rings = 1 },
  { material = "water", radius = 0.63, rings = 3 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
    # 案内管（中は減速材）
    "gt": """zones = [
  { material = "water", radius = 0.561, rings = 3 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
    # 計装管（中は減速材）
    "it": """zones = [
  { material = "water", radius = 0.559, rings = 3 },
  { material = "zircaloy", radius = 0.605, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
    # Pyrex の棒の毒物の高さ（内管の中と、Pyrex の内外の隙間は真空。表 6）
    "py": """zones = [
  { material = "void", radius = 0.214, rings = 1 },
  { material = "ss304", radius = 0.231, rings = 1 },
  { material = "void", radius = 0.241, rings = 1 },
  { material = "pyrex", radius = 0.427, rings = 2 },
  { material = "void", radius = 0.437, rings = 1 },
  { material = "ss304", radius = 0.484, rings = 1 },
  { material = "water", radius = 0.561, rings = 1 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
    # Pyrex の棒の下の端栓（SS304 の中実の棒、半径は被覆管の外側）
    "pl": """zones = [
  { material = "ss304", radius = 0.484, rings = 2 },
  { material = "water", radius = 0.561, rings = 1 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
    # Pyrex の棒のプレナム（被覆管の中は真空）
    "pp": """zones = [
  { material = "void", radius = 0.437, rings = 1 },
  { material = "ss304", radius = 0.484, rings = 1 },
  { material = "water", radius = 0.561, rings = 1 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]""",
}

# 2D の格子の層: (名前, Pyrex の位置のセル, 支持格子, 説明)。3A は Pyrex の位置も案内管
GRIDS = {
    None: None,
    "zircaloy": frame_thickness(*MID_GRID[1:4]),
    "inconel": frame_thickness(*END_GRID[1:4]),
}
LATTICES = {
    "3A": [
        ("fuel", "gt", None, "支持格子のない燃料"),
        ("zgrid", "gt", "zircaloy", "中間の支持格子（Zircaloy-4）の高さの燃料"),
        ("igrid", "gt", "inconel", "下の端の支持格子（Inconel-718）の高さの燃料"),
    ],
    "3B": [
        ("poison", "py", None, "Pyrex の棒の毒物の高さの燃料"),
        (
            "poison-zgrid",
            "py",
            "zircaloy",
            "毒物の高さの、中間の支持格子（Zircaloy-4）の高さの燃料",
        ),
        (
            "empty",
            "gt",
            None,
            "Pyrex の棒の下の端栓より下の燃料（案内管は空。平板の媒質だけに使う）",
        ),
        (
            "plug",
            "pl",
            None,
            "Pyrex の棒の下の端栓の高さの燃料（平板の媒質だけに使う）",
        ),
        ("plenum", "pp", None, "Pyrex の棒のプレナムの高さの燃料"),
        (
            "empty-igrid",
            "gt",
            "inconel",
            "下の端の支持格子（Inconel-718）の高さの、端栓より下の燃料",
        ),
        (
            "plug-igrid",
            "pl",
            "inconel",
            "下の端の支持格子の高さの、Pyrex の棒の下の端栓の高さの燃料",
        ),
        ("poison-igrid", "py", "inconel", "下の端の支持格子の高さの、毒物の高さの燃料"),
    ],
}


def lattice(problem: str, name: str, insert: str, grid: str | None, what: str) -> str:
    cell_of = {".": "pin", "I": "it", "T": "gt", "P": insert}
    rows = ",\n".join(
        '  "' + " ".join(f"{cell_of[c]:<3}" for c in row).rstrip() + '"'
        for row in LAYOUT
    )
    used = sorted({cell_of[c] for row in LAYOUT for c in row}, key=list(CELLS).index)
    names = ["fuel", "water", "zircaloy"]
    if insert != "gt":
        names.append("ss304")
    if insert == "py":
        names.append("pyrex")
    if grid == "inconel":
        names.append("inconel")
    frame = (
        f'\nframe = {{ thickness = {GRIDS[grid]:.6f}, material = "{grid}" }}'
        if grid
        else ""
    )
    cells = "\n\n".join(f"[cell.{c}]\n{CELLS[c]}\nsectors = 8{frame}" for c in used)
    medium = f'medium = "medium-{problem}-{name}.msgpack"\n' if grid is None else ""
    grid_note = ""
    if grid:
        g = END_GRID if grid == "inconel" else MID_GRID
        grid_note = (
            f"\n# 支持格子は {g[1]:g} g / {g[2]} g/cm³ / {g[3]} cm を 289 のセルに等しく分け、セルの外周の厚さ "
            f"{GRIDS[grid]:.6f} cm の枠にする（表 3。2Q と同じ置き方）"
        )
    return f"""# このファイルは p3_inputs.py が書く。手で直さない。
# VERA Problem {problem}（{PROBLEMS[problem]["description"]}。表 P3-1、P3-3）の、{what}の単一集合体の格子
# （外周は鏡面反射）。セルの分け方とトラックは examples/vera/p2 の集合体と同じ（ADR-0111）。Hf とギャップの He は除き
# （ADR-0054）、Pyrex の棒の中の He も真空（ADR-0120）。{grid_note}
library = "../vera.msgpack"
few-group = {{ output = "lattice-{problem}/{name}.msgpack", bounds = [1e7, 0.625, 1e-5], pins = true }}
{medium}pitch = 1.26
azimuthal = 30
spacing = 0.03
gap = {{ width = 0.04, material = "water" }}
rows = [
{rows},
]

{materials(problem, names)}

{cells}
"""


def decimal(x: Fraction) -> str:
    """有限小数で書ける分数を、丸めずに 10 進で書く"""
    for digits in range(14):
        if (x * 10**digits).denominator == 1:
            return f"{float(x):.{max(digits, 1)}f}"
    raise ValueError(x)


def columns(t: Fraction) -> tuple[int, Fraction]:
    """厚さを、幅 `WIDTH` に近く小数 7 桁で割り切れる幅の列に分ける"""
    n0 = max(1, round(float(t) / WIDTH))
    for d in range(60):
        for n in (n0 + d, n0 - d):
            if n >= 1 and (t / n * 10**7).denominator == 1:
                return n, t / n
    raise ValueError(t)


def areas(problem: str) -> dict[str, float]:
    """集合体の断面（ピッチ 21.50 cm の正方形）の中の、層の物質の面積 (cm²)。仕様書の表 1、2、3、6、9"""
    gt = 24 * math.pi * (0.602**2 - 0.561**2)
    it = math.pi * (0.605**2 - 0.559**2)
    rods = 264
    out = {
        "assembly": 21.50 * 21.50,
        "tubes": gt + it,
        # 端栓は被覆管の外半径の Zircaloy の棒（端栓の半径は仕様書にない。ADR-0128）
        "plugs": rods * math.pi * 0.475**2,
        "clad": rods * math.pi * (0.475**2 - 0.418**2),
        # プレナムの He は真空（ADR-0120）。ばねは仕様書に質量がなく、入れない
        "plenum gas": rods * math.pi * 0.418**2,
        # 上の端の支持格子の Inconel-718: 質量 / 密度 / 高さ
        "grid": END_GRID[1] / END_GRID[2] / float(END_GRID[3]),
        "pyrex clad": 0.0,
        "pyrex gas": 0.0,
        "thimble": 0.0,
    }
    if problem == "3B":
        out["pyrex clad"] = 16 * math.pi * (0.484**2 - 0.437**2)
        out["pyrex gas"] = 16 * math.pi * 0.437**2
        out["thimble"] = 8 * math.pi * 0.538**2
    return out


def mixture(problem: str, parts: dict[str, float]) -> list[tuple[str, str]]:
    """面積 (cm²) の成分を体積の割合にし、残りを水にする（割合の和がちょうど 1 になるよう、水は 1 − ほかの和）"""
    total = areas(problem)["assembly"]
    written = [(m, f"{v / total:.8f}") for m, v in parts.items() if v > 0]
    rest = Fraction(1) - sum(Fraction(f) for _, f in written)
    return written + [("water", f"{float(rest):.8f}")]


def structure(problem: str, top: Fraction, bottom: Fraction) -> tuple[str, list]:
    """燃料の外の、仕様書の高さ [bottom, top) の層の物質の成分（名前と、面積か割合）"""
    a = areas(problem)
    z = (bottom + top) / 2
    parts: dict[str, float] = {}
    if z < NOZZLE_BOTTOM_TOP or z > NOZZLE_TOP_BOTTOM:
        raise ValueError(z)
    # 案内管と計装管（ボトムノズルからトップノズルまで）
    parts["zircaloy"] = a["tubes"]
    if ROD_BOTTOM <= z < FUEL_BOTTOM or PLENUM_TOP <= z < ROD_TOP:
        parts["zircaloy"] += a["plugs"]
    elif FUEL_TOP <= z < PLENUM_TOP:
        parts["zircaloy"] += a["clad"]
        parts["void"] = a["plenum gas"]
    if TOP_GRID_CENTRE - END_GRID[3] / 2 <= z < TOP_GRID_CENTRE + END_GRID[3] / 2:
        parts["inconel"] = a["grid"]
    if problem == "3B":
        if POISON_TOP <= z < NOZZLE_TOP_BOTTOM:
            parts["ss304"] = a["pyrex clad"]
            parts["void"] = parts.get("void", 0.0) + a["pyrex gas"]
        if THIMBLE_BOTTOM <= z < THIMBLE_TOP:
            parts["ss304"] = parts.get("ss304", 0.0) + a["thimble"]
    return mixture(problem, parts)


def structure_name(bottom: Fraction) -> str:
    """燃料の外の構造の層の物質とセルの名前（下端の高さの `.` を `_` にしたもの。TOML のキーに `.` を使えないため）"""
    return "s" + decimal(bottom).replace(".", "_")


def physical_layers(problem: str) -> list[tuple[Fraction, Fraction, str]]:
    """平板の物理的な層 (下端, 上端, 物質の名前)（仕様書の高さ）。燃料は燃料のノードごと"""
    out = [
        (-OFFSET, Fraction(-5), "water"),
        (Fraction(-5), Fraction(0), "plate"),
        (Fraction(0), NOZZLE_BOTTOM_TOP, "bnozzle"),
    ]
    below = [NOZZLE_BOTTOM_TOP, ROD_BOTTOM, FUEL_BOTTOM]
    above = [FUEL_TOP, PLENUM_TOP, ROD_TOP, NOZZLE_TOP_BOTTOM]
    above += [TOP_GRID_CENTRE - END_GRID[3] / 2, TOP_GRID_CENTRE + END_GRID[3] / 2]
    if problem == "3B":
        above += [THIMBLE_BOTTOM, THIMBLE_TOP]
    for lo, hi in itertools.pairwise(below):
        out.append((lo, hi, structure_name(lo)))
    for lo, hi, name in fuel_nodes(problem):
        out.append((lo, hi, f"medium-{name}"))
    above = sorted(set(above))
    for lo, hi in itertools.pairwise(above):
        out.append((lo, hi, structure_name(lo)))
    out += [
        (NOZZLE_TOP_BOTTOM, NOZZLE_TOP_TOP, "tnozzle"),
        (NOZZLE_TOP_TOP, PLATE_TOP, "plate"),
        (PLATE_TOP, PLATE_TOP + 50, "water"),
    ]
    return out


def fuel_nodes(problem: str) -> list[tuple[Fraction, Fraction, str]]:
    """燃料の炉心のノード (下端, 上端, 平板の燃料の媒質の層の名前)。表 P3-5 の層の境界に挿入物の端を足したもの"""
    bounds = {FUEL_BOTTOM, *(Fraction(t) for t in LEVEL_TOPS)}
    if problem == "3B":
        bounds |= {PLUG_BOTTOM, POISON_BOTTOM, POISON_TOP}
    bounds = sorted(bounds)
    out = []
    for lo, hi in itertools.pairwise(bounds):
        out.append((lo, hi, insert_at(problem, (lo + hi) / 2)))
    return out


def insert_at(problem: str, z: Fraction) -> str:
    """燃料の高さ z の、支持格子のない格子の層の名前"""
    if problem == "3A":
        return "fuel"
    if z < PLUG_BOTTOM:
        return "empty"
    if z < POISON_BOTTOM:
        return "plug"
    if z < POISON_TOP:
        return "poison"
    return "plenum"


def grid_at(z: Fraction) -> str | None:
    """燃料の高さ z にある支持格子"""
    for name, *_, height, centres in (END_GRID, MID_GRID):
        for c in centres:
            if abs(z - Fraction(c)) < height / 2:
                return name
    return None


def fuel_layer(problem: str, z: Fraction) -> str:
    """燃料の高さ z の 2D の格子の層の名前"""
    base = insert_at(problem, z)
    grid = grid_at(z)
    suffix = {None: "", "zircaloy": "zgrid", "inconel": "igrid"}[grid]
    if problem == "3A":
        return suffix or "fuel"
    return f"{base}-{suffix}" if suffix else base


# 平板のノード（炉心の軸方向の反射体の層）: (名前, 仕様書の高さの下端, 上端)。ADR-0128 の表と同じ分け方
REFLECTOR_NODES = [
    ("lower-water", -OFFSET, Fraction(-5)),
    ("lower-plate", Fraction(-5), Fraction(0)),
    ("bottom-nozzle", Fraction(0), NOZZLE_BOTTOM_TOP),
    ("lower-gap-plug", NOZZLE_BOTTOM_TOP, FUEL_BOTTOM),
    ("plenum", FUEL_TOP, PLENUM_TOP),
    ("upper-plug-gap", PLENUM_TOP, NOZZLE_TOP_BOTTOM),
    ("top-nozzle", NOZZLE_TOP_BOTTOM, NOZZLE_TOP_TOP),
    ("upper-plate", NOZZLE_TOP_TOP, PLATE_TOP),
    ("upper-water", PLATE_TOP, PLATE_TOP + 50),
]


def nodes(problem: str) -> list[tuple[str, Fraction, Fraction, bool]]:
    """平板と炉心のノード (名前, 下端, 上端, 反射体か)（仕様書の高さ、下から）"""
    out = [(n, lo, hi, True) for n, lo, hi in REFLECTOR_NODES[:4]]
    for k, (lo, hi, _) in enumerate(fuel_nodes(problem)):
        out.append((f"fuel-{k + 1:02}", lo, hi, False))
    out += [(n, lo, hi, True) for n, lo, hi in REFLECTOR_NODES[4:]]
    return out


# ノズルと炉心板の鋼と水の割合（ADR-0128 の決定。examples/reflector/axial_slab.py と同じ）: 炉心板は 50% ずつ（表 5）、
# ノズルの水の割合は表 P3-3 の H-1 の比（ボトムノズル 0.72077、トップノズル 0.80853。3A と 3B で同じ）
FIXED_MIXTURES = {
    "plate": [("ss304", "0.5"), ("water", "0.5")],
    "bnozzle": [("ss304", "0.27923"), ("water", "0.72077")],
    "tnozzle": [("ss304", "0.19147"), ("water", "0.80853")],
}


def slab(problem: str) -> str:
    layers = physical_layers(problem)
    widths, cells = [], []
    edges = []  # 物理的な層ごとの (最初の列, 最後の列, 下端, 上端)
    start = 1
    for lo, hi, cell in layers:
        n, w = columns(hi - lo)
        widths.append(f"{{ size = {decimal(w)}, count = {n} }}")
        cells += [cell] * n
        edges.append((start, start + n - 1, lo, hi))
        start += n
    node_lines = []
    for name, lo, hi, reflector in nodes(problem):
        first = next(e[0] for e in edges if e[2] == lo)
        last = next(e[1] for e in edges if e[3] == hi)
        node_lines.append(
            f'  {{ output = "axial-{problem}/{name}.msgpack", columns = [{first}, {last}], rows = [1, 1]'
            + (", reflector = true }" if reflector else " }")
        )
    mixtures = []
    for name, parts in FIXED_MIXTURES.items():
        mixtures.append((name, parts))
    for lo, hi, cell in layers:
        if cell.startswith("s") and cell not in [m[0] for m in mixtures]:
            mixtures.append((cell, structure(problem, hi, lo)))
    mix_text = "\n\n".join(
        f"[material.{name}]\nmix = ["
        + ", ".join(f'{{ material = "{m}", volume = {v} }}' for m, v in parts)
        + "]"
        for name, parts in mixtures
    )
    media = sorted({c for _, _, c in layers if c.startswith("medium-")})
    media_text = "\n\n".join(
        f'[material.{m}]\nmedium = "medium-{problem}-{m.removeprefix("medium-")}.msgpack"'
        for m in media
    )
    used = []
    for c in cells:
        if c not in used:
            used.append(c)
    cell_text = "\n\n".join(
        f'[cell.{c}]\nzones = [{{ material = "{c}", rings = 1 }}]' for c in used
    )
    layer_lines = "\n".join(
        f"#   {decimal(lo)} 〜 {decimal(hi)} cm: {c}" for lo, hi, c in layers
    )
    names = ["water", "zircaloy", "ss304", "inconel"]
    return f"""# このファイルは p3_inputs.py が書く。手で直さない。
# VERA Problem {problem}（{PROBLEMS[problem]["description"]}）の軸方向の反射体の平板（ADR-0128）: 集合体の高さ全体
# （下の水の外端から上の水の外端まで）を格子の x（左が炉心の −z）に並べた1行の平板。左右は真空、上下は鏡面反射。
# 燃料は支持格子のない単一集合体の格子（lattice-{problem}-*.toml）が書いた媒質で、燃料の中の支持格子は入れない（ADR-0128）。
# 層の物質は、層の断面（ピッチ 21.50 cm の正方形）の中で体積で混ぜる（`s<下端の高さ>` は燃料の外の構造の層。成分は面積から:
# 案内管と計装管、端栓（被覆管の外半径の Zircaloy の棒）、プレナム（被覆管と真空）、上の端の支持格子（Inconel-718 の
# 質量 / 密度 / 高さ）、3B の Pyrex の棒のプレナム（SS304 の被覆管と真空）とシンブルプラグ（SS304 の棒））。自己遮蔽は
# しない（#0008）。物理的な層（仕様書の高さ、cm。下部炉心板の上面から）:
{layer_lines}
library = "../vera.msgpack"
self-shielding = false
pitch = {ROW}
azimuthal = 32
spacing = 0.05
boundary = {{ left = "vacuum", right = "vacuum" }}
few-group = {{ bounds = [1e7, 0.625, 1e-5], axial = true, nodes = [
{",\n".join(node_lines)},
] }}
widths = [
  {(",\n  ").join(widths)},
]
heights = [{{ size = {ROW} }}]
rows = ["{" ".join(cells)}"]

{media_text}

{materials(problem, names)}

{mix_text}

{cell_text}
"""


def core_files(problem: str) -> dict[str, str]:
    """炉心の4つの入力（ADR-0105）: f* の方法（`core-<問題>`）と、平板の定数だけの等価性（`core-<問題>-slab`）"""
    ns = nodes(problem)
    top = decimal(ns[-1][2] + OFFSET)
    mesh = ", ".join(decimal(lo + OFFSET) for _, lo, _, _ in ns) + f", {top}"
    files = {}
    for kind in ["fstar", "slab"]:
        layers = []
        for name, lo, hi, reflector in ns:
            if kind == "fstar" and not reflector:
                layers.append(
                    f"lattice-{problem}/{fuel_layer(problem, (lo + hi) / 2)}.msgpack"
                )
            else:
                layers.append(f"axial-{problem}/{name}.msgpack")
        stem = f"core-{problem}" + ("-slab" if kind == "slab" else "")
        what = (
            "燃料は 2D の格子の定数、軸方向の反射体は平板の定数（f* の方法。ADR-0128）"
            if kind == "fstar"
            else "すべての層を平板の定数にした、定数を作った体系の等価性の確認（ピッチは平板の行の高さ）"
        )
        pins = (
            f'pin-powers = "{stem}-pins.msgpack"\naxial-pin-powers = "{stem}-layers.msgpack"\n'
            if kind == "fstar"
            else ""
        )
        files[f"{stem}.toml"] = (
            f"# このファイルは p3_inputs.py が書く。VERA Problem {problem} の炉心: {what}。\n"
            "# 径方向は鏡面反射（集合体 1 体）、z の両端は真空で、平板の定数の真空の辺の比 J/φ（ADR-0101、ADR-0128）。\n"
            "# z は下の水の外端から（仕様書の高さ + 55 cm）\n"
            f'cycle = "{stem}-cycle.toml"\n'
            f"axial-mesh = [{mesh}]\n"
            "radial-nodes = 1\n"
            'boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], '
            'z = ["vacuum", "vacuum"] }\n' + pins
        )
        files[f"{stem}-cycle.toml"] = (
            "# このファイルは p3_inputs.py が書く\n"
            f'plant = "{stem}-plant.toml"\n'
            f'inventory = "{stem}-assemblies.toml"\n'
            f'loading = ["{problem}"]\n'
        )
        pitch = "21.5" if kind == "fstar" else ROW
        files[f"{stem}-plant.toml"] = (
            "# このファイルは p3_inputs.py が書く。集合体 1 体（ピッチ (cm) は集合体のピッチ 21.50 cm、等価性の確認では平板の\n"
            "# 行の高さ）\n"
            f'pitch = {pitch}\nheight = {top}\nmap = ["*"]\n'
        )
        files[f"{stem}-assemblies.toml"] = (
            "# このファイルは p3_inputs.py が書く\n"
            f"[design.{problem}]\nbounds = [{mesh}]\nlayers = [\n"
            + "".join(f'  "{layer}",\n' for layer in layers)
            + f']\n\n[assembly.{problem}]\ndesign = "{problem}"\n'
        )
    return files


def code(crate: str) -> list[str]:
    """Rust の段階の deps のうち、段階が呼ぶ実行ファイル（../../../dvc.yaml の build が作る。ADR-0159）"""
    return [f"../../../bin/{crate}"]


def stage(name: str, cmd: str, deps: list[str], outs: list[str]) -> str:
    return (
        f"  {name}:\n    cmd: {cmd}\n    deps:\n"
        + "".join(f"      - {d}\n" for d in deps)
        + "    outs:\n"
        + "".join(f"      - {o}\n" for o in outs)
    )


def pipeline() -> str:
    """段階: ライブラリ（../dvc.yaml の xs）→ 2D の格子（層ごと）→ 平板 → 炉心（等価性と f* の方法）→ 比較"""
    stages = []
    for problem in PROBLEMS:
        media = []
        fuel_files = []
        for name, _, grid, _ in LATTICES[problem]:
            stem = f"lattice-{problem}-{name}"
            outs = [f"lattice-{problem}/{name}.msgpack", f"{stem}.out"]
            if grid is None:
                media.append(f"medium-{problem}-{name}.msgpack")
                outs.append(media[-1])
            fuel_files.append(outs[0])
            stages.append(
                stage(
                    stem,
                    f"mkdir -p lattice-{problem} && ../../../bin/lattice {stem}.toml > {stem}.out",
                    [f"{stem}.toml", "../vera.msgpack", *code("lattice")],
                    outs,
                )
            )
        # 平板の少数群定数はディレクトリに書く（DVC は段階の前に出力を消すので、作ってから実行する）。MOC の k と
        # ノード平均の中性子束を比べるため、ログを debug で受ける
        axial = f"axial-{problem}"
        stages.append(
            stage(
                axial,
                f"mkdir -p {axial} && RUST_LOG=debug ../../../bin/lattice {axial}.toml > {axial}.out 2> {axial}.log",
                [f"{axial}.toml", "../vera.msgpack", *media, *code("lattice")],
                [axial, f"{axial}.out", f"{axial}.log"],
            )
        )
        for kind in ["slab", "fstar"]:
            stem = f"core-{problem}" + ("-slab" if kind == "slab" else "")
            inputs = [
                f"{stem}{s}.toml" for s in ["", "-cycle", "-plant", "-assemblies"]
            ]
            outs = [f"{stem}.out", f"{stem}.log"]
            if kind == "fstar":
                outs += [f"{stem}-pins.msgpack", f"{stem}-layers.msgpack"]
            used = [
                f
                for f in fuel_files
                if kind == "fstar"
                and f in core_files(problem)[f"{stem}-assemblies.toml"]
            ]
            stages.append(
                stage(
                    stem,
                    f"RUST_LOG=debug ../../../bin/core {stem}.toml > {stem}.out 2> {stem}.log",
                    [*inputs, axial, *used, *code("core")],
                    outs,
                )
            )
        stages.append(
            stage(
                f"compare-{problem}-slab",
                f"python3 ../../reflector/compare_axial.py {axial}.log core-{problem}-slab.log > compare-{problem}-slab.out",
                [
                    "../../reflector/compare_axial.py",
                    f"{axial}.log",
                    f"core-{problem}-slab.log",
                ],
                [f"compare-{problem}-slab.out"],
            )
        )
        stages.append(
            stage(
                f"compare-{problem}",
                f"python3 ../compare_p3.py {problem} ../../../references/vera/spec.txt core-{problem}.out "
                f"core-{problem}-pins.msgpack core-{problem}-layers.msgpack > compare-{problem}.out",
                [
                    "../compare_p3.py",
                    "../compare_keno.py",
                    "../../../references/vera/spec.txt",
                    f"core-{problem}.out",
                    f"core-{problem}-pins.msgpack",
                    f"core-{problem}-layers.msgpack",
                ],
                [f"compare-{problem}.out"],
            )
        )
        # OpenMC の参照解は DVC のステージの出力でない（AWS の VM の実行 m3 から取り出した references/ のファイル）ので、
        # 比べる 4 つのファイルを deps に書く。表を描くツール（examples/verification_tables.py）が JSON を読む（ADR-0176）
        openmc = f"../../../references/openmc/aws/m3/openmc-{problem}"
        stages.append(
            stage(
                f"compare-{problem}-openmc",
                f"python3 ../compare_p3.py --json compare-{problem}-openmc.json {problem} ../../../references/vera/spec.txt "
                f"core-{problem}.out core-{problem}-pins.msgpack core-{problem}-layers.msgpack {openmc} > compare-{problem}-openmc.out",
                [
                    "../compare_p3.py",
                    "../compare_keno.py",
                    "../../../references/vera/spec.txt",
                    *(
                        f"{openmc}/{f}"
                        for f in [
                            "k.json",
                            "layers.msgpack",
                            "pins.msgpack",
                            "pins-std.msgpack",
                        ]
                    ),
                    f"core-{problem}.out",
                    f"core-{problem}-pins.msgpack",
                    f"core-{problem}-layers.msgpack",
                ],
                [f"compare-{problem}-openmc.out", f"compare-{problem}-openmc.json"],
            )
        )
    return (
        "# このファイルは p3_inputs.py が書く。手で直さない。\n"
        "# VERA Problem 3A、3B を xlc で解く段階: 多群ライブラリ（../dvc.yaml の xs）→ 燃料の高さの 2D の格子（層ごとの\n"
        "# 少数群定数と、平板の燃料の媒質）→ 軸方向の反射体の平板（axial-<問題>）→ 炉心（平板の定数だけの等価性\n"
        "# core-<問題>-slab と、f* の方法 core-<問題>）→ 比較（平板の MOC との等価性、KENO-VI、OpenMC）。実行するのは `dvc repro`\n"
        "# （ADR-0159）。deps の書き方は examples/m0/dvc.yaml と同じ（段階の入力と、段階が呼ぶ実行ファイル）\n"
        "stages:\n" + "".join(stages)
    )


def main() -> None:
    (HERE / "dvc.yaml").write_text(pipeline())
    for problem in PROBLEMS:
        for name, insert, grid, what in LATTICES[problem]:
            (HERE / f"lattice-{problem}-{name}.toml").write_text(
                lattice(problem, name, insert, grid, what)
            )
        (HERE / f"axial-{problem}.toml").write_text(slab(problem))
        for name, text in core_files(problem).items():
            (HERE / name).write_text(text)


if __name__ == "__main__":
    main()
