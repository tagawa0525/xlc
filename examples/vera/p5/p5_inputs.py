"""VERA Problem 5-2D（炉心の高さの中央の 2D の炉心、零出力）を xlc で解く入力と段階を書く。

このファイルのある場所に次を書く。寸法と組成の出典は下の定数のコメントにあり、書いた入力のコメントにも写す。

- `assembly-<種類>.toml`: 集合体の種類（濃縮度と、Pyrex の棒の本数か D バンクの制御棒）ごとの単一集合体の格子（外周は
  鏡面反射）。少数群定数（形状関数つき）を `lattice/<種類>.msgpack` に書く
- `reflector-I.toml`、`reflector-corner.toml`: 径方向の反射体の I 型と角型の格子（ADR-0095、ADR-0132）。燃料は外周の
  集合体の代表の 3.10%、Pyrex の棒なし（ADR-0084）。ノードごとの少数群定数を `lattice/` に書く
- `small-moc.toml`: 3 × 3 の 3.10% の集合体と反射体の輪の 2D の全体を、真空の外周で MOC で解く格子（最小の入力の参照）
- `core-<ケース>.toml` ほか（ADR-0105 の 4 つ）: 全炉心（ADR-0090）の炉心の計算。ケースは次のとおり
  - `5A`、`5B`、`5C`: 5A-2D（制御棒なし）、5B-2D（D バンクに AIC）、5C-2D（D バンクに B4C）
  - `5A-Ly`: 5A の L 型の位置の I 型の定数の f* の辺を、もう一方の燃料の面に向けたもの（L 型の近似の大きさの見積もり）
  - `5A-Inext`: 5A の、角型の位置に隣る I 型の位置に、角型の格子の右と上の定数でなく I 型の定数を置いたもの
  - `small-E`、`small-I`: `small-moc.toml` と同じ体系を、角型の位置に隣る位置に角型の格子の右と上の定数を置いて、
    または I 型の定数を置いて解いたもの
  - `small-step-Lx`、`small-step-Ly`: 階段の外形（L 型、角型、`-`）を持つ小さな炉心（集合体 21 体）。L 型の位置の
    f* の辺の向きを変えた 2 つ
- `dvc.yaml`: 段階（ADR-0123）

炉心の図（`map`、`loading`）は上の行が北（1 行）、左の列が西（R 列）。全炉心の装荷は、図 M3-1 の 1/4（H〜A 列、8〜15 行）を
H 列と 8 行の中心を通る 2 つの対称面で鏡映して広げる（openmc_core_2d.py と同じ解釈）。xlc は対称面がピンの中心を
通る半分の集合体を表せない（#0067）ので、全炉心を解く（ADR-0090）。

使い方: python3 p5_inputs.py
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent

TEMPERATURE = 565
# 少数群の境界 (eV)（examples/reflector と同じ 2 群）
BOUNDS = "[1e7, 0.625, 1e-5]"
# 集合体の格子の方位角と間隔（ADR-0111）。反射体の格子と最小の入力の MOC も同じにする
AZIMUTHAL = 30
SPACING = 0.03
# 炉心の 2D の高さ (cm)。z の両面は鏡面反射なので、値は解に効かない
HEIGHT = "21.5"
PITCH = "21.5"

# 表 M3-4 の個数密度 (/(barn cm))。冷却材は 0.743 g/cc、ホウ素 1300 ppm（表 M3-1。assembly-2A.toml と同じ）
FUELS = {
    "21": "U-234 = 4.04814e-06, U-235 = 4.88801e-04, U-236 = 2.23756e-06, U-238 = 2.23844e-02, O-16 = 4.57591e-02",
    "26": "U-234 = 5.09503e-06, U-235 = 6.06709e-04, U-236 = 2.76809e-06, U-238 = 2.22663e-02, O-16 = 4.57617e-02",
    "31": "U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02",
}
ENRICHMENT = {"21": "2.11%", "26": "2.619%", "31": "3.10%"}
WATER = "H-1 = 4.96224e-02, O-16 = 2.48112e-02, B-10 = 1.07070e-05, B-11 = 4.30971e-05"

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

# 天然の C は C-12 98.93%、C-13 1.07% に分ける（ADR-0079）
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

PYREX = "B-10 = 9.63266e-04, B-11 = 3.90172e-03, O-16 = 4.67761e-02, Si-28 = 1.81980e-02, Si-29 = 9.24474e-04, Si-30 = 6.10133e-04"

AIC = """Ag-107 = 2.36159e-02
Ag-109 = 2.19403e-02
Cd-106 = 3.41523e-05
Cd-108 = 2.43165e-05
Cd-110 = 3.41250e-04
Cd-111 = 3.49720e-04
Cd-112 = 6.59276e-04
Cd-113 = 3.33873e-04
Cd-114 = 7.84957e-04
Cd-116 = 2.04641e-04
In-113 = 3.44262e-04
In-115 = 7.68050e-03"""

B4C = """B-10 = 1.52689e-02
B-11 = 6.14591e-02
C-12 = 1.89768e-02
C-13 = 2.05247e-04"""

# 17×17 の配置（行は上から）。`G` は案内管、`I` は計装管（examples/vera/p2/assembly-2A.toml と同じ）
LAYOUT = [
    ".................",
    ".................",
    ".....G..G..G.....",
    "...G.........G...",
    ".................",
    "..G..G..G..G..G..",
    ".................",
    ".................",
    "..G..G..I..G..G..",
    ".................",
    ".................",
    "..G..G..G..G..G..",
    ".................",
    "...G.........G...",
    ".....G..G..G.....",
    ".................",
    ".................",
]
# 図 5 の Pyrex の棒の配置（本数ごとの、1/8 の中の位置。中央からの行と列の数。openmc_model.py の PYREX_OCTANTS と同じ）
PYREX_OCTANTS = {
    8: {(3, 0), (5, 5)},
    12: {(3, 0), (6, 3)},
    16: {(3, 0), (5, 5), (6, 3)},
    20: {(3, 0), (5, 5), (6, 0), (6, 3)},
    24: {(3, 0), (3, 3), (5, 5), (6, 0), (6, 3)},
}

# 図 M3-1 の 1/4 炉心: 行（8〜15）ごとに H 列から（濃縮度, Pyrex の棒の本数）（openmc_core_2d.py の LOADING と同じ）
LOADING = (
    (("21", 0), ("26", 20), ("21", 0), ("26", 20), ("21", 0), ("26", 20), ("21", 0), ("31", 12)),
    (("26", 20), ("21", 0), ("26", 24), ("21", 0), ("26", 20), ("21", 0), ("31", 24), ("31", 0)),
    (("21", 0), ("26", 24), ("21", 0), ("26", 20), ("21", 0), ("26", 16), ("21", 0), ("31", 8)),
    (("26", 20), ("21", 0), ("26", 20), ("21", 0), ("26", 20), ("21", 0), ("31", 16), ("31", 0)),
    (("21", 0), ("26", 20), ("21", 0), ("26", 20), ("26", 0), ("26", 24), ("31", 0)),
    (("26", 20), ("21", 0), ("26", 16), ("21", 0), ("26", 24), ("31", 12), ("31", 0)),
    (("21", 0), ("31", 24), ("21", 0), ("31", 16), ("31", 0), ("31", 0)),
    (("31", 12), ("31", 0), ("31", 8), ("31", 0)),
)  # fmt: skip
# D バンクの位置（1/4 炉心の行, 列。H8 が (0, 0)。図 M3-1 の右）
BANK_D = {(0, 0), (0, 4), (4, 0), (4, 4)}
# 制御棒（表 8。吸収体の半径、被覆管は SS304 の 0.386〜0.484 cm、吸収体と被覆管の間は真空。ADR-0120）
RODS = {"aic": ("AIC", AIC, 0.382), "b4c": ("B4C", B4C, 0.373)}
CASES = {"5A": None, "5B": "aic", "5C": "b4c"}
# 全炉心の列の名前（西から。I、O、Q は使わない）
COLUMNS = "RPNMLKJHGFEDCBA"

# 方向（反時計回りの順）: 0 = +x（東）、1 = +y（北）、2 = −x（西）、3 = −y（南）。図の (行, 列) の動き
STEPS = [(0, 1), (-1, 0), (0, -1), (1, 0)]


def octant_images(octant: set[tuple[int, int]]) -> set[tuple[int, int]]:
    """1/8 の中の位置を、8 つの対称操作で移した集合体の位置（行, 列。0 始まり）にする。"""
    return {
        (8 + si * x, 8 + sj * y)
        for a, b in octant
        for x, y in ((a, b), (b, a))
        for si in (1, -1)
        for sj in (1, -1)
    }


def assembly_types() -> dict[str, tuple[str, int, str | None]]:
    """5-2D に現れる集合体の種類: 名前 → (濃縮度, Pyrex の棒の本数, 制御棒)"""
    types = {}
    for r, row in enumerate(LOADING):
        for c, (fuel, pyrex) in enumerate(row):
            types[assembly_name(fuel, pyrex, None)] = (fuel, pyrex, None)
            if (r, c) in BANK_D:
                if pyrex:
                    raise SystemExit("Pyrex rods in a rodded assembly")
                for rod in RODS:
                    types[assembly_name(fuel, 0, rod)] = (fuel, 0, rod)
    return dict(sorted(types.items()))


def assembly_name(fuel: str, pyrex: int, rod: str | None) -> str:
    return fuel + (f"-p{pyrex}" if pyrex else "") + (f"-{rod}" if rod else "")


def materials(names: list[str]) -> str:
    t = TEMPERATURE
    blocks = []
    for name in names:
        if name.startswith("fuel"):
            blocks.append(
                f"[material.{name}]\ntemperature = {t}\nnuclides = {{ {FUELS[name[4:]]} }}"
            )
        elif name == "water":
            blocks.append(
                f"[material.water]\ntemperature = {t}\nnuclides = {{ {WATER} }}"
            )
        elif name == "pyrex":
            blocks.append(
                f"[material.pyrex]\ntemperature = {t}\nnuclides = {{ {PYREX} }}"
            )
        else:
            body = {"zircaloy": ZIRCALOY, "ss304": SS304, "aic": AIC, "b4c": B4C}[name]
            blocks.append(
                f"[material.{name}]\ntemperature = {t}\n\n[material.{name}.nuclides]\n{body}"
            )
    return "\n\n".join(blocks)


def pin_cell(fuel: str) -> str:
    return f"""[cell.pin{fuel}]
zones = [
  {{ material = "fuel{fuel}", radius = 0.4096, rings = 5 }},
  {{ material = "void", radius = 0.418, rings = 1 }},
  {{ material = "zircaloy", radius = 0.475, rings = 1 }},
  {{ material = "water", radius = 0.63, rings = 3 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8"""


TUBE_CELLS = {
    "gt": """# 案内管（中は減速材）
[cell.gt]
zones = [
  { material = "water", radius = 0.561, rings = 3 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8""",
    "it": """# 計装管（中は減速材。シンブルは入れない。仕様書の問題 5-2D）
[cell.it]
zones = [
  { material = "water", radius = 0.559, rings = 3 },
  { material = "zircaloy", radius = 0.605, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8""",
    "py": """# Pyrex の棒を入れた案内管（表 6。内管の中と、Pyrex の内外の隙間は真空。ADR-0120）
[cell.py]
zones = [
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
]
sectors = 8""",
}


def rod_cell(rod: str) -> str:
    title, _, radius = RODS[rod]
    return f"""# {title} の制御棒を入れた案内管（表 8。吸収体と被覆管の間は真空。ADR-0120）
[cell.cr]
zones = [
  {{ material = "{rod}", radius = {radius}, rings = 5 }},
  {{ material = "void", radius = 0.386, rings = 1 }},
  {{ material = "ss304", radius = 0.484, rings = 1 }},
  {{ material = "water", radius = 0.561, rings = 1 }},
  {{ material = "zircaloy", radius = 0.602, rings = 1 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8"""


REFLECTOR_CELLS = """# 反射体の均質なマス（ピッチに合わないマスは長方形のセルになる）
[cell.w]
zones = [{ material = "water", rings = 1 }]

[cell.ss]
zones = [{ material = "ss304", rings = 1 }]"""


def pin_rows(fuel: str, pyrex: int = 0, rod: str | None = None) -> list[list[str]]:
    """集合体の 17 × 17 のセルの名前（行は上から）"""
    poisoned = octant_images(PYREX_OCTANTS[pyrex]) if pyrex else set()
    rows = []
    for i, line in enumerate(LAYOUT):
        row = []
        for j, mark in enumerate(line):
            if mark == "I":
                row.append("it")
            elif mark == "G":
                row.append("cr" if rod else ("py" if (i, j) in poisoned else "gt"))
            else:
                if (i, j) in poisoned:
                    raise SystemExit(f"Pyrex rod at {(i, j)} outside the guide tubes")
                row.append(f"pin{fuel}")
        rows.append(row)
    if sum(row.count("py") for row in rows) != pyrex:
        raise SystemExit(f"bad Pyrex pattern for {pyrex} rods")
    return rows


def rows_toml(rows: list[list[str]]) -> str:
    width = max(len(name) for row in rows for name in row)
    return (
        "rows = [\n"
        + "".join(
            '  "' + " ".join(name.ljust(width) for name in row).rstrip() + '",\n'
            for row in rows
        )
        + "]"
    )


def assembly(name: str, fuel: str, pyrex: int, rod: str | None) -> str:
    what = f"濃縮度 {ENRICHMENT[fuel]}"
    if pyrex:
        what += f"、Pyrex の棒 {pyrex} 本（図 5 の配置）"
    if rod:
        what += f"、24 本の案内管に {RODS[rod][0]} の制御棒（D バンク）"
    names = [f"fuel{fuel}", "zircaloy", "water"]
    cells = [pin_cell(fuel), TUBE_CELLS["gt"], TUBE_CELLS["it"]]
    if pyrex or rod:
        names.append("ss304")
    if pyrex:
        names.append("pyrex")
        cells.append(TUBE_CELLS["py"])
    if rod:
        names.append(rod)
        cells.append(rod_cell(rod))
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
# VERA Problem 5-2D の集合体（{what}）の単一集合体の格子（外周は鏡面反射）。燃料、被覆管、減速材は 565 K、
# 減速材 0.743 g/cc、ホウ素 1300 ppm（表 M3-1、M3-4）。ほかは examples/vera/p2 の 2A、2E、2G、2H と同じ形。少数群定数は
# 形状関数と角の不連続因子つき（ADR-0096）で、炉心の燃料のノードの定数（f* の方法）にする。
library = "../vera.msgpack"
few-group = {{ output = "lattice/{name}.msgpack", bounds = {BOUNDS}, pins = true }}
pitch = 1.26
azimuthal = {AZIMUTHAL}
spacing = {SPACING}
gap = {{ width = 0.04, material = "water" }}
{rows_toml(pin_rows(fuel, pyrex, rod))}

{materials(names)}

{chr(10).join(cells)}
"""


def columns_of(side: str) -> list[str]:
    """燃料の集合体 1 体の幅（余白 0.04、ピン 17、余白 0.04）と、反射体 1 つの幅（燃料側から隙間の水 0.19、
    SS304 のバッフル 2.85 cm（0.95 × 3）、水 18.46 cm（0.71 × 26））。仕様書の 1.13 節の表 16。
    `side` は "fuel"、"after"（燃料の右か下の反射体）、"before"（左か上）"""
    if side == "fuel":
        return ["{ size = 0.04 }", "{ size = 1.26, count = 17 }", "{ size = 0.04 }"]
    after = [
        "{ size = 0.19 }",
        "{ size = 0.95, count = 3 }",
        "{ size = 0.71, count = 26 }",
    ]
    return after if side == "after" else after[::-1]


def cell_marks(side: str) -> list[str]:
    """`columns_of` のセルの並びの種類: 余白 m、ピン p、隙間 g、バッフル b、外の水 w"""
    if side == "fuel":
        return ["m"] + ["p"] * 17 + ["m"]
    after = ["g"] + ["b"] * 3 + ["w"] * 26
    return after if side == "after" else after[::-1]


def core_lattice_rows(
    fuel_rows: list[list[str]], kinds: list[list[str]]
) -> list[list[str]]:
    """集合体と反射体を並べた格子のセルの名前。`kinds` は位置ごとの "F"（燃料）、"R"（反射体）で、行は上から。
    バッフルは燃料の外周から隙間 0.19 cm をおいた帯で、燃料の階段の形に沿い、角では帯が直角に会う
    （openmc_core_2d.py の、集合体の外周を正方形で太らせた形と同じ）。"""
    n_rows, n_cols = len(kinds), len(kinds[0])
    col_sides = side_list(
        [any(kinds[i][j] == "F" for i in range(n_rows)) for j in range(n_cols)]
    )
    row_sides = side_list(
        [any(kinds[i][j] == "F" for j in range(n_cols)) for i in range(n_rows)]
    )
    # セルの中心の座標（集合体のピッチを 1 とする）と、燃料の正方形の和からの距離で、隙間とバッフルを決める
    widths = {"m": 0.04, "p": 1.26, "g": 0.19, "b": 0.95, "w": 0.71}

    def centres(side_list: list[str]) -> list[tuple[int, float, str]]:
        out = []
        for k, side in enumerate(side_list):
            x = 0.0
            for mark in cell_marks(side):
                w = widths[mark]
                out.append((k, x + w / 2, mark))
                x += w
        return out

    xs, ys = centres(col_sides), centres(row_sides)
    fuel = {(i, j) for i in range(n_rows) for j in range(n_cols) if kinds[i][j] == "F"}
    pitch = 21.5

    def distance(i: int, y: float, j: int, x: float) -> float:
        """セルの中心から燃料の集合体の正方形の和までの距離（チェビシェフ距離。正方形で太らせた形の境界の判定）"""
        best = float("inf")
        for fi, fj in fuel:
            dx = max(0.0, abs((j * pitch + x) - (fj * pitch + pitch / 2)) - pitch / 2)
            dy = max(0.0, abs((i * pitch + y) - (fi * pitch + pitch / 2)) - pitch / 2)
            best = min(best, max(dx, dy))
        return best

    rows = []
    for i, y, _ in ys:
        row = []
        for j, x, _ in xs:
            if kinds[i][j] == "F":
                # 燃料の集合体の中: 余白は水、ほかはピンの配置
                yi = [m for k, _, m in ys if k == i]
                xj = [m for k, _, m in xs if k == j]
                a = sum(1 for k, yy, _ in ys if k == i and yy < y)
                b = sum(1 for k, xx, _ in xs if k == j and xx < x)
                if yi[a] == "m" or xj[b] == "m":
                    row.append("w")
                else:
                    row.append(fuel_rows[a - 1][b - 1])
            else:
                d = distance(i, y, j, x)
                row.append("ss" if 0.19 < d < 0.19 + 2.85 else "w")
        rows.append(row)
    return rows


def grid_widths(sides: list[str]) -> str:
    return "[\n" + "".join(f"  {c},\n" for s in sides for c in columns_of(s)) + "]"


def side_list(fuel_at: list[bool]) -> list[str]:
    out = []
    for k in range(len(fuel_at)):
        if fuel_at[k]:
            out.append("fuel")
        elif k + 1 < len(fuel_at) and fuel_at[k + 1]:
            out.append("before")
        else:
            out.append("after")
    return out


def block_lattice(
    kinds: list[list[str]], boundary: str, few_group: str, comment: str, extra: str = ""
) -> str:
    """3.10% の集合体と反射体を `kinds` のとおり並べた格子"""
    n_rows, n_cols = len(kinds), len(kinds[0])
    col_fuel = [any(kinds[i][j] == "F" for i in range(n_rows)) for j in range(n_cols)]
    row_fuel = [any(kinds[i][j] == "F" for j in range(n_cols)) for i in range(n_rows)]
    rows = core_lattice_rows(pin_rows("31"), kinds)
    names = ["fuel31", "zircaloy", "water", "ss304"]
    cells = [pin_cell("31"), TUBE_CELLS["gt"], TUBE_CELLS["it"], REFLECTOR_CELLS]
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
{comment}
library = "../vera.msgpack"
pitch = 1.26
azimuthal = {AZIMUTHAL}
spacing = {SPACING}
boundary = {boundary}
{few_group}{extra}widths = {grid_widths(side_list(col_fuel))}
heights = {grid_widths(side_list(row_fuel))}
{rows_toml(rows)}

{materials(names)}

{chr(10).join(cells)}
"""


def node_ranges(sides: list[str]) -> list[tuple[int, int]]:
    """位置ごとの格子のセルの列（1 始まりの最初と最後）"""
    out, start = [], 1
    for s in sides:
        n = len(cell_marks(s))
        out.append((start, start + n - 1))
        start += n
    return out


def reflector_lattices() -> dict[str, str]:
    common = (
        "# 燃料は外周の集合体の代表として 3.10%、Pyrex の棒なし（examples/vera/p2/assembly-2A.toml と同じ物質とセル。\n"
        "# ADR-0084）。反射体は燃料側から、集合体の外周からの隙間の水 0.19 cm、SS304 のバッフル 2.85 cm、水 18.46 cm\n"
        "# （仕様書の 1.13 節の表 16。集合体のピッチ 21.50 cm の幅）。ノードごとの形状関数と角の不連続因子を付ける（ADR-0131）。"
    )
    out = {}
    kinds = [["F", "R"]]
    cols = node_ranges(["fuel", "after"])
    rws = node_ranges(["fuel"])
    fg = (
        f"few-group = {{ bounds = {BOUNDS}, nodes = [\n"
        f'  {{ output = "lattice/reflector-I-fuel.msgpack", columns = [{cols[0][0]}, {cols[0][1]}], rows = [{rws[0][0]}, {rws[0][1]}] }},\n'
        f'  {{ output = "lattice/reflector-I.msgpack", columns = [{cols[1][0]}, {cols[1][1]}], rows = [{rws[0][0]}, {rws[0][1]}], reflector = true }},\n'
        "], pins = true }\n"
    )
    out["reflector-I.toml"] = block_lattice(
        kinds,
        '{ right = "vacuum" }',
        fg,
        "# 径方向の反射体の I 型の格子（ADR-0095）: VERA Problem 5-2D の燃料の集合体の右に反射体を置く。左、上、下は鏡面\n"
        "# 反射、右は真空（examples/reflector/reflector-I.toml と同じ形）。\n" + common,
    )
    kinds = [["R", "R"], ["F", "R"]]
    cols = node_ranges(["fuel", "after"])
    rws = node_ranges(["before", "fuel"])
    nodes = [
        ("corner-fuel", cols[0], rws[1], False),
        ("corner-right", cols[1], rws[1], True),
        ("corner-top", cols[0], rws[0], True),
        ("corner", cols[1], rws[0], True),
    ]
    fg = (
        f"few-group = {{ bounds = {BOUNDS}, nodes = [\n"
        + "".join(
            f'  {{ output = "lattice/reflector-{n}.msgpack", columns = [{c[0]}, {c[1]}], rows = [{r[0]}, {r[1]}]'
            + (", reflector = true" if refl else "")
            + " },\n"
            for n, c, r, refl in nodes
        )
        + "], pins = true }\n"
    )
    out["reflector-corner.toml"] = block_lattice(
        kinds,
        '{ right = "vacuum", top = "vacuum" }',
        fg,
        "# 径方向の反射体の角型の格子（ADR-0095、ADR-0132）: VERA Problem 5-2D の燃料の集合体の右と上に反射体を置き、バッフルは\n"
        "# 集合体の外側の角で L 字に折る。左と下は鏡面反射、右と上は真空（examples/reflector/reflector-corner.toml と同じ形）。\n"
        "# ノードは集合体、右と上の反射体（燃料に面で接する）、右上の反射体（燃料に角の点だけで接する）。\n"
        + common,
    )
    return out


def small_moc() -> str:
    kinds = [list(line) for line in SMALL_KINDS]
    return block_lattice(
        kinds,
        '{ left = "vacuum", right = "vacuum", top = "vacuum", bottom = "vacuum" }',
        "",
        "# 最小の入力の参照: 3.10% の集合体（Pyrex の棒なし）3 × 3 と反射体の輪（I 型の位置 12、角型の位置 4）の 2D の全体を、\n"
        "# 外周を真空にして MOC で解く。反射体の寸法とバッフルの形は reflector-I.toml、reflector-corner.toml と同じ。\n"
        "# core の small-E、small-I と比べる（セルごとの核分裂率）。",
        'pin-powers = "small-moc-pins.msgpack"\n',
    )


# 最小の入力の 3 × 3 の集合体と反射体の輪（F は燃料、R は反射体）
SMALL_KINDS = ["RRRRR", "RFFFR", "RFFFR", "RFFFR", "RRRRR"]


def fuel_set(counts: list[int]) -> set[tuple[int, int]]:
    """1/4 の行ごとの燃料の数（中心の行から。中心の集合体を含む）を、中心で鏡映した全体の位置（行, 列。中心が (0, 0)）"""
    return {
        (sr * r, sc * c)
        for r, n in enumerate(counts)
        for c in range(n)
        for sr in (1, -1)
        for sc in (1, -1)
    }


def turn_of(fr: int, to: int) -> int:
    return (to - fr) % 4


def rotate_diagonal(d: tuple[int, int], k: int) -> tuple[int, int]:
    """(x, y) の向きを反時計回りに 90° × k 回す"""
    x, y = d
    for _ in range(k):
        x, y = -y, x
    return (x, y)


def classify(fuel: set[tuple[int, int]], size: int) -> dict[tuple[int, int], tuple]:
    """図の位置（行, 列。中心が (0, 0)、行は下が正）ごとの種類: ("F",)、("I", 燃料の方向)、("L", 方向 1, 方向 2)、
    ("C", 燃料の対角の向き (x, y))、("-",)"""
    out = {}
    h = size // 2
    for i in range(-h, h + 1):
        for j in range(-h, h + 1):
            if (i, j) in fuel:
                out[(i, j)] = ("F",)
                continue
            faces = [d for d, (di, dj) in enumerate(STEPS) if (i + di, j + dj) in fuel]
            if len(faces) == 1:
                out[(i, j)] = ("I", faces[0])
            elif len(faces) == 2:
                out[(i, j)] = ("L", *faces)
            elif len(faces) > 2:
                raise SystemExit(
                    f"reflector at {(i, j)} faces fuel on {len(faces)} sides"
                )
            else:
                diag = [
                    (dj, -di)
                    for di in (1, -1)
                    for dj in (1, -1)
                    if (i + di, j + dj) in fuel
                ]
                if len(diag) > 1:
                    raise SystemExit(
                        f"reflector at {(i, j)} touches fuel at {len(diag)} corners"
                    )
                out[(i, j)] = ("C", diag[0]) if diag else ("-",)
    return out


def reflector_map(
    fuel: set[tuple[int, int]], size: int, l_rule: str, next_to_corner: str
) -> list[list[str]]:
    """炉心の形の `map` の記号（行は上から）。燃料は `*`、反射体は種類と回転（ADR-0105）。

    - I 型の位置は I 型の格子の反射体（燃料が −x）を、燃料の方向に回す
    - 角型の位置は角型の格子の右上の反射体（燃料が −x −y の対角）を回す
    - 角型の位置に隣る I 型の位置は、`next_to_corner` が "E" なら角型の格子の右（E。燃料が −x、角型が +y）か上（N。
      燃料が −y、角型が +x）の反射体を、燃料と角型の位置の向きに合うほうで置く（ADR-0132 の帰結）。"I" なら I 型
    - L 型の位置は I 型の定数（ADR-0132。#0068 まで）。f* の辺は、`l_rule` が "x" なら北東と南西の象限で x の向きの燃料、
      北西と南東の象限で y の向きの燃料に向ける（炉心の中心のまわりの 4 回の回転で同じ向き）。"y" はその逆
    """
    kinds = classify(fuel, size)
    h = size // 2
    rows = []
    for i in range(-h, h + 1):
        row = []
        for j in range(-h, h + 1):
            kind = kinds[(i, j)]
            if kind[0] == "F":
                row.append("*")
            elif kind[0] == "-":
                row.append("-")
            elif kind[0] == "C":
                k = next(k for k in range(4) if rotate_diagonal((-1, -1), k) == kind[1])
                row.append(turned("C", k))
            elif kind[0] == "L":
                d1, d2 = kind[1], kind[2]
                x_dir = d1 if d1 in (0, 2) else d2
                y_dir = d2 if x_dir == d1 else d1
                # 象限: 燃料は中心の側にあるので、燃料の向きの反対が象限の向き
                east, north = x_dir == 2, y_dir == 3
                first = east == north  # 北東か南西
                face = x_dir if (first == (l_rule == "x")) else y_dir
                row.append(turned("I", turn_of(2, face)))
            else:
                d = kind[1]
                corners = [
                    e
                    for e, (di, dj) in enumerate(STEPS)
                    if kinds.get((i + di, j + dj), ("-",))[0] == "C"
                ]
                if next_to_corner == "E" and corners:
                    if len(corners) > 1:
                        raise SystemExit(f"I site at {(i, j)} next to two corner sites")
                    c = corners[0]
                    k = turn_of(2, d)
                    if (1 + k) % 4 == c:
                        row.append(turned("E", k))
                    else:
                        k = turn_of(3, d)
                        if (0 + k) % 4 != c:
                            raise SystemExit(f"no corner-lattice node fits {(i, j)}")
                        row.append(turned("N", k))
                else:
                    row.append(turned("I", turn_of(2, d)))
        rows.append(row)
    return rows


def turned(name: str, k: int) -> str:
    return name if k == 0 else f"{name}@{90 * k}"


REFLECTORS = {
    "I": ("lattice/reflector-I.msgpack", "I 型の格子の反射体（燃料が −x）"),
    "E": (
        "lattice/reflector-corner-right.msgpack",
        "角型の格子の右の反射体（燃料が −x、角型が +y）",
    ),
    "N": (
        "lattice/reflector-corner-top.msgpack",
        "角型の格子の上の反射体（燃料が −y、角型が +x）",
    ),
    "C": (
        "lattice/reflector-corner.msgpack",
        "角型の格子の右上の反射体（燃料が −x −y の対角）",
    ),
}


def grid_text(rows: list[list[str]]) -> str:
    width = max(len(s) for row in rows for s in row)
    return (
        "[\n"
        + "".join(
            '  "' + " ".join(s.ljust(width) for s in row).rstrip() + '",\n'
            for row in rows
        )
        + "]"
    )


def plant(map_rows: list[list[str]], comment: str) -> str:
    used = sorted({s.split("@")[0] for row in map_rows for s in row} & set(REFLECTORS))
    blocks = "".join(
        f'\n# {REFLECTORS[n][1]}\n[reflector.{n}]\nbounds = [0.0, {HEIGHT}]\nlayers = ["{REFLECTORS[n][0]}"]\n'
        for n in used
    )
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
{comment}
pitch = {PITCH}
height = {HEIGHT}
map = {grid_text(map_rows)}
{blocks}"""


def case(name: str, comment: str, outline: bool = True) -> str:
    """計算のケース。`outline` は炉心の形に `-` があること（外形の面の条件を書く。ADR-0132）"""
    face = ', outline = "vacuum"' if outline else ""
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
{comment}
cycle = "core-{name}-cycle.toml"
axial-mesh = [0.0, {HEIGHT}]
radial-nodes = 1
pin-powers = "core-{name}-pins.msgpack"
boundary = {{ x = ["vacuum", "vacuum"], y = ["vacuum", "vacuum"], z = ["reflective", "reflective"]{face} }}
"""


def cycle(name: str, plant_file: str, inventory: str, loading: list[list[str]]) -> str:
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
plant = "{plant_file}"
inventory = "{inventory}"
loading = {grid_text(loading)}
"""


def assembly_id(i: int, j: int) -> str:
    """全炉心の位置（中心が (0, 0)、行は下が正）の集合体の ID（列の名前と行の番号。H08 が中心）"""
    return f"{COLUMNS[j + 7]}{i + 8:02d}"


def core_5_2d(rod: str | None) -> tuple[dict[tuple[int, int], str], str]:
    """5-2D の全炉心の、位置ごとの集合体の種類と台帳"""
    designs = {}
    for i in range(-7, 8):
        for j in range(-7, 8):
            r, c = abs(i), abs(j)
            if c >= len(LOADING[r]):
                continue
            fuel, pyrex = LOADING[r][c]
            designs[(i, j)] = assembly_name(
                fuel, pyrex, rod if (r, c) in BANK_D else None
            )
    return designs, inventory(designs)


def inventory(designs: dict[tuple[int, int], str]) -> str:
    used = sorted(set(designs.values()))
    text = "# このファイルは p5_inputs.py が書く。手で直さない。\n"
    text += "# 集合体の台帳。設計は集合体の種類の単一集合体の格子の定数（f* の方法）。ID は全炉心の位置（列の名前と行の番号）\n"
    text += "".join(
        f'\n[design.{d}]\nbounds = [0.0, {HEIGHT}]\nlayers = ["lattice/{d}.msgpack"]\n'
        for d in used
    )
    text += "".join(
        f'\n[assembly.{assembly_id(i, j)}]\ndesign = "{d}"\n'
        for (i, j), d in sorted(designs.items())
    )
    return text


def loading_of(designs: dict[tuple[int, int], str], size: int) -> list[list[str]]:
    h = size // 2
    return [
        [assembly_id(i, j) if (i, j) in designs else "." for j in range(-h, h + 1)]
        for i in range(-h, h + 1)
    ]


# 5-2D の 1/4 の行ごとの燃料の数（8〜15 行。H8 を含む。図 M3-1）
COUNTS_5 = [len(row) for row in LOADING]
# 階段の外形の小さな炉心の 1/4 の行ごとの燃料の数（中心を含む）。集合体 21 体
COUNTS_STEP = [3, 3, 2]


def core_files() -> dict[str, str]:
    files = {}
    size = 17
    fuel = fuel_set(COUNTS_5)
    plants = {
        "": (
            "x",
            "E",
            "L 型の位置の f* の辺は北東と南西の象限で x の向き、角型の位置に隣る I 型の位置は角型の格子の右と上の反射体",
        ),
        "-Ly": (
            "y",
            "E",
            "L 型の位置の f* の辺は北東と南西の象限で y の向き（5A-Ly）、角型の位置に隣る位置は core-plant.toml と同じ",
        ),
        "-Inext": (
            "x",
            "I",
            "L 型の位置は core-plant.toml と同じ、角型の位置に隣る I 型の位置も I 型の格子の反射体（5A-Inext）",
        ),
    }
    for suffix, (rule, nxt, what) in plants.items():
        files[f"core-plant{suffix}.toml"] = plant(
            reflector_map(fuel, size, rule, nxt),
            "# VERA Problem 5-2D の全炉心の形（ADR-0105、ADR-0132）: 集合体 193 体（*）、燃料に面か角の点で接する反射体の輪、\n"
            "# その外は `-`。反射体の位置は燃料との接し方で、I 型、L 型（I 型の定数。#0068 まで）、角型に分ける。\n# "
            + what
            + "。",
        )
    for name, rod in CASES.items():
        designs, inv = core_5_2d(rod)
        files[f"core-{name}-assemblies.toml"] = inv
        what = {
            None: "制御棒なし",
            "aic": "D バンク（H8、D8、H12、D12 と鏡映の位置の 16 体）に AIC の制御棒",
            "b4c": "D バンク（H8、D8、H12、D12 と鏡映の位置の 16 体）に B4C の制御棒",
        }[rod]
        variants = [("", "")] + (
            [("-Ly", "-Ly"), ("-Inext", "-Inext")] if name == "5A" else []
        )
        for v, p in variants:
            files[f"core-{name}{v}-cycle.toml"] = cycle(
                name,
                f"core-plant{p}.toml",
                f"core-{name}-assemblies.toml",
                loading_of(designs, size),
            )
            files[f"core-{name}{v}.toml"] = case(
                f"{name}{v}",
                f"# VERA Problem {name}-2D（{what}）の全炉心の 2D の計算。外周の直方体と外形の面は真空（ADR-0132）、z は鏡面反射。\n"
                "# 制御棒の位置はまだ入力にない（ADR-0105 の `rods`）ので、制御棒の入った集合体は台帳で別の設計にする。",
            )
    # 最小の入力: 3 × 3 の集合体と反射体の輪
    small = {
        (i - 2, j - 2)
        for i, line in enumerate(SMALL_KINDS)
        for j, k in enumerate(line)
        if k == "F"
    }
    small_designs = {p: "31" for p in small}
    files["core-small-assemblies.toml"] = inventory(small_designs)
    for nxt in ["E", "I"]:
        files[f"core-small-{nxt}-plant.toml"] = plant(
            reflector_map(small, 5, "x", nxt),
            "# 最小の入力（small-moc.toml の体系）の炉心の形。角型の位置に隣る I 型の位置は、"
            + ("角型の格子の右と上の反射体" if nxt == "E" else "I 型の格子の反射体")
            + "。",
        )
        files[f"core-small-{nxt}-cycle.toml"] = cycle(
            f"small-{nxt}",
            f"core-small-{nxt}-plant.toml",
            "core-small-assemblies.toml",
            loading_of(small_designs, 5),
        )
        files[f"core-small-{nxt}.toml"] = case(
            f"small-{nxt}",
            "# 最小の入力（small-moc.toml の体系）を core で解く。外周の直方体は真空、z は鏡面反射。",
            outline=False,
        )
    step = fuel_set(COUNTS_STEP)
    step_designs = {p: "31" for p in step}
    files["core-small-step-assemblies.toml"] = inventory(step_designs)
    for rule in ["x", "y"]:
        name = f"small-step-L{rule}"
        files[f"core-{name}-plant.toml"] = plant(
            reflector_map(step, 7, rule, "E"),
            "# 階段の外形の小さな炉心（3.10% の集合体 21 体。1/4 の行ごとの燃料の数 3、3、2）の形。L 型の位置の f* の辺は北東と\n"
            f"# 南西の象限で {rule} の向き。角型の位置に隣る I 型の位置は角型の格子の右と上の反射体。",
        )
        files[f"core-{name}-cycle.toml"] = cycle(
            name,
            f"core-{name}-plant.toml",
            "core-small-step-assemblies.toml",
            loading_of(step_designs, 7),
        )
        files[f"core-{name}.toml"] = case(
            name, "# 階段の外形の小さな炉心を core で解く（外形と対称性の確認）。"
        )
    return files


def code(crate: str) -> list[str]:
    """Rust の段階の deps のうち、バイナリを作るコード（examples/m0/dvc.yaml と同じ並び。ADR-0123）"""
    return [
        f"../../../crates/{crate}/src",
        f"../../../crates/{crate}/Cargo.toml",
        "../../../crates/numeric/src",
        "../../../crates/numeric/Cargo.toml",
        "../../../crates/format/src",
        "../../../crates/format/Cargo.toml",
        "../../../Cargo.toml",
        "../../../Cargo.lock",
        "../../../rust-toolchain.toml",
        "../../../.cargo/config.toml",
        "../../../flake.lock",
    ]


def stage(name: str, cmd: str, deps: list[str], outs: list[str]) -> str:
    return (
        f"  {name}:\n    cmd: {cmd}\n    deps:\n"
        + "".join(f"      - {d}\n" for d in deps)
        + "    outs:\n"
        + "".join(f"      - {o}\n" for o in outs)
    )


def layer_files(text: str) -> list[str]:
    return sorted(
        {line.split('"')[1] for line in text.splitlines() if line.startswith("layers")}
    )


def pipeline(files: dict[str, str]) -> str:
    stages = []
    for name in assembly_types():
        stem = f"assembly-{name}"
        stages.append(
            stage(
                stem,
                f"mkdir -p lattice && cargo run --release -p xlc-lattice -- {stem}.toml > {stem}.out",
                [f"{stem}.toml", "../vera.msgpack", *code("lattice")],
                [f"lattice/{name}.msgpack", f"{stem}.out"],
            )
        )
    for stem, outs in [
        ("reflector-I", ["reflector-I-fuel", "reflector-I"]),
        (
            "reflector-corner",
            [
                "reflector-corner-fuel",
                "reflector-corner-right",
                "reflector-corner-top",
                "reflector-corner",
            ],
        ),
    ]:
        stages.append(
            stage(
                stem,
                f"mkdir -p lattice && cargo run --release -p xlc-lattice -- {stem}.toml > {stem}.out",
                [f"{stem}.toml", "../vera.msgpack", *code("lattice")],
                [*(f"lattice/{o}.msgpack" for o in outs), f"{stem}.out"],
            )
        )
    stages.append(
        stage(
            "small-moc",
            "cargo run --release -p xlc-lattice -- small-moc.toml > small-moc.out",
            ["small-moc.toml", "../vera.msgpack", *code("lattice")],
            ["small-moc-pins.msgpack", "small-moc.out"],
        )
    )
    cases = [
        f[5:-5]
        for f in files
        if f.startswith("core-")
        and f.count("-") >= 1
        and not any(
            f.endswith(s) for s in ("-cycle.toml", "-plant.toml", "-assemblies.toml")
        )
        and not f.startswith("core-plant")
    ]
    for name in sorted(cases):
        cyc = f"core-{name}-cycle.toml"
        cycle_text = files[cyc]
        plant_file = cycle_text.split('plant = "')[1].split('"')[0]
        inv_file = cycle_text.split('inventory = "')[1].split('"')[0]
        constants = sorted(
            set(layer_files(files[plant_file]) + layer_files(files[inv_file]))
        )
        stages.append(
            stage(
                f"core-{name}",
                f"RUST_LOG=debug cargo run --release -p xlc-core -- core-{name}.toml > core-{name}.out 2> core-{name}.log",
                [
                    f"core-{name}.toml",
                    cyc,
                    plant_file,
                    inv_file,
                    *constants,
                    *code("core"),
                ],
                [f"core-{name}.out", f"core-{name}.log", f"core-{name}-pins.msgpack"],
            )
        )
    for name in ["small-E", "small-I"]:
        stages.append(
            stage(
                f"compare-{name}",
                f"python3 ../compare_p5.py moc small-moc.out small-moc-pins.msgpack core-{name}.out core-{name}-pins.msgpack > compare-{name}.out",
                [
                    "../compare_p5.py",
                    "small-moc.out",
                    "small-moc-pins.msgpack",
                    f"core-{name}.out",
                    f"core-{name}-pins.msgpack",
                ],
                [f"compare-{name}.out"],
            )
        )
    for name in ["small-step-Lx", "small-step-Ly"]:
        stages.append(
            stage(
                f"compare-{name}",
                f"python3 ../compare_p5.py symmetry core-{name}.out core-{name}-pins.msgpack > compare-{name}.out",
                ["../compare_p5.py", f"core-{name}.out", f"core-{name}-pins.msgpack"],
                [f"compare-{name}.out"],
            )
        )
    for name in ["5A", "5A-Ly", "5A-Inext", "5B", "5C"]:
        problem = name.split("-")[0]
        stages.append(
            stage(
                f"compare-{name}",
                f"python3 ../compare_p5.py keno {problem} ../../../references/vera/spec.txt core-{name}.out core-{name}-pins.msgpack > compare-{name}.out",
                [
                    "../compare_p5.py",
                    "../../../references/vera/spec.txt",
                    f"core-{name}.out",
                    f"core-{name}-pins.msgpack",
                ],
                [f"compare-{name}.out"],
            )
        )
    return (
        "# このファイルは p5_inputs.py が書く。手で直さない。\n"
        "# VERA Problem 5-2D を xlc で解く段階: 多群ライブラリ（../dvc.yaml の xs）→ 単一集合体の格子（集合体の種類ごと）、\n"
        "# 反射体の I 型と角型の格子、最小の入力の MOC → 炉心（5A、5B、5C と、近似の見積もりと最小の入力の版）→ 比較（KENO-VI、\n"
        "# 最小の入力の MOC、対称性）。走らせるのは `dvc repro`（ADR-0123）。deps の書き方は examples/m0/dvc.yaml と同じ\n"
        "stages:\n" + "".join(stages)
    )


def main() -> None:
    files = core_files()
    for name, (fuel, pyrex, rod) in assembly_types().items():
        files[f"assembly-{name}.toml"] = assembly(name, fuel, pyrex, rod)
    files.update(reflector_lattices())
    files["small-moc.toml"] = small_moc()
    files["dvc.yaml"] = pipeline(files)
    for name, text in files.items():
        (HERE / name).write_text(text)


if __name__ == "__main__":
    main()
