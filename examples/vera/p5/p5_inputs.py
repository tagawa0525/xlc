"""VERA Problem 5-2D（炉心の高さの中央の 2D の炉心、零出力）を xlc で解く入力と段階を書く。

このファイルのある場所に次を書く。寸法と組成の出典は下の定数のコメントにあり、書いた入力のコメントにもコピーする。

- `assembly-<種類>.toml`: 集合体の種類（濃縮度と、Pyrex の棒の本数か D バンクの制御棒）ごとの単一集合体の格子（外周は
  鏡面反射）。炉心の入力（最小の入力と f* の比較を含む）が読む種類だけ。少数群定数（形状関数つき）を
  `lattice/<種類>.msgpack` に書く
- `reflector-I.toml`、`reflector-corner.toml`: 径方向の反射体の I 型と角型の格子（ADR-0095、ADR-0132）。燃料は外周の
  集合体の代表の 3.10%、Pyrex の棒なし（ADR-0084）。ノードごとの少数群定数を `lattice/` に書く。反射体のノードの定数は
  炉心のすべての反射体の位置に、燃料のノードの定数は外周のこの種類の位置に置く
- `reflector-<形>-<種類>.toml`: 外周の位置に代表と違う種類の集合体（5-2D では Pyrex のある 3.10%）があるときの、その種類の
  燃料の反射体の格子。燃料のノードの定数を、その位置に置く（ADR-0153）
- `small-moc.toml`: 3 × 3 の 3.10% の集合体と反射体の輪の 2D の全体を、真空の外周で MOC で解く格子（最小の入力の参照）
- `core-<ケース>.toml` ほか（ADR-0105 の 4 つ）: 全炉心（ADR-0090）の炉心の計算。燃料のノードの定数は、反射体に 1 面で
  面する位置は I 型、隣り合う 2 面で面する位置は角型の反射体の格子の燃料のノードのもの、ほかは単一集合体のもの（ADR-0153）。
  ケースは次のとおり
  - `5A`、`5B`、`5C`: 5A-2D（制御棒なし）、5B-2D（D バンクに AIC）、5C-2D（D バンクに B4C）
  - `5A-Ly`: 5A の L 型の位置の I 型の定数の f* の辺を、もう一方の燃料の面に向けたもの（L 型の近似の大きさの見積もり）
  - `5A-Inext`: 5A の、角型の位置に隣る I 型の位置に、角型の格子の右と上の定数でなく I 型の定数を置いたもの
  - `small-E`、`small-I`: `small-moc.toml` と同じ体系を、角型の位置に隣る位置に角型の格子の右と上の定数を置いて、
    または I 型の定数を置いて解いたもの
  - `lattice-I-equivalence`、`lattice-I-fstar`、`lattice-corner-equivalence`、`lattice-corner-fstar`: 反射体の格子と同じ
    体系を、燃料のノードに格子の燃料のノードの定数（等価性）か単一集合体の定数（f* の方法）を置いて解いたもの
  - `small-step-Lx`、`small-step-Ly`: 階段の外形（L 型、角型、`-`）を持つ小さな炉心（集合体 21 体）。L 型の位置の
    f* の辺の向きを変えた 2 つ
- `dvc.yaml`: 段階（ADR-0159）

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
# OpenMC の参照解と比べるケース → 問題（参照解のあるものだけ。5B-2D と 5C-2D は #0073）
OPENMC = {"5A": "5A", "5A-Ly": "5A", "5A-Inext": "5A"}
# 反射体の定数を作る燃料（外周の集合体の代表。ADR-0084）。I 型と角型の格子の反射体のノードの定数を、炉心のすべての
# 反射体の位置に置く。燃料のノードの定数は、外周のこの種類の位置に置く（ADR-0153）
REPRESENTATIVE = "31"
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

# バッフルは厚さ 2.85 cm の 1 列（行）のセルにし、セルの中を 3 × 3 に分ける。構造材の塊（ADR-0086）はセルの中で
# まとめるので、厚さを 3 列に分けると中の列のセルが SS304 に囲まれ、露出した面のない塊になる（#0085）
[cell.ss]
zones = [{ material = "ss304", rings = 1 }]
grid = 3"""


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


def assembly_parts(
    fuel: str, pyrex: int, rod: str | None
) -> tuple[str, list[str], list[str]]:
    """集合体の種類の説明、物質の名前、セル"""
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
    return what, names, cells


def assembly(name: str, fuel: str, pyrex: int, rod: str | None) -> str:
    what, names, cells = assembly_parts(fuel, pyrex, rod)
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
# VERA Problem 5-2D の集合体（{what}）の単一集合体の格子（外周は鏡面反射）。燃料、被覆管、減速材は 565 K、
# 減速材 0.743 g/cc、ホウ素 1300 ppm（表 M3-1、M3-4）。ほかは examples/vera/p2 の 2A、2E、2G、2H と同じ形。少数群定数は
# 形状関数と角の不連続因子つき（ADR-0096）で、炉心の反射体に面しない位置の燃料のノードの定数にする（ADR-0153）。
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
    SS304 のバッフル 2.85 cm、水 18.46 cm（0.71 × 26））。仕様書の 1.13 節の表 16。
    `side` は "fuel"、"after"（燃料の右か下の反射体）、"before"（左か上）"""
    if side == "fuel":
        return ["{ size = 0.04 }", "{ size = 1.26, count = 17 }", "{ size = 0.04 }"]
    after = [
        "{ size = 0.19 }",
        "{ size = 2.85 }",
        "{ size = 0.71, count = 26 }",
    ]
    return after if side == "after" else after[::-1]


def cell_marks(side: str) -> list[str]:
    """`columns_of` のセルの並びの種類: 余白 m、ピン p、隙間 g、バッフル b、外の水 w"""
    if side == "fuel":
        return ["m"] + ["p"] * 17 + ["m"]
    after = ["g", "b"] + ["w"] * 26
    return after if side == "after" else after[::-1]


def core_lattice_rows(
    fuel_rows: list[list[str]], kinds: list[list[str]]
) -> list[list[str]]:
    """集合体と反射体を並べた格子のセルの名前。`kinds` は位置ごとの "F"（燃料）、"R"（反射体）で、行は上から。
    バッフルは燃料の外周から隙間 0.19 cm をおいた帯で、角では帯が直角に会う（openmc_core_2d.py の、集合体の外周を
    正方形で太らせた形と同じ）。燃料の行と列の両方にある反射体の位置（L 型）は、列の幅と行の高さが燃料のピンの
    ものになり、バッフルを仕様の寸法で置けないのでエラーにする（#0068）"""
    n_rows, n_cols = len(kinds), len(kinds[0])
    for i in range(n_rows):
        for j in range(n_cols):
            if (
                kinds[i][j] == "R"
                and any(kinds[i][c] == "F" for c in range(n_cols))
                and any(kinds[r][j] == "F" for r in range(n_rows))
            ):
                raise SystemExit(
                    f"reflector at {(i, j)} is in a fuel row and a fuel column (#0068)"
                )
    # 反射体の列（行）は、燃料の片側にだけ並ぶこと（燃料に挟まれた反射体の列には、隙間とバッフルを両側に置けない）
    for line in [*kinds, *zip(*kinds)]:
        fuel_at = [k for k, mark in enumerate(line) if mark == "F"]
        if fuel_at and any(line[k] != "F" for k in range(fuel_at[0], fuel_at[-1] + 1)):
            raise SystemExit(f"reflector between fuel in {''.join(line)}")
    col_sides = side_list(
        [any(kinds[i][j] == "F" for i in range(n_rows)) for j in range(n_cols)]
    )
    row_sides = side_list(
        [any(kinds[i][j] == "F" for j in range(n_cols)) for i in range(n_rows)]
    )
    # セルの中心の座標（集合体のピッチを 1 とする）と、燃料の正方形の和からの距離で、隙間とバッフルを決める
    widths = {"m": 0.04, "p": 1.26, "g": 0.19, "b": 2.85, "w": 0.71}

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
    kinds: list[list[str]],
    boundary: str,
    few_group: str,
    comment: str,
    extra: str = "",
    design: str = REPRESENTATIVE,
) -> str:
    """集合体の種類 `design` と反射体を `kinds` のとおり並べた格子"""
    n_rows, n_cols = len(kinds), len(kinds[0])
    col_fuel = [any(kinds[i][j] == "F" for i in range(n_rows)) for j in range(n_cols)]
    row_fuel = [any(kinds[i][j] == "F" for j in range(n_cols)) for i in range(n_rows)]
    fuel, pyrex, rod = assembly_types()[design]
    rows = core_lattice_rows(pin_rows(fuel, pyrex, rod), kinds)
    _, names, cells = assembly_parts(fuel, pyrex, rod)
    if "ss304" not in names:
        names.append("ss304")
    cells.append(REFLECTOR_CELLS)
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


def reflector_stem(shape: str, design: str) -> str:
    """反射体の格子の名前（`reflector-<形>`、代表の燃料でなければ `reflector-<形>-<集合体の種類>`）"""
    return f"reflector-{shape}" + ("" if design == REPRESENTATIVE else f"-{design}")


# 反射体の格子の、形ごとのノード（名前の末尾、ノードの位置の列と行の番号、反射体か）。名前は格子の名前に続ける
LATTICE_NODES = {
    "I": [("-fuel", 0, 0, False), ("", 1, 0, True)],
    "corner": [
        ("-fuel", 0, 1, False),
        ("-right", 1, 1, True),
        ("-top", 0, 0, True),
        ("", 1, 0, True),
    ],
}


def reflector_lattice(shape: str, design: str) -> str:
    fuel, pyrex, rod = assembly_types()[design]
    what = assembly_parts(fuel, pyrex, rod)[0] + (
        "" if pyrex or rod else "、Pyrex の棒なし"
    )
    stem = reflector_stem(shape, design)
    common = (
        f"# 燃料は{what}。燃料のノードの定数を、外周のこの種類の位置に置く（ADR-0153）。\n"
        "# 反射体は燃料側から、集合体の外周からの隙間の水 0.19 cm、SS304 のバッフル 2.85 cm、水 18.46 cm（仕様書の 1.13 節の\n"
        "# 表 16。集合体のピッチ 21.50 cm の幅）。ノードごとの形状関数と角の不連続因子を付ける（ADR-0131）。"
    )
    if design == REPRESENTATIVE:
        common += "\n# 反射体のノードの定数は、炉心のすべての反射体の位置に置く（ADR-0084、ADR-0153）。"
    if shape == "I":
        kinds = [["F", "R"]]
        cols, rws = node_ranges(["fuel", "after"]), node_ranges(["fuel"])
        boundary = '{ right = "vacuum" }'
        title = (
            "# 径方向の反射体の I 型の格子（ADR-0095）: VERA Problem 5-2D の燃料の集合体の右に反射体を置く。左、上、下は鏡面\n"
            "# 反射、右は真空（examples/reflector/reflector-I.toml と同じ形）。\n"
        )
    else:
        kinds = [["R", "R"], ["F", "R"]]
        cols, rws = node_ranges(["fuel", "after"]), node_ranges(["before", "fuel"])
        boundary = '{ right = "vacuum", top = "vacuum" }'
        title = (
            "# 径方向の反射体の角型の格子（ADR-0095、ADR-0132）: VERA Problem 5-2D の燃料の集合体の右と上に反射体を置き、バッフルは\n"
            "# 集合体の外側の角で L 字に折る。左と下は鏡面反射、右と上は真空（examples/reflector/reflector-corner.toml と同じ形）。\n"
            "# ノードは集合体、右と上の反射体（燃料に面で接する）、右上の反射体（燃料に角の点だけで接する）。\n"
        )
    fg = (
        f"few-group = {{ bounds = {BOUNDS}, nodes = [\n"
        + "".join(
            f'  {{ output = "lattice/{stem}{n}.msgpack", columns = [{cols[c][0]}, {cols[c][1]}], rows = [{rws[r][0]}, {rws[r][1]}]'
            + (", reflector = true" if refl else "")
            + " },\n"
            for n, c, r, refl in LATTICE_NODES[shape]
        )
        + "], pins = true }\n"
    )
    return block_lattice(kinds, boundary, fg, title + common, design=design)


def reflector_lattices(lattices: list[tuple[str, str]]) -> dict[str, str]:
    """反射体の格子: 代表の燃料（ADR-0084）の I 型と角型と、外周のほかの種類の燃料の、その位置の形の格子（ADR-0153）"""
    return {
        f"{reflector_stem(shape, design)}.toml": reflector_lattice(shape, design)
        for shape, design in lattices
    }


def square_moc(name: str) -> str:
    n = SQUARES[name]
    return block_lattice(
        square_kinds(n),
        '{ left = "vacuum", right = "vacuum", top = "vacuum", bottom = "vacuum" }',
        "",
        f"# 正方形の炉心の参照: 3.10% の集合体（Pyrex の棒なし）{n} × {n} と反射体の輪の 2D の全体を、外周を真空にして\n"
        "# MOC で解く。反射体の寸法とバッフルの形は reflector-I.toml、reflector-corner.toml と同じ。\n"
        f"# core の {name}-E、{name}-I と比べる（セルごとの核分裂率）。",
        f'pin-powers = "{name}-moc-pins.msgpack"\n',
    )


# 正方形の炉心（集合体 n × n と反射体の輪）: 最小の入力（3 × 3）。5 × 5 は自己遮蔽で 82 GB を超えて止まった（#0086）
SQUARES = {"small": 3}


def square_kinds(n: int) -> list[list[str]]:
    """F は燃料、R は反射体"""
    return [
        ["F" if 0 < i <= n and 0 < j <= n else "R" for j in range(n + 2)]
        for i in range(n + 2)
    ]


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
                if (faces[0] - faces[1]) % 2 == 0:
                    raise SystemExit(
                        f"reflector at {(i, j)} lies between fuel on opposite sides"
                    )
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


def case(
    name: str,
    comment: str,
    outline: bool = True,
    xy: str = 'x = ["vacuum", "vacuum"], y = ["vacuum", "vacuum"]',
) -> str:
    """計算のケース。`outline` は炉心の形に `-` があること（外形の面の条件を書く。ADR-0132）、`xy` は x と y の境界"""
    face = ', outline = "vacuum"' if outline else ""
    return f"""# このファイルは p5_inputs.py が書く。手で直さない。
{comment}
cycle = "core-{name}-cycle.toml"
axial-mesh = [0.0, {HEIGHT}]
radial-nodes = 1
pin-powers = "core-{name}-pins.msgpack"
boundary = {{ {xy}, z = ["reflective", "reflective"]{face} }}
"""


def cycle(plant_file: str, inventory: str, loading: list[list[str]]) -> str:
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


def periphery(fuel: set[tuple[int, int]]) -> dict[tuple[int, int], tuple[str, int]]:
    """反射体に面する燃料の位置（中心が (0, 0)、行は下が正）ごとの、反射体の格子の形と回転（反時計回りの 90° の回数）。
    1 面で面する位置は I 型（格子の燃料は +x に反射体）、隣り合う 2 面で面する位置は角型（+x と +y に反射体）。
    反射体の位置に角の点だけで接する燃料と、反射体に面しない燃料は含めない（単一集合体の定数。ADR-0153）"""
    out = {}
    for i, j in fuel:
        faces = [d for d, (di, dj) in enumerate(STEPS) if (i + di, j + dj) not in fuel]
        if len(faces) == 1:
            out[(i, j)] = ("I", faces[0])
        elif len(faces) == 2 and (faces[1] - faces[0]) % 4 in (1, 3):
            first = faces[0] if (faces[1] - faces[0]) % 4 == 1 else faces[1]
            (ai, aj), (bi, bj) = STEPS[first], STEPS[(first + 1) % 4]
            # 角型の格子は 2 面の間の対角も反射体なので、対角が燃料の位置（凹んだ外形）には置けない
            if (i + ai + bi, j + aj + bj) in fuel:
                raise SystemExit(
                    f"fuel at {(i, j)} faces the reflector on sides {faces} around fuel"
                )
            out[(i, j)] = ("corner", first)
        elif faces:
            raise SystemExit(f"fuel at {(i, j)} faces the reflector on sides {faces}")
    return out


def site_design(design: str, shape: str | None) -> str:
    """台帳の設計の名前: 反射体に面しない位置は集合体の種類、面する位置は `<種類>-<形>`"""
    return design if shape is None else f"{design}-{shape}"


def inventory(designs: dict[tuple[int, int], str]) -> str:
    """台帳。反射体に面する位置の設計は、その種類の燃料とその位置の形の反射体の格子の燃料のノードの定数（ADR-0153）"""
    edge = periphery(set(designs))
    sites = {
        p: site_design(d, edge[p][0] if p in edge else None) for p, d in designs.items()
    }
    layers = {}
    for p, d in designs.items():
        if p in edge:
            layers[sites[p]] = f"lattice/{reflector_stem(edge[p][0], d)}-fuel.msgpack"
        else:
            layers[d] = f"lattice/{d}.msgpack"
    text = "# このファイルは p5_inputs.py が書く。手で直さない。\n"
    text += (
        "# 集合体の台帳。ID は全炉心の位置（列の名前と行の番号）。設計は、反射体に面しない位置は集合体の種類の単一集合体の\n"
        "# 格子の定数、反射体に 1 面で面する位置（<種類>-I）と 2 面で面する位置（<種類>-corner）は、その種類の燃料の I 型と角型の\n"
        "# 反射体の格子の燃料のノードの定数（ADR-0153。燃料の向きは装荷で回す）\n"
    )
    text += "".join(
        f'\n[design.{d}]\nbounds = [0.0, {HEIGHT}]\nlayers = ["{layers[d]}"]\n'
        for d in sorted(layers)
    )
    text += "".join(
        f'\n[assembly.{assembly_id(i, j)}]\ndesign = "{d}"\n'
        for (i, j), d in sorted(sites.items())
    )
    return text


def loading_of(designs: dict[tuple[int, int], str], size: int) -> list[list[str]]:
    """装荷。反射体に面する位置の集合体は、反射体の格子の燃料の反射体の向きを、その位置の反射体の向きに回す"""
    edge = periphery(set(designs))
    h = size // 2
    return [
        [
            turned(assembly_id(i, j), edge[(i, j)][1] if (i, j) in edge else 0)
            if (i, j) in designs
            else "."
            for j in range(-h, h + 1)
        ]
        for i in range(-h, h + 1)
    ]


def used_lattices(files: dict[str, str]) -> tuple[list[str], list[tuple[str, str]]]:
    """炉心の入力（炉心の形と台帳の `layers`）が読む、単一集合体の格子の種類と、反射体の格子の（形, 集合体の種類）"""
    used = {f for text in files.values() for f in layer_files(text)}
    singles = sorted(d for d in assembly_types() if f"lattice/{d}.msgpack" in used)
    lattices = sorted(
        (shape, d)
        for shape, nodes in LATTICE_NODES.items()
        for d in assembly_types()
        if any(
            f"lattice/{reflector_stem(shape, d)}{n}.msgpack" in used for n, *_ in nodes
        )
    )
    made = {f"lattice/{d}.msgpack" for d in singles} | {
        f"lattice/{reflector_stem(shape, d)}{n}.msgpack"
        for shape, d in lattices
        for n, *_ in LATTICE_NODES[shape]
    }
    if missing := sorted(used - made):
        raise SystemExit(
            f"the core inputs read lattices that no stage writes: {missing}"
        )
    return singles, lattices


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
            "aic": "D バンク（H8、D8、H12、D12 と鏡映の位置の 9 体）に AIC の制御棒",
            "b4c": "D バンク（H8、D8、H12、D12 と鏡映の位置の 9 体）に B4C の制御棒",
        }[rod]
        variants = [("", "")] + (
            [("-Ly", "-Ly"), ("-Inext", "-Inext")] if name == "5A" else []
        )
        for v, p in variants:
            files[f"core-{name}{v}-cycle.toml"] = cycle(
                f"core-plant{p}.toml",
                f"core-{name}-assemblies.toml",
                loading_of(designs, size),
            )
            files[f"core-{name}{v}.toml"] = case(
                f"{name}{v}",
                f"# VERA Problem {name}-2D（{what}）の全炉心の 2D の計算。外周の直方体と外形の面は真空（ADR-0132）、z は鏡面反射。\n"
                "# 制御棒の位置はまだ入力にない（ADR-0105 の `rods`）ので、制御棒の入った集合体は台帳で別の設計にする。",
            )
    # 正方形の炉心: 集合体 n × n と反射体の輪
    for square, n in SQUARES.items():
        h = n // 2 + 1
        fuel_sites = {(i, j) for i in range(-h + 1, h) for j in range(-h + 1, h)}
        designs = {p: "31" for p in fuel_sites}
        files[f"core-{square}-assemblies.toml"] = inventory(designs)
        for nxt in ["E", "I"]:
            files[f"core-{square}-{nxt}-plant.toml"] = plant(
                reflector_map(fuel_sites, n + 2, "x", nxt),
                f"# 正方形の炉心（{square}-moc.toml の体系）の炉心の形。角型の位置に隣る I 型の位置は、"
                + ("角型の格子の右と上の反射体" if nxt == "E" else "I 型の格子の反射体")
                + "。",
            )
            files[f"core-{square}-{nxt}-cycle.toml"] = cycle(
                f"core-{square}-{nxt}-plant.toml",
                f"core-{square}-assemblies.toml",
                loading_of(designs, n + 2),
            )
            files[f"core-{square}-{nxt}.toml"] = case(
                f"{square}-{nxt}",
                f"# 正方形の炉心（{square}-moc.toml の体系）を core で解く。外周の直方体は真空、z は鏡面反射。",
                outline=False,
            )
    # 反射体の格子そのもの: 燃料のノードに、格子のノードの定数（等価性）か単一集合体の定数（f* の方法）を置く
    designs = {
        "reflector-I-fuel": "I 型の格子の燃料のノード",
        "reflector-corner-fuel": "角型の格子の燃料のノード",
        "31": "単一集合体",
    }
    files["core-lattice-assemblies.toml"] = (
        "# このファイルは p5_inputs.py が書く。手で直さない。\n"
        "# 反射体の格子の燃料のノードの定数と、単一集合体の定数の台帳\n"
        + "".join(
            f'\n# {what}\n[design.{d}]\nbounds = [0.0, {HEIGHT}]\nlayers = ["lattice/{d}.msgpack"]\n\n[assembly.{d}]\ndesign = "{d}"\n'
            for d, what in designs.items()
        )
    )
    shapes = {
        "I": (
            [["*", "I"]],
            'x = ["reflective", "vacuum"], y = ["reflective", "reflective"]',
        ),
        "corner": (
            [["N", "C"], ["*", "E"]],
            'x = ["reflective", "vacuum"], y = ["reflective", "vacuum"]',
        ),
    }
    for shape, (map_rows, xy) in shapes.items():
        files[f"core-lattice-{shape}-plant.toml"] = plant(
            map_rows,
            f"# reflector-{shape}.toml の格子をノードに分けた炉心の形（燃料の位置と、格子の反射体のノード）。",
        )
        for kind, fuel in (("equivalence", f"reflector-{shape}-fuel"), ("fstar", "31")):
            name = f"lattice-{shape}-{kind}"
            loading = [["." for _ in row] for row in map_rows]
            loading[-1][0] = fuel
            files[f"core-{name}-cycle.toml"] = cycle(
                f"core-lattice-{shape}-plant.toml",
                "core-lattice-assemblies.toml",
                loading,
            )
            files[f"core-{name}.toml"] = case(
                name,
                f"# reflector-{shape}.toml の格子と同じ体系を、燃料のノードに"
                + (
                    "格子の燃料のノードの定数（等価性）"
                    if kind == "equivalence"
                    else "単一集合体の定数（f* の方法）"
                )
                + "を置いて解く。",
                outline=False,
                xy=xy,
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
            f"core-{name}-plant.toml",
            "core-small-step-assemblies.toml",
            loading_of(step_designs, 7),
        )
        files[f"core-{name}.toml"] = case(
            name, "# 階段の外形の小さな炉心を core で解く（外形と対称性の確認）。"
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


def layer_files(text: str) -> list[str]:
    """`layers = ["a", "b"]` の行が並べる、すべての層のファイル"""
    return sorted(
        {
            f
            for line in text.splitlines()
            if line.startswith("layers")
            for f in line.split('"')[1::2]
        }
    )


def pipeline(
    files: dict[str, str], singles: list[str], lattices: list[tuple[str, str]]
) -> str:
    stages = []
    for name in singles:
        stem = f"assembly-{name}"
        stages.append(
            stage(
                stem,
                f"mkdir -p lattice && ../../../bin/lattice {stem}.toml > {stem}.out",
                [f"{stem}.toml", "../vera.msgpack", *code("lattice")],
                [f"lattice/{name}.msgpack", f"{stem}.out"],
            )
        )
    for shape, design in lattices:
        stem = reflector_stem(shape, design)
        stages.append(
            stage(
                stem,
                f"mkdir -p lattice && ../../../bin/lattice {stem}.toml > {stem}.out",
                [f"{stem}.toml", "../vera.msgpack", *code("lattice")],
                [
                    *(f"lattice/{stem}{n}.msgpack" for n, *_ in LATTICE_NODES[shape]),
                    f"{stem}.out",
                ],
            )
        )
    for square in SQUARES:
        stages.append(
            stage(
                f"{square}-moc",
                f"../../../bin/lattice {square}-moc.toml > {square}-moc.out",
                [f"{square}-moc.toml", "../vera.msgpack", *code("lattice")],
                [f"{square}-moc-pins.msgpack", f"{square}-moc.out"],
            )
        )
    cases = [
        f[5:-5]
        for f in files
        if f.startswith("core-")
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
                f"RUST_LOG=debug ../../../bin/core core-{name}.toml > core-{name}.out 2> core-{name}.log",
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
    for name in [f"{s}-{v}" for s in SQUARES for v in "EI"]:
        moc = name.split("-")[0] + "-moc"
        stages.append(
            stage(
                f"compare-{name}",
                f"python3 ../compare_p5.py moc {moc}.out {moc}-pins.msgpack core-{name}.out core-{name}-pins.msgpack > compare-{name}.out",
                [
                    "../compare_p5.py",
                    "../compare_k.py",
                    f"{moc}.out",
                    f"{moc}-pins.msgpack",
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
                [
                    "../compare_p5.py",
                    "../compare_k.py",
                    f"core-{name}.out",
                    f"core-{name}-pins.msgpack",
                ],
                [f"compare-{name}.out"],
            )
        )
    for shape in ["I", "corner"]:
        for kind in ["equivalence", "fstar"]:
            name = f"lattice-{shape}-{kind}"
            stages.append(
                stage(
                    f"compare-{name}",
                    f"python3 ../compare_p5.py k reflector-{shape}.out core-{name}.out > compare-{name}.out",
                    [
                        "../compare_p5.py",
                        "../compare_k.py",
                        f"reflector-{shape}.out",
                        f"core-{name}.out",
                    ],
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
                    "../compare_k.py",
                    "../../../references/vera/spec.txt",
                    f"core-{name}.out",
                    f"core-{name}-pins.msgpack",
                ],
                [f"compare-{name}.out"],
            )
        )
    # OpenMC の参照解（../dvc.yaml の openmc-<問題>-2D-extract の出力。examples/vera/verification.md の「Problem 5-2D の参照解」）
    for name, problem in OPENMC.items():
        results = f"openmc-{problem}-2D-results"
        stages.append(
            stage(
                f"compare-{name}-openmc",
                f"python3 ../compare_p5.py openmc core-{name}.out core-{name}-pins.msgpack {results} > compare-{name}-openmc.out",
                [
                    "../compare_p5.py",
                    "../compare_k.py",
                    results,
                    f"core-{name}.out",
                    f"core-{name}-pins.msgpack",
                ],
                [f"compare-{name}-openmc.out"],
            )
        )
    return (
        "# このファイルは p5_inputs.py が書く。手で直さない。\n"
        "# VERA Problem 5-2D を xlc で解く段階: 多群ライブラリ（../dvc.yaml の xs）→ 単一集合体の格子（集合体の種類ごと）、\n"
        "# 反射体の I 型と角型の格子、最小の入力の MOC → 炉心（5A、5B、5C と、近似の見積もりと最小の入力の版）→ 比較（KENO-VI、\n"
        "# OpenMC、最小の入力の MOC、対称性）。実行するのは `dvc repro`（ADR-0159）。deps の書き方は examples/m0/dvc.yaml と同じ\n"
        "stages:\n" + "".join(stages)
    )


def main() -> None:
    files = core_files()
    # 格子は、炉心の入力が読むものだけを作る（5-2D の 31-p8 は外周にだけあり、単一集合体の格子は使わない）
    singles, lattices = used_lattices(files)
    for name in singles:
        files[f"assembly-{name}.toml"] = assembly(name, *assembly_types()[name])
    files.update(reflector_lattices(lattices))
    for square in SQUARES:
        files[f"{square}-moc.toml"] = square_moc(square)
    files["dvc.yaml"] = pipeline(files, singles, lattices)
    for name, text in files.items():
        (HERE / name).write_text(text)


if __name__ == "__main__":
    main()
