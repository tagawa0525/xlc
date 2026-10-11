"""VERA Problem 3A の軸方向の反射体の平板（ADR-0128）の格子計算と炉心の計算の入力を書く。

集合体の高さ全体（下の水の外端から上の水の外端まで、518.937 cm）を格子の x（左が炉心の −z）に並べた1行の平板の
入力（`axial-3A-<版>.toml`。少数群定数はディレクトリ `axial-<版>/` に書く）と、それを 1 × 1 × N のノードで解く炉心の
入力（`core-axial-*.toml`）を、このファイルのある場所に書く。層の寸法と組成の出典と導出は下の `LAYERS` と
`mixtures` のコメントにあり、書いた入力のコメントにも写す。

版:
- `16`: 燃料を 16 ノード（22.86 cm）に分け、列の幅を 1 cm 程度にする（`WIDTH`）
- `32`: 燃料を 32 ノード（11.43 cm）に分ける（列は `16` と同じ。f* の燃料のノードの高さへの依存）
- `fine`: `16` の列を半分の幅にする（MOC の平坦線源の近似の収束）

使い方: python3 axial_slab.py
"""

import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent

# 列の幅の目安 (cm)。格子計算の CMFD の粗メッシュはセルごとで、平板では列の数だけの1次元の問題になる。その固有値問題を
# べき乗法の反復ごとに Gauss-Seidel を1回だけ行って解くため、幅 0.25 cm（2058 列）では 10⁵ 回で収束しなかった（#0069）
WIDTH = 1.0

# 格子の行の高さ (cm)。炉心の計算はピッチを同じ値にして、ノード平均の中性子束の規格化（全体の生成率 1）を平板の格子と
# そろえる（格子は奥行き 1 cm あたり、炉心はピッチ × ピッチの断面）
ROW = "1.0"

# 物理的な層（−z から）: (名前, 厚さ (cm), セルの名前)。高さは下の炉心板の上面（集合体の座面）からの cm で、仕様書の
# 1.1 節の表 1、1.3 節の表 3、1.4 節の表 4 と表 5。燃料棒は 6.053 + 4.228 = 10.281 cm から 385.1 cm、端栓は 1.67 cm、
# 燃料は 11.951〜377.711 cm、プレナムは 16.0 cm（377.711〜393.711 cm）、上の端の支持格子（Inconel、高さ 3.866 cm）は
# 中心 388.2 cm で 386.267〜390.133 cm、肩のギャップは 395.381〜397.51 cm（2.129 cm）
LAYERS = [
    ("lower water", "50.0", "water"),
    ("lower core plate", "5.0", "plate"),
    ("bottom nozzle", "6.053", "bnozzle"),
    ("lower gap", "4.228", "gap"),
    ("lower end plug", "1.67", "plug"),
    ("fuel", "365.76", "fuel"),
    ("plenum", "8.556", "plenum"),
    ("plenum with the top end grid", "3.866", "pgrid"),
    ("plenum", "3.578", "plenum"),
    ("upper end plug", "1.67", "plug"),
    ("shoulder gap", "2.129", "gap"),
    ("top nozzle", "8.827", "tnozzle"),
    ("upper core plate", "7.6", "plate"),
    ("upper water", "50.0", "water"),
]

# ノード（core の層）: (名前, 物理的な層の番号)。燃料（層 5）は版ごとに分け、ほかは反射体のノード
NODES = [
    ("lower-water", [0]),
    ("lower-plate", [1]),
    ("bottom-nozzle", [2]),
    ("lower-gap-plug", [3, 4]),
    ("fuel", [5]),
    ("plenum", [6, 7, 8]),
    ("upper-plug-gap", [9, 10]),
    ("top-nozzle", [11]),
    ("upper-plate", [12]),
    ("upper-water", [13]),
]


def areas() -> dict[str, float]:
    """集合体の断面（ピッチ 21.50 cm の正方形）の中の、層の物質の面積 (cm²)。仕様書の表 1、表 2、表 3"""
    pitch = 21.50
    gt = 24 * math.pi * (0.602**2 - 0.561**2)  # 案内管 24 本の Zircaloy
    it = math.pi * (0.605**2 - 0.559**2)  # 計装管の Zircaloy
    rods = 264
    return {
        "assembly": pitch * pitch,
        "tubes": gt + it,
        # 端栓は被覆管の外半径の Zircaloy の棒（端栓の半径は仕様書にない。ADR-0128）
        "plugs": rods * math.pi * 0.475**2,
        "clad": rods * math.pi * (0.475**2 - 0.418**2),
        # プレナムの He は真空（ADR-0120）。プレナムのばねは仕様書に質量がなく、入れない
        "plenum gas": rods * math.pi * 0.418**2,
        # 上の端の支持格子の Inconel-718: 1017 g / 8.19 g/cm³ を高さ 3.866 cm で割った断面積（仕様書の表 3 と 2 節）
        "grid": 1017.0 / 8.19 / 3.866,
    }


def fraction(value: float) -> str:
    """体積の割合を小数 8 桁で書く"""
    return f"{value:.8f}"


def mixtures() -> list[tuple[str, str, list[tuple[str, str]]]]:
    """混ぜた物質: (名前, コメント, [(成分, 体積の割合)])。水の割合は 1 − ほかの割合を同じ桁で書き、和をちょうど 1 にする"""
    a = areas()
    total = a["assembly"]

    def with_water(parts: list[tuple[str, float]]) -> list[tuple[str, str]]:
        written = [(name, fraction(v / total)) for name, v in parts]
        rest = Fraction(1) - sum(Fraction(f) for _, f in written)
        return written + [("water", f"{float(rest):.8f}")]

    return [
        (
            "plate",
            "炉心板: SS304 と冷却材を 50% ずつ（仕様書の表 5）",
            [("ss304", "0.5"), ("water", "0.5")],
        ),
        (
            "bnozzle",
            (
                "ボトムノズル: 水の割合は仕様書の表 P3-3 の H-1 の比 0.72077（6.25 kg / 8.00 g/cm³ / (21.5² × 6.053 cm) の鋼の"
                "割合からは 0.72078）"
            ),
            [("ss304", "0.27923"), ("water", "0.72077")],
        ),
        (
            "tnozzle",
            (
                "トップノズル: 水の割合は仕様書の表 P3-3 の H-1 の比 0.80853（6.25 kg / 8.00 g/cm³ / (21.5² × 8.827 cm) の鋼の"
                "割合からも 0.80853）"
            ),
            [("ss304", "0.19147"), ("water", "0.80853")],
        ),
        (
            "gap",
            f"下のギャップと肩のギャップ: 案内管 24 本と計装管の Zircaloy（{a['tubes']:.5f} cm²）と水",
            with_water([("zircaloy", a["tubes"])]),
        ),
        (
            "plug",
            f"端栓: 半径 0.475 cm の Zircaloy の棒 264 本（{a['plugs']:.5f} cm²）、管と水",
            with_water([("zircaloy", a["plugs"] + a["tubes"])]),
        ),
        (
            "plenum",
            f"プレナム: 被覆管（{a['clad']:.5f} cm²）、He を真空（{a['plenum gas']:.5f} cm²）、管と水",
            with_water(
                [("zircaloy", a["clad"] + a["tubes"]), ("void", a["plenum gas"])]
            ),
        ),
        (
            "pgrid",
            f"上の端の支持格子の高さのプレナム: プレナムに Inconel-718（{a['grid']:.5f} cm²）を足し、水を減らす",
            with_water(
                [
                    ("zircaloy", a["clad"] + a["tubes"]),
                    ("void", a["plenum gas"]),
                    ("inconel", a["grid"]),
                ]
            ),
        ),
    ]


def columns(thickness: str, fine: bool) -> tuple[int, Fraction]:
    """厚さを、幅 `WIDTH` に近く小数 6 桁で割り切れる幅の列に分ける（`fine` なら列の数を倍にする）。列の幅の和が
    厚さに丸めの範囲で等しく、炉心の z の幅と比べられる（ADR-0128）"""
    t = Fraction(thickness)
    n0 = max(1, round(float(t) / WIDTH))
    for d in range(50):
        for n in (n0 + d, n0 - d):
            if n >= 1 and (t / n * 10**6).denominator == 1:
                return (2 * n, t / (2 * n)) if fine else (n, t / n)
    raise ValueError(thickness)


def decimal(x: Fraction) -> str:
    """有限小数で書ける分数を、丸めずに 10 進で書く"""
    for digits in range(12):
        if (x * 10**digits).denominator == 1:
            return f"{float(x):.{max(digits, 1)}f}"
    raise ValueError(x)


def fuel_nodes(variant: str) -> int:
    return 32 if variant == "32" else 16


def slab(variant: str) -> tuple[str, list[tuple[str, Fraction, Fraction]]]:
    """平板の格子計算の入力と、ノードの（名前、下端、上端）(cm、下の水の外端から)"""
    fine = variant == "fine"
    n_fuel = fuel_nodes(variant)
    # 物理的な層を、燃料はノードの高さごとに分けて、列に分ける
    pieces = []  # (層の番号, 厚さ, セル)
    for i, (_, thickness, cell) in enumerate(LAYERS):
        if cell == "fuel":
            each = Fraction(thickness) / n_fuel
            pieces += [(i, each, cell)] * n_fuel
        else:
            pieces.append((i, Fraction(thickness), cell))
    widths, cells, piece_columns = [], [], []
    start = 1
    for _, thickness, cell in pieces:
        n, w = columns(decimal(thickness), fine)
        widths.append(f"{{ size = {decimal(w)}, count = {n} }}")
        cells += [cell] * n
        piece_columns.append((start, start + n - 1, thickness))
        start += n
    # ノードの列の範囲
    nodes = []  # (名前, 最初の列, 最後の列, 厚さ, 反射体か)
    k = 0
    for name, layers in NODES:
        if layers == [5]:
            for f in range(n_fuel):
                first, last, t = piece_columns[k]
                nodes.append((f"fuel-{f + 1:02}", first, last, t, False))
                k += 1
            continue
        first = piece_columns[k][0]
        thickness = Fraction(0)
        for _ in layers:
            thickness += piece_columns[k][2]
            last = piece_columns[k][1]
            k += 1
        nodes.append((name, first, last, thickness, True))

    node_lines = ",\n".join(
        f'  {{ output = "axial-{variant}/{name}.msgpack", columns = [{first}, {last}], rows = [1, 1]'
        + (", reflector = true }" if reflector else " }")
        for name, first, last, _, reflector in nodes
    )
    layer_lines = "\n".join(
        f"#   {name}: {thickness} cm" for name, thickness, _ in LAYERS
    )
    mixture_text = "\n\n".join(
        f"# {comment}\n[material.{name}]\nmix = ["
        + ", ".join(f'{{ material = "{m}", volume = {v} }}' for m, v in parts)
        + "]"
        for name, comment, parts in mixtures()
    )
    text = f"""# このファイルは axial_slab.py が書く（版 {variant}）。手で直さない。
# VERA Problem 3A の軸方向の反射体の平板（ADR-0128）: 集合体の高さ全体（下の水の外端から上の水の外端まで 518.937 cm）を
# 格子の x（左が炉心の −z）に並べた1行の平板。左右は真空、上下は鏡面反射。燃料は単一集合体の格子計算
# （assembly-3A.toml）が書いた媒質。層の物質は、層の断面（ピッチ 21.50 cm の正方形）の中で体積で混ぜる。自己遮蔽は
# しない（真空の辺のある格子の Dancoff 係数は #0008 で未決）。燃料の中の支持格子は入れない（#0011）。
# 物理的な層（−z から。仕様書の 1.1 節の表 1、1.3 節の表 3、1.4 節の表 4、表 5）:
{layer_lines}
# 燃料は {n_fuel} ノード。列の幅は層ごとに {WIDTH} cm {"の半分 " if fine else ""}程度。
library = "reflector.msgpack"
self-shielding = false
pitch = {ROW}
azimuthal = 32
spacing = 0.05
boundary = {{ left = "vacuum", right = "vacuum" }}
few-group = {{ bounds = [1e7, 0.625, 1e-5], axial = true, nodes = [
{node_lines},
] }}
widths = [
  {(",\n  ").join(widths)},
]
heights = [{{ size = {ROW} }}]
rows = ["{" ".join(cells)}"]

[material.fuel]
medium = "medium-3A.msgpack"

{MATERIALS}

{mixture_text}

[cell.fuel]
zones = [{{ material = "fuel", rings = 1 }}]
""" + "\n".join(
        f'\n[cell.{name}]\nzones = [{{ material = "{name}", rings = 1 }}]'
        for name in [
            "water",
            "plate",
            "bnozzle",
            "gap",
            "plug",
            "plenum",
            "pgrid",
            "tnozzle",
        ]
    )
    bounds = []
    z = Fraction(0)
    for name, _, _, thickness, _ in nodes:
        bounds.append((name, z, z + thickness))
        z += thickness
    return text + "\n", bounds


# 層の物質（600 K）。燃料の集合体（assembly-3A.toml）と同じ Zircaloy-4 と冷却材、reflector-I.toml と同じ SS304
MATERIALS = """# 冷却材（仕様書の表 P3-3 の 3A、1300 ppm。0.743 g/cm³）
[material.water]
temperature = 600
nuclides = { H-1 = 4.96224e-02, O-16 = 2.48112e-02, B-10 = 1.07070e-05, B-11 = 4.30971e-05 }

# Zircaloy-4（仕様書の表 P3-3 の被覆管。Hf は assembly-3A.toml と同じく入れない）
[material.zircaloy]
temperature = 600

[material.zircaloy.nuclides]
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

# SS304（仕様書の表 P3-3。天然の C は C-12 98.93%、C-13 1.07% に分ける。ADR-0079）
[material.ss304]
temperature = 600

[material.ss304.nuclides]
C-12 = 3.17461e-04
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
Ni-64 = 7.21770e-05

# Inconel-718（仕様書の表 P3-3）
[material.inconel]
temperature = 600

[material.inconel.nuclides]
Si-28 = 4.04885e-03
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


def core_files(
    variant: str, bounds: list[tuple[str, Fraction, Fraction]]
) -> dict[str, str]:
    """炉心の計算の入力: 平板の定数だけの設計（等価性）と、燃料を単一集合体の定数にした設計（f* の方法）"""
    mesh = ", ".join(decimal(b[1]) for b in bounds) + f", {decimal(bounds[-1][2])}"
    height = decimal(bounds[-1][2])
    files = {}
    designs = []
    for kind in ["slab", "fstar"]:
        layer_bounds, layers = [], []
        for name, low, high in bounds:
            if kind == "fstar" and name.startswith("fuel"):
                # 燃料の層は単一集合体の定数を1つ（z の厚さを持たず、どのノードにも置ける）
                if layers and layers[-1] == "fg-3A.msgpack":
                    continue
                layers.append("fg-3A.msgpack")
            else:
                layers.append(f"axial-{variant}/{name}.msgpack")
            layer_bounds.append(decimal(low))
        layer_bounds.append(height)
        designs.append(
            f"[design.{kind}-{variant}]\nbounds = [{', '.join(layer_bounds)}]\nlayers = [\n"
            + "".join(f'  "{layer}",\n' for layer in layers)
            + f']\n\n[assembly.{kind.upper()}-{variant}]\ndesign = "{kind}-{variant}"\n'
        )
        what = (
            "定数を作った体系の等価性（燃料のノードも平板の定数）"
            if kind == "slab"
            else "f* の方法（燃料は単一集合体の定数 fg-3A.msgpack）"
        )
        files[f"core-axial-{kind}-{variant}-cycle.toml"] = (
            f"# このファイルは axial_slab.py が書く。{what}の装荷\n"
            'plant = "core-axial-plant.toml"\n'
            f'inventory = "core-axial-assemblies-{variant}.toml"\n'
            f'loading = ["{kind.upper()}-{variant}"]\n'
        )
        files[f"core-axial-{kind}-{variant}.toml"] = (
            f"# このファイルは axial_slab.py が書く。平板 axial-3A-{variant}.toml の {len(bounds)} ノードを 1 × 1 × {len(bounds)} で解く\n"
            f"# （{what}）。径方向は反射、z の両端は真空で、平板の定数の真空の辺の比 J/φ（ADR-0101、ADR-0128）\n"
            f'cycle = "core-axial-{kind}-{variant}-cycle.toml"\n'
            f"axial-mesh = [{mesh}]\n"
            "radial-nodes = 1\n"
            'boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], '
            'z = ["vacuum", "vacuum"] }\n'
        )
    files[f"core-axial-assemblies-{variant}.toml"] = (
        f"# このファイルは axial_slab.py が書く。平板 axial-3A-{variant}.toml の層の設計\n"
        + "\n".join(designs)
    )
    return files


def main() -> None:
    plant_height = None
    for variant in ["16", "32", "fine"]:
        text, bounds = slab(variant)
        (HERE / f"axial-3A-{variant}.toml").write_text(text)
        height = decimal(bounds[-1][2])
        if plant_height not in (None, height):
            raise SystemExit(f"heights differ: {plant_height} vs {height}")
        plant_height = height
        for name, content in core_files(variant, bounds).items():
            (HERE / name).write_text(content)
    (HERE / "core-axial-plant.toml").write_text(
        "# このファイルは axial_slab.py が書く。軸方向の反射体の平板の炉心: 1 体だけ。ピッチは平板の格子の行の高さ\n"
        "# （ノード平均の中性子束の規格化をそろえる）\n"
        f'pitch = {ROW}\nheight = {plant_height}\nmap = ["*"]\n'
    )


if __name__ == "__main__":
    main()
