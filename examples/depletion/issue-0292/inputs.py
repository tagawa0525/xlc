"""#0292: 燃料の種類ごとの履歴の効きを核種に分ける入力を書く。

燃料は 3 種類で、どれも 2A の燃料のピンの形状（../issue-0284/lattice-2A.toml。燃料、ギャップの真空、被覆管、減速材）と状態
（燃料、被覆管、減速材 565 K、ホウ素 1300 ppm）。減速材の密度の履歴は、2A-density（../p2/lattice-2A-density.toml）と同じ
0.743 → 0.661 g/cc:

- u: 3.1% のウランのピン 1 本（鏡面反射）
- mox: LWR MOX ピンセル MOX-1 の組成（../lattice-mox1-rings1.toml。ADR-0206）のピン 1 本（鏡面反射）
- gd: 3×3 の中央に Gd のピン（../p9-gd/fuel-2O.toml の 2O の Gd のピンの組成と 8 つの等面積のリングの物質）、周りに u の燃料の
  ピン 8 本。Gd のピンを 1 本で鏡面反射にすると、周りの燃料の中性子の来ない、Gd のピンだけの格子になる（スペクトルが集合体の
  中の Gd のピンから離れる）ので、周りに燃料を置く。Gd と燃料の本数の比 1:8 は、2O の 24:240 より Gd が多い

`python3 inputs.py` は燃料ごとの <燃料>/nominal.toml と <燃料>/density.toml（減速材 0.743 と 0.661 g/cc の燃焼。出力は <燃料>/out/）を
書き、<燃料>/out/nominal-state.toml（ベース計算の履歴の組成から基準条件（0.743）へのブランチ計算で、履歴の輸送と同じ問題。`micro` は追う 31 核種）
を書く。nominal-state は冷却計算の 1 群の断面積（ADR-0238）で輸送しない核種の吸収を見る（untransported_absorption.py）。燃焼の後に
`python3 inputs.py --branches <燃料>` が、2 つの履歴から次を <燃料>/out/ に書く（どれも 0.661 g/cc の 1 状態の
ブランチ計算で、`BRANCH_BURNUP` の燃焼点を解く）:

- from-nominal.toml: ベース計算（0.743）の履歴の組成から（瞬時ブランチ。ADR-0187）。`micro` に下の核種のすべてを書き、核種ごとの σᵢ と
  均質化した数密度 N̄ᵢ を持つ
- history.toml: 0.661 の履歴の組成から（履歴ブランチの組成。k∞ は履歴の k∞ と同じ問題）。同じ `micro`
- tracked.toml: ベース計算の組成のうち、ミクロ燃焼で追う 31 核種（`TRACKED`）だけを 0.661 の履歴の値にした組成（hybrid-tracked.msgpack）から
- only-<核種>.toml: ベース計算の組成のうち、その核種だけを 0.661 の履歴の値にした組成（hybrid-<核種>.msgpack）から

置き換える核種は、輸送する核種（履歴の `transportable`）のうち、2 つの履歴の燃焼の領域の数密度が `BRANCH_BURNUP` のどこかで違うもの。
燃焼度 0 は from-nominal と history だけが解く（2 つの組成が同じ点で、ブランチ計算の設定の違いを見る検査）。
"""

import argparse
import copy
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "history"))
from make_hybrid import hybrid

FUELS = ("u", "mox", "gd")
NOMINAL_DENSITY = 0.743
VARIED_DENSITY = 0.661
# 燃焼の点は ../p2/lattice-2A.toml と同じ 40 点
DEPLETION_BURNUP = [
    0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 18, 20, 22, 24, 26, 28, 30,
    32, 34, 36, 38, 40, 44, 48, 52, 56, 60,
]  # fmt: skip
# 低、中、高の燃焼度。2、6、10 は Gd の残る点（2O の HZP の履歴の k∞ の山は 8.0〜8.5 MWd/kg。../p9-gd/fuel-2O.toml）
BRANCH_BURNUP = [2, 6, 10, 20, 40, 60]
# nominal-state の燃焼度の点（../issue-0279/core-micro/table-d0.toml の 17 点と同じ。輸送しない核種の吸収の比べ方を 2A と揃える）
COOLING_BURNUP = [0, 0.1, 0.5, 1, 2, 4, 6, 8, 10, 14, 20, 26, 32, 40, 48, 56, 60]
# ADR-0035 の元素の輸送する核種 26 と、その輸送する先行核 5（ADR-0216。../issue-0279/core-micro/table-d0.toml の micro）
TRACKED = [
    "U-234", "U-235", "U-236", "U-238", "Np-237", "Np-239", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242",
    "Am-241", "Am-243", "Xe-135", "Sm-147", "Sm-149", "Sm-150", "Sm-151", "Sm-152",
    "Gd-152", "Gd-154", "Gd-155", "Gd-156", "Gd-157", "Gd-158", "Gd-160",
    "Cm-242", "Cm-244", "Eu-153", "I-135", "Pm-148m",
]  # fmt: skip

UO2 = "{ U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02 }"
MOX = (
    "{ U-235 = 5.1050e-05, U-238 = 2.0370e-02, Pu-238 = 4.6690e-05, Pu-239 = 1.4650e-03, Pu-240 = 5.6910e-04, "
    "Pu-241 = 2.7130e-04, Pu-242 = 1.4130e-04, Am-241 = 3.0280e-05, O-16 = 4.5880e-02 }"
)
GD = (
    "{ U-234 = 3.18096e-06, U-235 = 3.90500e-04, U-236 = 1.79300e-06, U-238 = 2.10299e-02, Gd-152 = 3.35960e-06, "
    "Gd-154 = 3.66190e-05, Gd-155 = 2.48606e-04, Gd-156 = 3.43849e-04, Gd-157 = 2.62884e-04, Gd-158 = 4.17255e-04, "
    "Gd-160 = 3.67198e-04, O-16 = 4.53705e-02 }"
)
# 2O の Gd のピンの 8 つの等面積のリングの外半径（../p9-gd/fuel-2O.toml の cell.gd）
GD_RADII = [
    0.1448154688, 0.2048, 0.2508277497, 0.2896309376, 0.3238172324, 0.3547240054, 0.3831457164, 0.4096,
]  # fmt: skip
GD_MATERIALS = [f"gd{i}" for i in range(len(GD_RADII))]

ZIRCALOY = """[material.zircaloy]
temperature = 565

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
"""

# 燃焼の減速材は 2A と 2A-density の入力の数密度（../issue-0284/lattice-2A.toml、../p2/lattice-2A-density.toml）。ブランチ計算は
# ホウ酸水の式で書く（branches の減速材の軸がホウ酸水の物質を求める）。2 つの書き方の差は数密度の相対 3e-6 以下
DEPLETION_WATER = {
    NOMINAL_DENSITY: "{ O-16 = 2.48112e-02, H-1 = 4.96224e-02, B-10 = 1.07070e-05, B-11 = 4.30971e-05 }",
    VARIED_DENSITY: "{ O-16 = 2.20729e-02, H-1 = 4.41459e-02, B-10 = 9.52537e-06, B-11 = 3.83408e-05 }",
}

PIN = """[cell.pin]
zones = [
  {{ material = "{fuel}", radius = 0.4096, rings = 5 }},
  {{ material = "void", radius = 0.418, rings = 1 }},
  {{ material = "zircaloy", radius = 0.475, rings = 1 }},
  {{ material = "water", radius = 0.63, rings = 3 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8
"""


def gd_cell() -> str:
    rings = "".join(
        f'  {{ material = "{m}", radius = {r}, rings = 1 }},\n'
        for m, r in zip(GD_MATERIALS, GD_RADII)
    )
    return (
        "[cell.gd]\nzones = [\n"
        + rings
        + """  { material = "void", radius = 0.418, rings = 1 },
  { material = "zircaloy", radius = 0.475, rings = 1 },
  { material = "water", radius = 0.63, rings = 3 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8
"""
    )


def fuel_materials(fuel: str) -> list[str]:
    """燃焼させる物質（入力の物質の名前）"""
    return {"u": ["fuel"], "mox": ["fuel"], "gd": ["fuel", *GD_MATERIALS]}[fuel]


def geometry(
    fuel: str, depth: int, water: str, library: str = "p2/depletion.msgpack"
) -> str:
    """ライブラリ、離散化、格子、物質、セル。`depth` は入力の置き場の examples/depletion からの深さ（パスの ../ の数）、`library` は
    examples/depletion からのライブラリのパス"""
    up = "../" * depth
    rows = (
        '["pin pin pin", "pin gd  pin", "pin pin pin"]' if fuel == "gd" else '["pin"]'
    )
    materials = [
        "[material.fuel]",
        "temperature = 565",
        f"nuclides = {MOX if fuel == 'mox' else UO2}",
        "",
    ]
    if fuel == "gd":
        for m in GD_MATERIALS:
            materials += [
                f"[material.{m}]",
                "temperature = 565",
                f"nuclides = {GD}",
                "",
            ]
    return "\n".join(
        [
            f'library = "{up}{library}"',
            "ufg-tolerance = 1e-4",
            "pitch = 1.26",
            "azimuthal = 30",
            "spacing = 0.03",
            f"rows = {rows}",
            "",
            *materials,
            ZIRCALOY,
            "[material.water]",
            "temperature = 565",
            water,
            "",
            PIN.format(fuel="fuel"),
            *([gd_cell()] if fuel == "gd" else []),
        ]
    )


def quoted(names: list[str]) -> str:
    return "[" + ", ".join(f'"{n}"' for n in names) + "]"


def depletion_input(fuel: str, density: float, name: str) -> str:
    head = f"# #0292: {fuel} の燃焼（減速材 {density} g/cc）。inputs.py が書く。手で直さない\n"
    return (
        head
        + geometry(fuel, 2, f"nuclides = {DEPLETION_WATER[density]}")
        + "\n".join(
            [
                "[depletion]",
                'chain = "../../chain.msgpack"',
                "power = 38.41",
                f"burnup = {DEPLETION_BURNUP}",
                f"materials = {quoted(fuel_materials(fuel))}",
                f'output = "out/{name}.msgpack"',
                "",
            ]
        )
    )


def branch_input(
    fuel: str,
    name: str,
    history: str,
    burnup: list[float],
    micro: list[str] | None,
    density: float = VARIED_DENSITY,
    library: str = "p2/depletion.msgpack",
    bounds: list[float] | None = None,
    note: str | None = None,
) -> str:
    """減速材 `density` の 1 状態のブランチ計算（examples/depletion からの深さ 3 に置く）。`bounds` は少数群の境界（省くと 2 群）、
    `note` は先頭の注釈（省くと inputs.py --branches が書く入力の注釈）"""
    water = f"water = {{ density = {density}, boron = 1300, b10 = 0.199 }}"
    groups = "1e7, 0.625, 1e-5" if bounds is None else ", ".join(map(repr, bounds))
    lines = [
        note
        or f"# #0292: {fuel} の {density} g/cc のブランチ計算（{history} の組成から）。inputs.py --branches が書く。手で直さない",
        f'few-group = {{ output = "{name}.msgpack", bounds = [{groups}] }}',
    ]
    body = geometry(fuel, 3, water, library)
    tail = [
        "[branches]",
        'boron = { ppm = [1300], materials = ["water"] }',
        f'moderator = {{ temperature = [565], density = [{density}], materials = ["water", "zircaloy"] }}',
        f"fuel = {{ temperature = [565], materials = {quoted(fuel_materials(fuel))} }}",
        "",
        "[history]",
        f'file = "{history}"',
        'chain = "../../../chain.msgpack"',
        f"burnup = {burnup}",
    ]
    if micro is not None:
        tail.append(f"micro = {quoted(micro)}")
    return "\n".join(lines) + "\n" + body + "\n".join(tail) + "\n"


def write_depletion() -> None:
    for fuel in FUELS:
        directory = HERE / fuel
        (directory / "out").mkdir(parents=True, exist_ok=True)
        for name, density in [
            ("nominal", NOMINAL_DENSITY),
            ("density", VARIED_DENSITY),
        ]:
            (directory / f"{name}.toml").write_text(
                depletion_input(fuel, density, name)
            )
        (directory / "out" / "nominal-state.toml").write_text(
            branch_input(
                fuel,
                "nominal-state",
                "nominal.msgpack",
                COOLING_BURNUP,
                TRACKED,
                NOMINAL_DENSITY,
            )
        )


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def differing_nuclides(nominal: dict, varied: dict, burnup: list[float]) -> list[str]:
    """輸送する核種のうち、2 つの履歴の燃焼の領域の数密度が `burnup` の点のどこかで違うもの（履歴の核種の順）"""
    if nominal["nuclides"] != varied["nuclides"]:
        raise SystemExit("the histories list different nuclides")
    if nominal["transportable"] != varied["transportable"]:
        raise SystemExit("the histories transport different nuclides")
    if [p["burnup"] for p in nominal["points"]] != [
        p["burnup"] for p in varied["points"]
    ]:
        raise SystemExit("the histories have different burnup points")
    points = {p["burnup"]: (p, q) for p, q in zip(nominal["points"], varied["points"])}
    missing = [b for b in burnup if b not in points]
    if missing:
        raise SystemExit(f"burnup {missing} are not points of the histories")
    index = {n: i for i, n in enumerate(nominal["nuclides"])}
    transportable = set(nominal["transportable"])
    return [
        n
        for n in nominal["nuclides"]
        if n in transportable
        and any(
            a[index[n]] != b[index[n]]
            for p, q in (points[x] for x in burnup)
            for a, b in zip(p["densities"], q["densities"])
        )
    ]


def write_branches(fuel: str) -> None:
    out = HERE / fuel / "out"
    nominal_path, varied_path = out / "nominal.msgpack", out / "density.msgpack"
    nominal, varied = read(nominal_path), read(varied_path)
    nuclides = differing_nuclides(nominal, varied, BRANCH_BURNUP)
    missing = [n for n in TRACKED if n not in nominal["nuclides"]]
    if missing:
        raise SystemExit(f"the histories have no {missing}")
    fresh = [0, *BRANCH_BURNUP]
    (out / "from-nominal.toml").write_text(
        branch_input(fuel, "from-nominal", "nominal.msgpack", fresh, nuclides)
    )
    (out / "history.toml").write_text(
        branch_input(fuel, "history", "density.msgpack", fresh, nuclides)
    )
    for name, listed in [("tracked", TRACKED)] + [(n, [n]) for n in nuclides]:
        tag = name if name == "tracked" else f"only-{name}"
        mixed = hybrid(copy.deepcopy(nominal), varied, listed)
        (out / f"hybrid-{name}.msgpack").write_bytes(msgpack.packb(mixed))
        (out / f"{tag}.toml").write_text(
            branch_input(fuel, tag, f"hybrid-{name}.msgpack", BRANCH_BURNUP, None)
        )
    (out / "nuclides.txt").write_text("\n".join(nuclides) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--branches", choices=FUELS)
    args = parser.parse_args()
    if args.branches is None:
        write_depletion()
    else:
        write_branches(args.branches)


if __name__ == "__main__":
    main()
