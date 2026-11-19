"""#0293: 燃焼の履歴の土台になる上流の前提（項目 2、6、7、8）を、ピンセルの小さな計算で見積もる入力を書く。

燃料は 3 種類で、形状は 2A の燃料のピン（../lattice-1C-rings1.toml と同じ寸法、燃料 1 つの燃焼の領域を MOC の 5 リング）。状態は
燃料 900 K（Problem 9 の基準条件の燃料の温度。ADR-0239）、被覆管と減速材 565 K、減速材 0.743 g/cc、ホウ素 1300 ppm（2A の燃焼の
OpenMC の参照の状態。ADR-0178）。減速材と被覆管を基準条件の 580 K にしないのは、ピンセルのライブラリ（../depletion.msgpack）の
S(α,β) を持たない核種の温度が 565、600、900、1200 K だけで、それ以外の温度は内挿できない（xlc_format::library の `lookup`）ため:

- u: 3.1% のウランのピン 1 本（鏡面反射）。ライブラリは ../depletion.msgpack
- mox: LWR MOX ピンセル MOX-1 の組成（../lattice-mox1-rings1.toml。ADR-0206）のピン 1 本。ライブラリは ../depletion.msgpack
- gd: 3×3 の中央に 2O の Gd のピン（../p9-gd/fuel-2O.toml の Gd の組成と 8 つの等面積のリング）、周りに u の燃料のピン 8 本。
  Gd のピン 1 本の鏡面反射では周りの燃料の中性子が来ず、集合体の中の Gd のピンのスペクトルから離れるので周りに燃料を置く（git の
  ブランチ feat/0279-off-history-comparison の issue-0292/gd と同じ形）。ピンセルのライブラリは Gd-155 のほかの Gd の同位体を
  持たないので、ライブラリは ../p2/depletion.msgpack（本番のライブラリ。Pm-147、Pm-148、Pm-149、Am-242m、Eu-151 を輸送しない。#0293 の項目 1）

`python3 inputs.py` は次を書く（出力は各ディレクトリの out/。`.gitignore` で追跡しない）:

- steps/<燃料>-40.toml、steps/<燃料>-80.toml（項目 2）: 燃焼度の点を本番の 40 点（../p9-product/lattice-2A.toml）にした燃焼と、
  隣り合う点の中点（0 と最初の点の中点を含む）を足して倍にした 80 点の燃焼。u と gd はさらに倍の 160 点（steps/<燃料>-160.toml）も
- sab/xs-<a|b>.toml と sab/nuclides-<a|b>/（項目 8）: H-1 の温度だけが違う 2 つのライブラリの入力。a は H-1 を評価済みの
  550、573.6、600 K、b は 550、600 K で作る（b の 573.6 K は 550 と 600 K から内挿する）。どちらも O-16、B-10、B-11 に 573.6 K を
  足す（S(α,β) を持たない核種は内挿できないので、減速材を 573.6 K にするため）。ほかの核種は ../xs-nuclides/ への symlink
- sab/<燃料>-<a|b>.toml（項目 8。u と mox）: 40 点の履歴の組成から、減速材だけを 573.6 K にしたブランチ計算（密度とホウ素は燃焼と同じ。
  被覆管は 565 K のまま）
- urr/<燃料>.toml（項目 6、7。u と mox）: 40 点の履歴の組成の、燃焼と同じ状態のブランチ計算。少数群を WIMS の 69 群のまま書き、`micro` に
  非分離共鳴を持つ核種と核分裂する核種を並べる（群ごとの微視的断面積と、核分裂率の熱の割合 `thermal_fission`）
"""

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

FUELS = ("u", "mox", "gd")
# 本番の燃焼度の点（../p9-product/lattice-2A.toml の 40 点）
BURNUP_40 = [
    0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 18, 20, 22, 24, 26, 28, 30,
    32, 34, 36, 38, 40, 44, 48, 52, 56, 60,
]  # fmt: skip
# ブランチ計算の燃焼点（新しい燃料、中、高の燃焼度）
BRANCH_BURNUP = [0, 20, 40, 60]
FUEL_TEMPERATURE = 900
MODERATOR_TEMPERATURE = 565
# 内挿を確かめる H-1 の評価済みの温度（ENDF/B-VIII.0 の tsl-HinH2O の 550 と 600 K の間）
SAB_TEMPERATURE = 573.6
WATER = "water = { density = 0.743, boron = 1300, b10 = 0.199 }"
POWER = 38.41
# ライブラリの群境界（WIMS の 69 群。../depletion.msgpack の group_bounds）。項目 6、7 の少数群をこのままにして、群ごとの値を書く
WIMS69 = [
    1e7, 6065500, 3679000, 2231000, 1353000, 821000, 500000, 302500, 183000, 111000,
    67340, 40850, 24780, 15030, 9118, 5530, 3519.1, 2239.45, 1425.1, 906.898,
    367.262, 148.728, 75.501, 48.052, 27.7, 15.968, 9.877, 4, 3.3, 2.6,
    2.1, 1.5, 1.3, 1.15, 1.123, 1.097, 1.071, 1.045, 1.02, 0.996,
    0.972, 0.95, 0.91, 0.85, 0.78, 0.625, 0.5, 0.4, 0.35, 0.32,
    0.3, 0.28, 0.25, 0.22, 0.18, 0.14, 0.1, 0.08, 0.067, 0.058,
    0.05, 0.042, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.005, 1e-5,
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

WATER_ZONES = """  { material = "void", radius = 0.418, rings = 1 },
  { material = "zircaloy", radius = 0.475, rings = 1 },
  { material = "water", radius = 0.63, rings = 3 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8
"""

# 項目 6、7 の `micro`: ピンセルのライブラリで輸送する核種のうち、非分離共鳴の範囲を持つものと核分裂するもの。どちらかに当たる
# かは estimates.py が ENDF から確かめ、足りなければ止まる
URR_MICRO = [
    "U-234", "U-235", "U-236", "U-238", "Np-237", "Np-239", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242",
    "Am-241", "Am-242m", "Am-243", "Cm-242", "Cm-244",
    "Mo-95", "Tc-99", "Ru-101", "Rh-103", "Ag-109", "Cs-133", "Nd-143", "Nd-145", "Pm-147", "Pm-148", "Pm-148m",
    "Pm-149", "Sm-147", "Sm-149", "Sm-150", "Sm-151", "Sm-152", "Eu-151", "Eu-153", "Gd-155", "I-135", "Xe-135",
    "Sn-112", "Sn-114", "Sn-115", "Sn-117", "Sn-118", "Sn-119",
]  # fmt: skip


def doubled(points: list[float]) -> list[float]:
    """燃焼度の点の隣り合う 2 点（0 と最初の点を含む）の中点を足して、点の数を倍にする"""
    previous = [0.0, *points[:-1]]
    out: list[float] = []
    for low, high in zip(previous, points):
        out += [round((low + high) / 2, 6), high]
    return out


def library(fuel: str, up: str) -> str:
    return f"{up}p2/depletion.msgpack" if fuel == "gd" else f"{up}depletion.msgpack"


def fuel_materials(fuel: str) -> list[str]:
    """燃焼させる物質（入力の物質の名前）"""
    return ["fuel", *GD_MATERIALS] if fuel == "gd" else ["fuel"]


def quoted(names: list[str]) -> str:
    return "[" + ", ".join(f'"{n}"' for n in names) + "]"


def numbers(values: list[float]) -> str:
    return "[" + ", ".join(f"{v:g}" for v in values) + "]"


def geometry(fuel: str, library_path: str) -> str:
    """ライブラリ、離散化、格子、物質、セル"""
    rows = (
        '["pin pin pin", "pin gd  pin", "pin pin pin"]' if fuel == "gd" else '["pin"]'
    )
    lines = [
        f'library = "{library_path}"',
        # 本番の格子計算の設定（../p9-product/lattice-2A.toml）
        "ufg-tolerance = 1e-4",
        "pitch = 1.26",
        "azimuthal = 30",
        "spacing = 0.03",
        f"rows = {rows}",
        "",
        "[material.fuel]",
        f"temperature = {FUEL_TEMPERATURE}",
        f"nuclides = {MOX if fuel == 'mox' else UO2}",
        "",
    ]
    if fuel == "gd":
        for m in GD_MATERIALS:
            lines += [
                f"[material.{m}]",
                f"temperature = {FUEL_TEMPERATURE}",
                f"nuclides = {GD}",
                "",
            ]
    lines += [
        ZIRCALOY,
        "[material.water]",
        f"temperature = {MODERATOR_TEMPERATURE}",
        WATER,
        "",
        '[cell.pin]\nzones = [\n  { material = "fuel", radius = 0.4096, rings = 5 },\n'
        + WATER_ZONES,
    ]
    if fuel == "gd":
        rings = "".join(
            f'  {{ material = "{m}", radius = {r}, rings = 1 }},\n'
            for m, r in zip(GD_MATERIALS, GD_RADII)
        )
        lines.append("[cell.gd]\nzones = [\n" + rings + WATER_ZONES)
    return "\n".join(lines)


def depletion_input(fuel: str, points: int) -> str:
    burnup = BURNUP_40
    while len(burnup) < points:
        burnup = doubled(burnup)
    assert len(burnup) == points
    return (
        f"# #0293 の項目 2: {fuel} の燃焼（燃焼度の点 {points}）。inputs.py が書く。手で直さない\n"
        + geometry(fuel, library(fuel, "../../"))
        + "\n[depletion]\n"
        + 'chain = "../../chain.msgpack"\n'
        + f"power = {POWER}\n"
        + f"burnup = {numbers(burnup)}\n"
        + f"materials = {quoted(fuel_materials(fuel))}\n"
        + f'output = "out/{fuel}-{points}.msgpack"\n'
    )


def branch_input(
    fuel: str,
    head: str,
    library_path: str,
    moderator: float,
    few_group: str,
    micro: list[str],
) -> str:
    """40 点の履歴（steps/out/<燃料>-40.msgpack）の組成からのブランチ計算。減速材の軸は減速材だけを変える"""
    return (
        head
        + few_group
        + "\n"
        + geometry(fuel, library_path)
        + "\n[branches]\n"
        + 'boron = { ppm = [1300], materials = ["water"] }\n'
        + f'moderator = {{ temperature = [{moderator:g}], density = [0.743], materials = ["water"] }}\n'
        + f"fuel = {{ temperature = [{FUEL_TEMPERATURE}], materials = {quoted(fuel_materials(fuel))} }}\n"
        + "\n[history]\n"
        + f'file = "../steps/out/{fuel}-40.msgpack"\n'
        + 'chain = "../../chain.msgpack"\n'
        + f"burnup = {numbers(BRANCH_BURNUP)}\n"
        + (f"micro = {quoted(micro)}\n" if micro else "")
    )


def library_input(variant: str) -> str:
    """../xs.toml のうち H-1、O-16、B-10、B-11 の温度だけを変えたライブラリの入力"""
    text = (HERE.parent / "xs.toml").read_text()
    h1 = "[550, 573.6, 600]" if variant == "a" else "[550, 600]"
    replacements = [
        ('output = "depletion.msgpack"', f'output = "out/depletion-{variant}.msgpack"'),
        ('[nuclide.H-1]\nendf = "n-001_H_001.endf"\ntemperatures = [550, 573.6, 600]',
         f'[nuclide.H-1]\nendf = "n-001_H_001.endf"\ntemperatures = {h1}'),
        ('[nuclide.O-16]\nendf = "n-008_O_016.endf"\ntemperatures = [565, 600, 900, 1200]',
         '[nuclide.O-16]\nendf = "n-008_O_016.endf"\ntemperatures = [565, 573.6, 600, 900, 1200]'),
        ('[nuclide.B-10]\nendf = "n-005_B_010.endf"\ntemperatures = [565, 600]',
         '[nuclide.B-10]\nendf = "n-005_B_010.endf"\ntemperatures = [565, 573.6, 600]'),
        ('[nuclide.B-11]\nendf = "n-005_B_011.endf"\ntemperatures = [565, 600]',
         '[nuclide.B-11]\nendf = "n-005_B_011.endf"\ntemperatures = [565, 573.6, 600]'),
    ]  # fmt: skip
    # 評価済みデータのディレクトリ（endf-directory と thermal-directory）は、置き場が 2 段深くなる
    text = text.replace('"../../references/', '"../../../../references/')
    for old, new in replacements:
        if text.count(old) != 1:
            sys.exit(f"../xs.toml: expected `{old}` once")
        text = text.replace(old, new)
    head = (
        f"# #0293 の項目 8: ../../xs.toml から inputs.py が書くライブラリ {variant}（H-1 の温度 {h1}。O-16、B-10、B-11 に 573.6 K を"
        " 足す）。手で直さない\n"
    )
    return head + text


# ライブラリ a で作り直す核種。b は H-1 だけを作り直し、O-16、B-10、B-11 は a のものを使う。ほかの核種は ../xs-nuclides/ のもの
REBUILT = ("H-1", "O-16", "B-10", "B-11")


def nuclide_names() -> list[str]:
    return [
        line[len("[nuclide.") : -1]
        for line in (HERE.parent / "xs.toml").read_text().splitlines()
        if line.startswith("[nuclide.")
    ]


def link_target(variant: str, name: str) -> Path | None:
    """ライブラリ `variant` の核種 `name` の symlink の先（作り直す核種は None）"""
    if name == "H-1" or (variant == "a" and name in REBUILT):
        return None
    if name in REBUILT:
        return HERE / "sab" / "nuclides-a" / f"{name}.msgpack"
    return HERE.parent / "xs-nuclides" / f"{name}.msgpack"


def link_nuclides(variant: str) -> None:
    """作り直さない核種の相対の symlink を sab/nuclides-<variant>/ に置く"""
    directory = HERE / "sab" / f"nuclides-{variant}"
    directory.mkdir(parents=True, exist_ok=True)
    for name in nuclide_names():
        target = link_target(variant, name)
        if target is None:
            continue
        link = directory / f"{name}.msgpack"
        relative = os.path.relpath(target, directory)
        if link.is_symlink():
            if os.readlink(link) == relative:
                continue
            link.unlink()
        link.symlink_to(relative)


def write(path: Path, text: str) -> None:
    """入力を書き、出力の置き場 out/（lattice は標準出力とログをシェルがそこに受ける）を作る"""
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.parent / "out").mkdir(exist_ok=True)
    path.write_text(text)


def main() -> None:
    for fuel in FUELS:
        # U は 80 点と 40 点の差が目安の 5 pcm を超えたので、さらに倍の 160 点で 2 次の収束を確かめる。Gd の本番の履歴
        # （../p9-gd/lattice-2O.toml の 72 点）は燃え尽きの間を 80 点と同じ 0.5 MWd/kg で刻むので、その誤差を 160 点で見る
        for points in (40, 80, 160) if fuel in ("u", "gd") else (40, 80):
            write(
                HERE / "steps" / f"{fuel}-{points}.toml", depletion_input(fuel, points)
            )
    for variant in ("a", "b"):
        write(HERE / "sab" / f"xs-{variant}.toml", library_input(variant))
        link_nuclides(variant)
        for fuel in ("u", "mox"):
            head = (
                f"# #0293 の項目 8: {fuel} の 40 点の履歴の組成から、減速材を {SAB_TEMPERATURE} K にしたブランチ計算（ライブラリ {variant}）。"
                "inputs.py が書く。手で直さない\n"
            )
            write(
                HERE / "sab" / f"{fuel}-{variant}.toml",
                branch_input(
                    fuel,
                    head,
                    f"out/depletion-{variant}.msgpack",
                    SAB_TEMPERATURE,
                    f'few-group = {{ output = "out/{fuel}-{variant}.msgpack", bounds = [1e7, 0.625, 1e-5] }}',
                    [],
                ),
            )
    for fuel in ("u", "mox"):
        head = (
            f"# #0293 の項目 6、7: {fuel} の 40 点の履歴の組成の、燃焼と同じ状態のブランチ計算（WIMS の 69 群のまま）。"
            "inputs.py が書く。手で直さない\n"
        )
        write(
            HERE / "urr" / f"{fuel}.toml",
            branch_input(
                fuel,
                head,
                library(fuel, "../../"),
                MODERATOR_TEMPERATURE,
                f'few-group = {{ output = "out/{fuel}.msgpack", bounds = {numbers(WIMS69)} }}',
                URR_MICRO,
            ),
        )


if __name__ == "__main__":
    main()
