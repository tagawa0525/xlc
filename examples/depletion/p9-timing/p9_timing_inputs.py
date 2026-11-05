"""Problem 9 の形の全炉心の燃焼の時間の内訳（#0156。examples/depletion/verification.md の「全炉心の燃焼の時間の内訳（Problem 9 の
形）」）のケースを書く。

使い方: python3 p9_timing_inputs.py [--check]

- 炉心の形と反射体は Problem 5 の 3D の全炉心（../../vera/p5/hzp/core-hzp-plant.toml をそのまま指す）。台帳は p5 の台帳
  （core-hzp-assemblies.toml。p5_inputs.py が書く）の設計ごとの軸方向の層のうち、格子計算の層（`lattice/` の、燃料のある層）を
  2A のミクロ燃焼のテーブル（`TABLE_2A`）に替え、ほかの層（ノズル、プレナム、水など）は p5 のテーブルのまま。装荷は p5 の装荷図
- テーブルの軸の制約: 2A のテーブルの減速材と燃料の温度は 565 K の 1 点で、p5 の反射体と非燃料の層のホウ素の軸は 1100〜1400 ppm
  （2A は 0〜1300 ppm）。このため減速材は表 P9-2 の入口の温度でなくすべて 565 K にし、ホウ素はどの解も 1100〜1300 ppm に
  入るように、ステップ 1 はホウ素を与え（燃焼点 0 は Xe のない解で、平衡 Xe の解と同じ目標の k では 200 ppm の窓に入らない）、
  ステップ 2 から探す。探索の目標の k は、ステップ 2 の状態にホウ素を与えた試行の燃焼点 1 の k（`TARGET_K` の注）
- ステップは表 P9-2 の場合 2 と 3（EFPD、出力の割合、バンク D の位置）。時間を測るのに要る最初の 2 ステップだけを書く
- `--check`: 書く内容がファイルと同じかを確かめ、違えば止める（ファイルを書かない）
"""

import sys
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
P5 = Path("../../vera/p5/hzp")
TABLE_2A = "../branches/table-2A-micro-boron-rod-axis.msgpack"

# 定格の出力密度 (W/g)。ほかの core の燃焼のケースと同じ
SPECIFIC_POWER = 38.41

# ステップ 1 のホウ素濃度 (ppm)。p5 と 2A のテーブルのホウ素の軸の共通の範囲（1100〜1300 ppm）の中
GIVEN_BORON = 1200.0

# 探索の範囲 (ppm)。p5 と 2A のテーブルのホウ素の軸の共通の範囲
SEARCH_RANGE = (1100.0, 1300.0)

# 探索の目標の k: ステップ 2 の状態にホウ素 GIVEN_BORON を与えた試行の燃焼点 1 の k（平衡 Xe）
TARGET_K = None

# 表 P9-2 の場合 2、3: (状態の名前, EFPD, 出力の割合, バンク D の位置 (ステップ))。場合 2、3 の値は端点（表の印 *）で、平均に直していない
STEPS = [
    ("case-2", 9.0, 0.657, 192),
    ("case-3", 32.0, 0.997, 219),
]

HEADER = "# このファイルは p9_timing_inputs.py が書く。手で直さない。\n"


def toml_list(values: list) -> str:
    return (
        "["
        + ", ".join(repr(v) if isinstance(v, float) else f'"{v}"' for v in values)
        + "]"
    )


def assemblies() -> str:
    """p5 の台帳の、格子計算の層を 2A のテーブルに替えた台帳"""
    p5 = tomllib.loads((HERE / P5 / "core-hzp-assemblies.toml").read_text())
    lines = [
        HEADER,
        "# p5 の台帳（../../vera/p5/hzp/core-hzp-assemblies.toml）の格子計算の層（燃料のある層）を 2A のミクロ燃焼のテーブルに替えた\n",
        "# 台帳。濃縮度、可燃性毒物、反射体に面する位置の違いは無視する（時間の内訳を測るため。#0156）\n",
    ]
    for name, design in p5["design"].items():
        layers = [
            TABLE_2A if layer.startswith("lattice/") else str(P5 / layer)
            for layer in design["layers"]
        ]
        lines.append(
            f"\n[design.{name}]\nbounds = {toml_list(design['bounds'])}\nlayers = [\n"
        )
        lines.extend(f'  "{layer}",\n' for layer in layers)
        lines.append("]\n")
    for name, assembly in p5["assembly"].items():
        lines.append(f'\n[assembly.{name}]\ndesign = "{assembly["design"]}"\n')
    return "".join(lines)


def cycle() -> str:
    """p5 の炉心の形と装荷図のサイクル"""
    p5 = tomllib.loads((HERE / P5 / "core-hzp-cycle.toml").read_text())
    rows = "".join(f'  "{row}",\n' for row in p5["loading"])
    return (
        HEADER
        + f'plant = "{P5 / p5["plant"]}"\n'
        + 'inventory = "core-p9-assemblies.toml"\n'
        + f"loading = [\n{rows}]\n"
    )


def state(name: str, bank_d: int, boron: str) -> str:
    banks = ", ".join(
        f"{b} = {bank_d if b == 'D' else 230}"
        for b in ["A", "B", "C", "D", "SA", "SB", "SC", "SD"]
    )
    return f"\n[state.{name}]\n{boron}\nb10 = 0.199\nmoderator = 565\nfuel = 565\nbanks = {{ {banks} }}\n"


def case(target_k: float | None) -> str:
    """燃焼のケース。`target_k` が None なら、ステップ 2 にもホウ素を与える（目標の k を決める試行）"""
    p5 = tomllib.loads((HERE / P5 / "core-hzp.toml").read_text())
    boundary = p5["boundary"]
    sides = ", ".join(
        f'{axis} = ["{boundary[axis][0]}", "{boundary[axis][1]}"]'
        for axis in ["x", "y", "z"]
    )
    given = f"boron = {GIVEN_BORON}"
    low, high = SEARCH_RANGE
    search = (
        given
        if target_k is None
        else f"search = {{ boron = [{low}, {high}], k = {target_k}, tolerance = 1.0e-6, max-solves = 20 }}"
    )
    text = (
        HEADER
        + "# Problem 9 の形の全炉心の燃焼（表 P9-2 の場合 2、3 の 2 ステップ）の時間の内訳を測る（#0156）。炉心と軸方向の\n"
        + "# メッシュは Problem 5 の 3D の全炉心（../../vera/p5/hzp/core-hzp.toml）と同じ。ステップ 1 はホウ素を与え、ステップ 2 は\n"
        + "# 臨界ホウ素を探す（範囲と目標の k は p9_timing_inputs.py の注）。Xe-135 と I-135 は出力のある解で平衡の値（ADR-0219）\n"
        + 'cycle = "core-p9-cycle.toml"\n'
        + f"axial-mesh = {toml_list(p5['axial-mesh'])}\n"
        + "radial-nodes = 1\n"
        + f'boundary = {{ {sides}, outline = "{boundary["outline"]}" }}\n'
        + state(STEPS[0][0], STEPS[0][3], given)
        + state(STEPS[1][0], STEPS[1][3], search)
        + "\n[depletion]\n"
        + f"specific-power = {SPECIFIC_POWER}\n"
        + 'history = "core-p9.msgpack"\n'
        + "# 格子計算の燃焼（../p2/lattice-2A.toml）と同じチェーン\n"
        + 'micro = { chain = "../chain.msgpack", xenon = "equilibrium" }\n'
        + "steps = [\n"
        + "".join(
            f'  {{ efpd = {efpd}, power = {power}, state = "{name}" }},\n'
            for name, efpd, power, _ in STEPS
        )
        + "]\n"
    )
    return text


def files() -> dict[str, str]:
    return {
        "core-p9-assemblies.toml": assemblies(),
        "core-p9-cycle.toml": cycle(),
        "core-p9.toml": case(TARGET_K),
    }


def main() -> None:
    check = sys.argv[1:] == ["--check"]
    if sys.argv[1:] not in ([], ["--check"]):
        sys.exit("usage: python3 p9_timing_inputs.py [--check]")
    stale = []
    for name, text in files().items():
        path = HERE / name
        if check:
            if not path.exists() or path.read_text() != text:
                stale.append(name)
        else:
            path.write_text(text)
    if stale:
        sys.exit(f"stale: {', '.join(stale)}; run python3 p9_timing_inputs.py")


if __name__ == "__main__":
    main()
