"""#0180 の、燃料の種類ごとの微分の摂動の自己遮蔽を固定する誤差（ADR-0249）の入力を書く。

燃焼の履歴（../../p2 の集合体の燃焼）の燃焼点から、履歴と同じ条件の 1 状態の
分岐を解き、Pu-240 の数密度の微分（ADR-0221）を持つテーブルを書く。種類ごとに 2 つ:

- <種類>-solved.toml: 摂動の組成で自己遮蔽を解き直す（`derivative-shielding = "solved"`）
- <種類>-fixed.toml: 状態の自己遮蔽のまま、燃焼の領域の物質だけを摂動の組成で作り直す（既定）

格子、物質、離散化は履歴の入力（`[depletion]` より前）と同じで、燃料の温度の軸 1 点（入力の燃料の温度）の分岐にする。比べるのは
compare.py。

使い方: python3 inputs.py（このディレクトリで）
"""

import sys
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
DEPLETION = HERE / "../.."

# 種類: （履歴の入力（examples/depletion から）、分岐を解く燃焼点 (MWd/kg)）。燃焼点は、Gd が残る 10 MWd/kg と燃焼の進んだ
# 40 MWd/kg。IFBA の 2M は自己遮蔽が 2A の約 3 倍なので 20 MWd/kg の 1 点。MOX のピンセルの履歴（../../lattice-mox1-rings1）は、
# その履歴のライブラリ（../../depletion.msgpack）が今の形式より古く読めず（`delayed` がない）、p2 のライブラリでは履歴の輸送する
# 核種と 900 K のデータの核種が違って分岐が止まるので、入れていない
CASES = {
    "2F": ("p2/lattice-2F.toml", [10, 40]),
    "2M": ("p2/lattice-2M.toml", [20]),
    "2O": ("p2/lattice-2O.toml", [10, 40]),
    "2P": ("p2/lattice-2P.toml", [10, 40]),
}
# 分岐の軸は 3 つとも要る（ホウ素の軸はホウ酸水で書いた物質）ので、履歴の入力の減速材の数密度をホウ酸水で書き直す（密度、
# ホウ素、B-10 の割合。H-1 と O-16 の数密度と相対で約 1e-6 違う。../../branches/branch-2A-micro.toml と同じ）
WATERS = {
    "nuclides = { O-16 = 2.48112e-02, H-1 = 4.96224e-02, B-10 = 1.07070e-05, B-11 = 4.30971e-05 }": (
        0.743,
        1300,
    ),
}


def listed(values) -> str:
    return "[" + ", ".join(f"{v:g}" for v in values) + "]"


def written(value) -> str:
    """TOML の値を、入力に書いたときの形で（整数の温度は整数のまま）"""
    return f"{value:g}" if isinstance(value, float) else str(value)


def branch(name: str, source: str, burnup: list, solved: bool) -> str:
    """履歴の入力 `source` の `[depletion]` より前（先頭の注釈を除く）に、テーブルの出力と、燃料の温度 1 点の分岐と `[history]` を
    足したもの"""
    path = DEPLETION / source
    text = path.read_text()
    start, end = text.find("library = "), text.find("[depletion]")
    if start < 0 or end < 0:
        sys.exit(f"{path} has no library line or [depletion]")
    data = tomllib.loads(text)
    output, depleted, library = (
        data["depletion"]["output"],
        data["depletion"]["materials"],
        data["library"],
    )
    materials = ", ".join(f'"{m}"' for m in depleted)
    temperatures = {m: v["temperature"] for m, v in data["material"].items()}
    body = text[start:end].replace(
        f'library = "{library}"',
        f'library = "{(path.parent / library).resolve().relative_to(HERE.resolve(), walk_up=True)}"\n'
        f'few-group = {{ output = "table-{name}-{"solved" if solved else "fixed"}.msgpack", '
        "bounds = [1e7, 0.625, 1e-5] }",
        1,
    )
    if "few-group = " not in body:
        sys.exit(f"{path} does not have the library line to add the table output")
    water = [(n, w) for n, w in WATERS.items() if n in body]
    if len(water) != 1:
        sys.exit(f"{path} does not have one water of WATERS")
    ((nuclides, (density, boron)),) = water
    body = body.replace(
        nuclides, f"water = {{ density = {density}, boron = {boron}, b10 = 0.199 }}"
    )
    fuel, moderator = written(temperatures["fuel"]), written(temperatures["water"])
    # 燃料の温度の軸は燃焼させる物質のすべてにかかるので、どれも履歴の入力で燃料と同じ温度であること（違えば、履歴と同じ条件の
    # 状態にならない）
    for material in depleted:
        if temperatures[material] != temperatures["fuel"]:
            sys.exit(f"{path}: {material} is not at the fuel temperature {fuel}")
    history = (path.parent / output).resolve().relative_to(HERE.resolve(), walk_up=True)
    shielding = 'derivative-shielding = "solved"\n' if solved else ""
    return (
        f"# このファイルは inputs.py が書く（手で直さない）。#0180 の燃料の種類ごとの微分の摂動の自己遮蔽（ADR-0249）: {source}\n"
        f"{body.rstrip()}\n\n"
        f'[branches]\nboron = {{ ppm = [{boron}], materials = ["water"] }}\n'
        f'moderator = {{ temperature = [{moderator}], density = [{density}], materials = ["water"] }}\n'
        f"fuel = {{ temperature = [{fuel}], materials = [{materials}] }}\n\n"
        f'[history]\nfile = "{history}"\nchain = "../../chain.msgpack"\nburnup = {listed(burnup)}\n'
        f'micro = ["Pu-240"]\nderivative = ["Pu-240"]\n{shielding}'
    )


def main() -> None:
    for name, (source, burnup) in CASES.items():
        for solved in (True, False):
            file = HERE / f"{name}-{'solved' if solved else 'fixed'}.toml"
            file.write_text(branch(name, source, burnup, solved))


if __name__ == "__main__":
    main()
