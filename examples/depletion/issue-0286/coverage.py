"""#0286 の、自由ガスから散乱の割合の小さい核種を外す（`ufg-free-gas.share`）と、外さない経路（`share = false`）の差を、
inputs.py の 9 ケースで測っていない条件で比べる入力を書く（ADR-0251 の「確認」）。

入力は coverage/ に `<ケース>-<変種>.toml` で書く。変種は `full`（`share = false`。どの核種も自由ガス）と `share-1e-4`。比べるのは
compare_coverage.py。

ケース:

- ピンセルの燃焼（`burn-uo2`、`burn-mox`）: 0 から 60 MWd/kg まで燃焼させ、燃焼点ごとの k∞ と数密度を比べる（変種ごとに別の履歴。
  閾値をまたぐ核種の段をここで見る）。形状は ../lattice-1C-rings4.toml（燃料 4 リングの燃焼の領域）、組成はウラン燃料が同じ入力、
  MOX が ../lattice-mox1-rings1.toml の MOX-1。ライブラリは ../p2/depletion.msgpack（分岐の燃料の温度 1200 K と減速材の 580、600 K を
  持つ）、燃料 900 K、被覆管と減速材 580 K、ホウ酸水 0.711901 g/cc、650 ppm（../issue-0180/branch/2A-60.toml の公称の状態）
- ピンセルの分岐（`branch-uo2`、`branch-mox`）: 燃焼の `full` の履歴の 0、20、40、60 MWd/kg から、公称、燃料の温度 1200 K、減速材の
  600 K（0.661165 g/cc）、その両方の 4 状態。ミクロ燃焼の 31 核種の微視的断面積も書く。どちらの変種も同じ組成（`full` の履歴）
- ピンセルの微分（`derivative-uo2`、`derivative-mox`）: 同じ履歴の 60 MWd/kg の公称の 1 状態で、Pu-240 と U-238 の数密度の微分を
  摂動の組成で自己遮蔽を解き直して求める（`derivative-shielding = "solved"`。一般化摂動の `shielded` と同じく、摂動の組成で外す
  核種を選び直す経路）
- 集合体の分岐（../issue-0180/fuels/inputs.py と同じ作り。燃焼の履歴の燃焼点から、履歴と同じ条件の分岐）: Gd の 2O の 10、40 MWd/kg
  （公称と減速材の 615 K（0.690 g/cc）の 2 状態。燃料の温度の分岐は、p2 の履歴が 600 K の燃料で H-1 を輸送し、ライブラリの H-1 に
  燃料の温度の格子点（900 K など）が無いので解けない）、Gd 24 本の 2P の 10、40 MWd/kg、Pyrex の 2F の
  40 MWd/kg、IFBA の 2M の 40 MWd/kg（公称の 1 状態）、本番の組成の 2A の 60 MWd/kg（../p9-product の履歴。../issue-0180/branch/2A-60.toml
  の公称と、燃料の温度 1200 K、減速材の 600 K（0.661165 g/cc）、その両方の 4 状態）。どれもミクロ燃焼の 31 核種を書く

使い方: python3 coverage.py（このディレクトリで）
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "coverage"
DEPLETION = HERE / ".."


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FUELS = load("fuels_inputs", DEPLETION / "issue-0180/fuels/inputs.py")

VARIANTS = {
    "full": "ufg-free-gas = { share = false }\n",
    "share-1e-4": "ufg-free-gas = { share = 1e-4 }\n",
}

# ADR-0035 の元素の輸送する核種 26 と、その輸送する先行核 5（ADR-0216。../issue-0180/branch/2A-60.toml と同じ）
MICRO = (
    "micro = [\n"
    '  "U-234", "U-235", "U-236", "U-238", "Np-237", "Np-239", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242",\n'
    '  "Am-241", "Am-243", "Xe-135", "Sm-147", "Sm-149", "Sm-150", "Sm-151", "Sm-152",\n'
    '  "Gd-152", "Gd-154", "Gd-155", "Gd-156", "Gd-157", "Gd-158", "Gd-160",\n'
    '  "Cm-242", "Cm-244", "Eu-153", "I-135", "Pm-148m",\n'
    "]\n"
)

# ピンセルの燃焼の点 (MWd/kg)。2 つの変種で同じ点なので、ステップの幅の誤差は差に入らない
BURNUP = [0.5, 2, 5, 10, 15, 20, 25, 30, 40, 50, 60]
# ピンセルの分岐を解く燃焼点
PIN_BRANCH_BURNUP = [0, 20, 40, 60]

# ピンセルの燃料の組成（燃料の 4 つの物質に同じものを置く）
PIN_FUELS = {
    "uo2": "nuclides = { U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02 }",
    "mox": (
        "nuclides = { U-235 = 5.1050e-05, U-238 = 2.0370e-02, Pu-238 = 4.6690e-05, Pu-239 = 1.4650e-03, Pu-240 = 5.6910e-04, "
        "Pu-241 = 2.7130e-04, Pu-242 = 1.4130e-04, Am-241 = 3.0280e-05, O-16 = 4.5880e-02 }"
    ),
}
# 集合体: （履歴の入力（examples/depletion から）、分岐を解く燃焼点、減速材の軸を 2 点にするか）
ASSEMBLIES = {
    "2O": ("p2/lattice-2O.toml", [10, 40], True),
    "2P": ("p2/lattice-2P.toml", [10, 40], False),
    "2F": ("p2/lattice-2F.toml", [40], False),
    "2M": ("p2/lattice-2M.toml", [40], False),
}


def replace_once(text: str, old: str, new: str, where: str) -> str:
    if text.count(old) != 1:
        sys.exit(f"{where}: expected {old!r} once")
    return text.replace(old, new)


def pin_lattice(fuel: str) -> str:
    """../lattice-1C-rings4.toml の `[depletion]` より前（先頭の注釈を除く）を、ライブラリと物質の温度とホウ酸水を替えたもの"""
    source = "lattice-1C-rings4.toml"
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") : text.find("[depletion]")]
    body = replace_once(
        body,
        'library = "depletion.msgpack"\n',
        'library = "../../p2/depletion.msgpack"\n',
        source,
    )
    old_fuel = PIN_FUELS["uo2"]
    if body.count(old_fuel) != 4:
        sys.exit(f"{source}: expected 4 fuel materials")
    body = body.replace(old_fuel, PIN_FUELS[fuel])
    body = replace_once(
        body,
        "[material.zircaloy]\ntemperature = 600\n",
        "[material.zircaloy]\ntemperature = 580\n",
        source,
    )
    body = replace_once(
        body,
        "[material.water]\ntemperature = 600\n"
        "nuclides = { H-1 = 4.41459e-02, O-16 = 2.20729e-02, B-10 = 9.52537e-06, B-11 = 3.83408e-05 }\n",
        "[material.water]\ntemperature = 580\nwater = { density = 0.711901, boron = 650, b10 = 0.199 }\n",
        source,
    )
    return body.rstrip() + "\n"


def listed(values) -> str:
    return "[" + ", ".join(f"{v:g}" for v in values) + "]"


def header(case: str, variant: str, what: str) -> str:
    return (
        "# このファイルは coverage.py が書く（手で直さない）。#0286 の自由ガスの核種の絞り込みの差（ADR-0251）: "
        f"ケース {case}、変種 {variant}。{what}\n"
    )


def with_variant(body: str, variant: str, where: str) -> str:
    return replace_once(
        body,
        'library = "../../p2/depletion.msgpack"\n',
        'library = "../../p2/depletion.msgpack"\n' + VARIANTS[variant],
        where,
    )


def burn(fuel: str, variant: str) -> str:
    case = f"burn-{fuel}"
    body = with_variant(pin_lattice(fuel), variant, case)
    return (
        header(case, variant, "ピンセルを 60 MWd/kg まで燃焼させる")
        + body
        + "\n[depletion]\n"
        'chain = "../../chain.msgpack"\n'
        "power = 38.41\n"
        f"burnup = {listed(BURNUP)}\n"
        'materials = ["fuel-0", "fuel-1", "fuel-2", "fuel-3"]\n'
        f'output = "{case}-{variant}.msgpack"\n'
        "few-group = { bounds = [1e7, 0.625, 1e-5] }\n"
    )


PIN_BRANCHES = (
    "[branches]\n"
    'boron = { ppm = [650], materials = ["water"] }\n'
    'moderator = { temperature = [580, 600], density = [0.711901, 0.661165], materials = ["water", "zircaloy"] }\n'
    'fuel = { temperature = [900, 1200], materials = ["fuel-0", "fuel-1", "fuel-2", "fuel-3"] }\n'
)
PIN_NOMINAL = (
    "[branches]\n"
    'boron = { ppm = [650], materials = ["water"] }\n'
    'moderator = { temperature = [580], density = [0.711901], materials = ["water", "zircaloy"] }\n'
    'fuel = { temperature = [900], materials = ["fuel-0", "fuel-1", "fuel-2", "fuel-3"] }\n'
)


def pin_branch(fuel: str, variant: str, derivative: bool) -> str:
    case = f"{'derivative' if derivative else 'branch'}-{fuel}"
    body = with_variant(pin_lattice(fuel), variant, case)
    body = replace_once(
        body,
        VARIANTS[variant],
        VARIANTS[variant]
        + f'few-group = {{ output = "table-{case}-{variant}.msgpack", bounds = [1e7, 0.625, 1e-5] }}\n',
        case,
    )
    history = (
        "[history]\n"
        f'file = "burn-{fuel}-full.msgpack"\n'
        'chain = "../../chain.msgpack"\n'
        f"burnup = {listed([60] if derivative else PIN_BRANCH_BURNUP)}\n" + MICRO
    )
    if derivative:
        history += 'derivative = ["Pu-240", "U-238"]\nderivative-shielding = "solved"\n'
    what = (
        "60 MWd/kg の公称の 1 状態で、Pu-240 と U-238 の微分の摂動の組成で自己遮蔽を解き直す"
        if derivative
        else "燃焼（full）の履歴からの 4 状態の分岐"
    )
    return (
        header(case, variant, what)
        + body
        + "\n"
        + (PIN_NOMINAL if derivative else PIN_BRANCHES)
        + "\n"
        + history
    )


def assembly(name: str, variant: str) -> str:
    source, burnup, axes = ASSEMBLIES[name]
    case = f"assembly-{name}"
    text = FUELS.branch(name, source, burnup, solved=False)
    text = text[text.find("library = ") :]
    text = replace_once(
        text, 'micro = ["Pu-240"]\nderivative = ["Pu-240"]\n', MICRO, source
    )
    text = replace_once(
        text,
        f"table-{name}-fixed.msgpack",
        f"table-{case}-{variant}.msgpack",
        source,
    )
    if axes:
        text = replace_once(
            text,
            'moderator = { temperature = [600], density = [0.743], materials = ["water"] }',
            'moderator = { temperature = [600, 615], density = [0.743, 0.690], materials = ["water"] }',
            source,
        )
    # FUELS.branch はパスを ../issue-0180/fuels から書き、coverage/ も examples/depletion から同じ深さなのでそのまま使う
    text = with_variant(text, variant, case)
    what = f"{source} の {listed(burnup)} MWd/kg からの{'2 状態' if axes else '公称の 1 状態'}の分岐"
    return header(case, variant, what) + text


def assembly_2a(variant: str) -> str:
    source = "issue-0180/branch/2A-60.toml"
    case = "assembly-2A"
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") :]
    body = replace_once(
        body,
        "table-2A-60.msgpack",
        f"table-{case}-{variant}.msgpack",
        source,
    )
    body = replace_once(
        body,
        "# 線形の補正の係数を σᵢ から数密度の微分にする核種（ADR-0221）\n"
        'derivative = ["Pu-240"]\n'
        "# 記録の値は摂動の組成で自己遮蔽を解き直して出した（ADR-0221。既定は状態の自己遮蔽のまま。ADR-0249）\n"
        'derivative-shielding = "solved"\n',
        "",
        source,
    )
    body = replace_once(
        body,
        "# 軸は公称の 1 点ずつ（1 状態）。密度は IF97 の 15.5132 MPa の値を小数 6 桁に丸めたもの（ADR-0239）\n"
        'moderator = { temperature = [580], density = [0.711901], materials = ["water", "zircaloy"] }\n',
        "# 公称と 1 点ずつ。密度は IF97 の 15.5132 MPa の値を小数 6 桁に丸めたもの（ADR-0239）\n"
        'moderator = { temperature = [580, 600], density = [0.711901, 0.661165], materials = ["water", "zircaloy"] }\n',
        source,
    )
    body = replace_once(
        body,
        "fuel = { temperature = [900], ",
        "fuel = { temperature = [900, 1200], ",
        source,
    )
    body = with_variant(body, variant, case)
    return header(case, variant, f"{source} の微分を除き、4 状態にしたもの") + body


def inputs() -> dict:
    """入力の名前（拡張子なし）と中身"""
    files = {}
    for variant in VARIANTS:
        for fuel in PIN_FUELS:
            files[f"burn-{fuel}-{variant}"] = burn(fuel, variant)
            files[f"branch-{fuel}-{variant}"] = pin_branch(fuel, variant, False)
            files[f"derivative-{fuel}-{variant}"] = pin_branch(fuel, variant, True)
        for name in ASSEMBLIES:
            files[f"assembly-{name}-{variant}"] = assembly(name, variant)
        files[f"assembly-2A-{variant}"] = assembly_2a(variant)
    return files


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for name, text in inputs().items():
        (OUT / f"{name}.toml").write_text(text)


if __name__ == "__main__":
    main()
