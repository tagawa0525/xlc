"""#0286 の、UFG の自己遮蔽の項目ごとの時間と差の入力を書く。

燃料の種類（ケース）ごとに、自己遮蔽の条件（変種）を替えた 1 状態の入力を `<ケース>-<変種>.toml` に書く。変種 `default` は
ADR-0251 の前の既定（ADR-0209、ADR-0243、ADR-0246 の自由ガスで、どの核種も自由ガス。ufg-tolerance = 1e-4）で、測ったときの
経路で解くように `share = false` を書く。ほかは `default` に 1 つの項目を足したもの（`share` の変種は `share = false` の代わりに
閾値を書く）。どのケースの入力も、測ったときの ADR-0252 の前の既定の衝突確率の求積 16 点（`ufg-cp-points = 16`）を書き、求積の点の
変種はその行を置き換える。比べるのは compare.py。

ケース:

- ピンセル（VERA 1C の形状。../issue-0274 の #0133 のピンセルの組成と同じ。燃料は 5 リングで一様）: 新しいウラン燃料の 900 K
  （`pin-uo2`）と 600 K、1200 K（`pin-uo2-600`、`pin-uo2-1200`）、2A の 40 MWd/kg の組成（`pin-uo2-40`）、MOX-1（`pin-mox`。
  ../lattice-mox1-rings1.toml の燃焼度 0 の組成）
- 集合体（燃焼の履歴の燃焼点から、履歴と同じ条件の 1 状態の分岐。../issue-0180/fuels/inputs.py と同じ作り）: 2A の 60 MWd/kg
  （`2A-60`。../issue-0180/branch/2A-60.toml の微分を除いたもの）、Gd の 2O の 10 MWd/kg（`2O-10`）、Pyrex の 2F の 10 MWd/kg
  （`2F-10`）、制御棒を入れて燃焼させた 2A の 20 MWd/kg（`2A-rods-20`。棒は入ったまま）

静止した標的にする核種の選び方の変種（`SHARE_VARIANTS`。ADR-0251 の検討した案）は、2A-60 だけで時間を比べる。

衝突確率の近似（`CP_VARIANTS`）は、上のケースに 2 つを足して比べる（変種の名前は `cp-` で始まる）: Gd の 2O の 5 MWd/kg（`2O-5`。
Gd のピンのリングの組成の差が 10 MWd/kg より大きい途中）と、燃焼した MOX のピンセル（`pin-mox-40`。mox-history.toml の
40 MWd/kg からの 1 状態の分岐。燃料はリングごとに組成が違う）。

燃焼した MOX の集合体（`mox-assembly-60`。../mox-assembly の履歴の 60 MWd/kg からの 1 状態の分岐）は、`MOX_ASSEMBLY_TIMING` の変種で
時間だけを比べる（差は coverage.py の `assembly-mox` の 4 状態で比べる）。

使い方: python3 inputs.py（このディレクトリで）
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEPLETION = HERE / ".."


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FUELS = load("fuels_inputs", DEPLETION / "issue-0180/fuels/inputs.py")

# 変種: 既定の入力の `library = ` の行の後に足す行
# （`static` は自由ガスを使わないので `share` を書けない）
# `share-*` の変種のうち、#0286 の最初の形（レサジーの積分の比。ADR-0251 の案 A）で測ったもの。既定の測り方（ADR-0251 の
# 格子点ごとの比の最大）と違うので、測り方を書く
LETHARGY_1E_4 = (
    'ufg-free-gas = { share = { threshold = 1e-4, measure = "lethargy" } }\n'
)
VARIANTS = {
    "default": "ufg-free-gas = { share = false }\n",
    "static": 'ufg-scattering = "static"\n',
    "heavy-static": "ufg-free-gas = { heavy = false, share = false }\n",
    "light-400kt": "ufg-free-gas = { upper = 400, share = false }\n",
    "sweeps-2": "ufg-free-gas = { sweeps = 2, share = false }\n",
    "sweeps-1": "ufg-free-gas = { sweeps = 1, share = false }\n",
    "tolerance-3e-4": "ufg-free-gas = { share = false }\n",
    "tolerance-1e-3": "ufg-free-gas = { share = false }\n",
    "table-8": "track-table = 8\nufg-free-gas = { share = false }\n",
    "share-1e-4": LETHARGY_1E_4,
    "share-1e-3": 'ufg-free-gas = { share = { threshold = 1e-3, measure = "lethargy" } }\n',
    "share-1e-2": 'ufg-free-gas = { share = { threshold = 1e-2, measure = "lethargy" } }\n',
}
# 静止した標的にする核種の選び方の変種（ADR-0251 の検討した案。`share-1e-4` からの時間の比を 2A-60 で測る）: 格子点ごとの割合の
# 最大で測る形と、アクチニドを割合によらず自由ガスに残す形
SHARE_VARIANTS = {
    "share-point": 'ufg-free-gas = { share = { threshold = 1e-4, measure = "point" } }\n',
    "share-actinides": 'ufg-free-gas = { share = { threshold = 1e-4, measure = "lethargy", keep-actinides = true } }\n',
}
SHARE_CASES = ("2A-60",)
# 衝突確率の近似の変種（#0286 の衝突確率の節。`cp-default` は ADR-0251 の前の既定で、測ったときの経路で解くように `share = false`
# を書く。ほかは `cp-default` に足す行で、`share` の変種は `share = false` の代わりに閾値を書く）: 弦の求積の点（`ufg-cp-points`）、
# 飛程の表の細かさ（`track-table`）、展開の上限（`track-table-deviation`）と、その組み合わせ、自由ガスの核種の絞り込みとの組み合わせ
CP_COMBINED = "ufg-cp-points = 8\ntrack-table = 8\ntrack-table-deviation = 0.03\n"
CP_FULL = "ufg-free-gas = { share = false }\n"
CP_VARIANTS = {
    "cp-default": CP_FULL,
    "cp-points-8": CP_FULL + "ufg-cp-points = 8\n",
    "cp-points-6": CP_FULL + "ufg-cp-points = 6\n",
    "cp-points-4": CP_FULL + "ufg-cp-points = 4\n",
    "cp-table-8": CP_FULL + "track-table = 8\n",
    "cp-deviation-0.03": CP_FULL + "track-table-deviation = 0.03\n",
    "cp-table-8-deviation-0.03": CP_FULL
    + "track-table = 8\ntrack-table-deviation = 0.03\n",
    "cp-combined": CP_FULL + CP_COMBINED,
    "cp-share-1e-4": LETHARGY_1E_4,
    "cp-share-1e-4-combined": LETHARGY_1E_4 + CP_COMBINED,
    "cp-share-1e-4-points-6": LETHARGY_1E_4 + "ufg-cp-points = 6\n",
}
# 衝突確率の近似だけで比べるケース: （履歴の入力（examples/depletion から）、分岐を解く燃焼点 (MWd/kg)）
CP_ASSEMBLIES = {"2O-5": ("p2/lattice-2O.toml", 5)}
# 燃焼した MOX のピンセルの分岐（mox-history.toml の燃焼点）
MOX_BURNUP = 40
# 燃焼した MOX の集合体（../mox-assembly の履歴）の 1 状態で時間を比べる変種: どの核種も自由ガスの 16 点（基準）、既定の候補の
# 選び方（ADR-0251）、求積 6 点（ADR-0252）、その両方
MOX_ASSEMBLY_CASE = "mox-assembly-60"
MOX_ASSEMBLY_BURNUP = 60
MOX_ASSEMBLY_VARIANTS = {
    "cp-share-point-points-6": SHARE_VARIANTS["share-point"] + "ufg-cp-points = 6\n",
}
MOX_ASSEMBLY_TIMING = (
    "default",
    "share-point",
    "cp-points-6",
    "cp-share-point-points-6",
)

# 間引きの許容差の変種（既定の 1e-4 を置き換える）
TOLERANCES = {"tolerance-3e-4": "3e-4", "tolerance-1e-3": "1e-3"}

PIN_SOURCES = {
    "pin-uo2": ("issue-0274/pin-1C-0-energy-1000.toml", None),
    "pin-uo2-600": ("issue-0274/pin-1C-0-energy-1000.toml", 600),
    "pin-uo2-1200": ("issue-0274/pin-1C-0-energy-1000.toml", 1200),
    "pin-uo2-40": ("issue-0274/pin-2A-40-energy-1000.toml", None),
}
# ../issue-0274 のピンセルの、自己遮蔽の条件の行（既定にするため除く）
PIN_SETTINGS = (
    'ufg-scattering = "free-gas"\n',
    "ufg-free-gas = { upper = 400, cells = 2, sweeps = 3, heavy = { energy = 1000 }, sab = true }\n",
    'ufg-transfer = "infinite-dilution"\n',
)
# 集合体: （履歴の入力（examples/depletion から）、分岐を解く燃焼点 (MWd/kg)）
ASSEMBLIES = {
    "2O-10": ("p2/lattice-2O.toml", 10),
    "2F-10": ("p2/lattice-2F.toml", 10),
    "2A-rods-20": ("p2/lattice-2A-rods.toml", 20),
}


def replace_once(text: str, old: str, new: str, where: str) -> str:
    if text.count(old) != 1:
        sys.exit(f"{where}: expected {old!r} once")
    return text.replace(old, new)


# 測ったときの既定（ADR-0252 の前）の求積の点。既定が変わっても、変種を測ったときと同じ 16 点で解く
EXACT_POINTS = "ufg-cp-points = 16\n"


def with_exact_points(text: str, where: str) -> str:
    """`library = ` の行の後に [`EXACT_POINTS`] を足す（履歴の入力が既に書いていればそのまま）"""
    if EXACT_POINTS.rstrip("\n") in text.splitlines():
        return text
    library = text.find("\nlibrary = ")
    library_end = text.find("\n", library + 1) + 1
    if library < 0 or library_end == 0:
        sys.exit(f"{where}: no library line")
    return text[:library_end] + EXACT_POINTS + text[library_end:]


def pin(name: str) -> str:
    source, temperature = PIN_SOURCES[name]
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") :]
    for line in PIN_SETTINGS:
        body = replace_once(body, line, "", source)
    body = replace_once(
        body,
        'library = "../p2/depletion.msgpack"\n',
        'library = "../p2/depletion.msgpack"\nufg-tolerance = 1e-4\n',
        source,
    )
    if temperature is not None:
        body = replace_once(
            body,
            "[material.fuel]\ntemperature = 900\n",
            f"[material.fuel]\ntemperature = {temperature}\n",
            source,
        )
    return (
        f"# {source} の組成と形状"
        + (f"（燃料 {temperature} K）\n" if temperature else "\n")
        + body
    )


def mox() -> str:
    source = "lattice-mox1-rings1.toml"
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") : text.find("[depletion]")]
    body = replace_once(
        body,
        'library = "depletion.msgpack"\n',
        'library = "../depletion.msgpack"\nufg-tolerance = 1e-4\n',
        source,
    )
    return f"# {source} の燃焼度 0 の組成と形状\n" + body.rstrip() + "\n"


def assembly_2a60() -> str:
    source = "issue-0180/branch/2A-60.toml"
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") :].replace('"../../', '"../')
    # micro の行の前の注釈（追う核種の説明）から除く
    start = body.find("# ADR-0035 の元素の輸送する核種")
    end = body.find('derivative-shielding = "solved"\n')
    if start < 0 or end < 0:
        sys.exit(f"{source}: no micro or derivative-shielding line")
    body = body[:start] + body[end + len('derivative-shielding = "solved"\n') :]
    body = replace_once(body, "table-2A-60.msgpack", "{table}", source)
    return f"# {source} の微分を除いたもの\n" + body


def assembly(name: str) -> str:
    source, burnup = {**ASSEMBLIES, **CP_ASSEMBLIES}[name]
    text = FUELS.branch(name, source, [burnup], solved=False)
    text = text[text.find("library = ") :]
    text = replace_once(
        text, 'micro = ["Pu-240"]\nderivative = ["Pu-240"]\n', "", source
    )
    text = replace_once(text, f"table-{name}-fixed.msgpack", "{table}", source)
    # FUELS.branch はパスを ../issue-0180/fuels から書く
    text = text.replace('"../../', '"../')
    return f"# {source} の {burnup} MWd/kg からの 1 状態の分岐\n" + text


def mox_burned() -> str:
    """mox-history.toml の `[depletion]` より前（UFG の格子の間引きを含む）に、少数群定数の出力と、履歴と同じ条件の 1 状態の分岐と `[history]` を足したもの"""
    source = "issue-0286/mox-history.toml"
    text = (HERE / "mox-history.toml").read_text()
    body = text[text.find("library = ") : text.find("[depletion]")]
    # 履歴の自由ガスの経路は変種の行が決めるので、履歴の入力の行を除く
    body = replace_once(
        body,
        "# 履歴（mox-history.msgpack）は ADR-0251 の前の既定（どの核種も自由ガス）で作った\n"
        "ufg-free-gas = { share = false }\n",
        "",
        source,
    )
    # 求積の点は、どのケースとも同じく [`with_exact_points`] の行で書く
    body = replace_once(
        body,
        "# 履歴は ADR-0252 の前の既定（衝突確率の求積 16 点）で作った\n" + EXACT_POINTS,
        "",
        source,
    )
    body = replace_once(
        body,
        'library = "../p2/depletion.msgpack"\n',
        'library = "../p2/depletion.msgpack"\n'
        'few-group = { output = "{pin-table}", bounds = [1e7, 0.625, 1e-5] }\n',
        source,
    )
    fuels = ", ".join(f'"fuel-{k}"' for k in range(5))
    return (
        f"# {source} の {MOX_BURNUP} MWd/kg からの 1 状態の分岐\n"
        + body.rstrip()
        + "\n\n[branches]\n"
        + 'boron = { ppm = [650], materials = ["water"] }\n'
        + 'moderator = { temperature = [580], density = [0.711901], materials = ["water", "zircaloy"] }\n'
        + f"fuel = {{ temperature = [900], materials = [{fuels}] }}\n\n"
        + '[history]\nfile = "mox-history.msgpack"\nchain = "../chain.msgpack"\n'
        + f"burnup = [{MOX_BURNUP}]\n"
    )


def mox_assembly() -> str:
    """../mox-assembly/lattice-mox-assembly.toml の `[depletion]` より前に、テーブルの出力と、履歴と同じ条件の 1 状態の分岐と `[history]` を
    足したもの"""
    source = "mox-assembly/lattice-mox-assembly.toml"
    text = (DEPLETION / source).read_text()
    body = text[text.find("library = ") : text.find("[depletion]")]
    # 履歴の自由ガスの経路は変種の行が決める。`ufg-hold-memory` は結果を変えず時間を変える（ADR-0213）ので除く。履歴と同じ
    # `ufg-cp-points = 16` は基準の側に残す（既定が変わっても基準を履歴と同じ 16 点に保つ）
    for line in ("ufg-free-gas = { share = false }\n", "ufg-hold-memory = 32\n"):
        body = replace_once(body, line, "", source)
    body = replace_once(
        body,
        'library = "../p2/depletion.msgpack"\n',
        'library = "../p2/depletion.msgpack"\n'
        'few-group = { output = "{table}", bounds = [1e7, 0.625, 1e-5] }\n',
        source,
    )
    return (
        f"# {source} の {MOX_ASSEMBLY_BURNUP} MWd/kg からの 1 状態の分岐\n"
        + body.rstrip()
        + "\n\n[branches]\n"
        + 'boron = { ppm = [650], materials = ["water"] }\n'
        + 'moderator = { temperature = [580], density = [0.711901], materials = ["water", "zircaloy"] }\n'
        + 'fuel = { temperature = [900], materials = ["mox-high", "mox-mid", "mox-low"] }\n\n'
        + '[history]\nfile = "../mox-assembly/lattice-mox-assembly.msgpack"\nchain = "../chain.msgpack"\n'
        + f"burnup = [{MOX_ASSEMBLY_BURNUP}]\n"
    )


def variant(base: str, case: str, name: str) -> str:
    """ケースの既定の入力 `base` に変種 `name` の行を足し、少数群定数（ピンセルは few-group、集合体はテーブル）の出力の名前を付ける"""
    text = base
    line = {**VARIANTS, **CP_VARIANTS, **SHARE_VARIANTS, **MOX_ASSEMBLY_VARIANTS}[name]
    if "ufg-cp-points = " in line:
        # 求積の点の変種は、どのケースも書く [`EXACT_POINTS`] の行を変種の行で置き換える
        text = replace_once(text, EXACT_POINTS, "", case)
    if name in TOLERANCES:
        text = replace_once(
            text,
            "ufg-tolerance = 1e-4\n",
            f"ufg-tolerance = {TOLERANCES[name]}\n",
            case,
        )
    head, rest = text.split("\n", 1)
    library_end = rest.find("\n", rest.find("library = ")) + 1
    rest = rest[:library_end] + line + rest[library_end:]
    output = f"fg-{case}-{name}.msgpack"
    if "{pin-table}" in rest:
        # 燃焼した MOX のピンセル（[`mox_burned`]）は 1 本なのでピン出力を書かない
        rest = replace_once(rest, '"{pin-table}"', f'"{output}"', case)
    elif "{table}" in rest:
        # 集合体はピン出力の形状関数もテーブルに書く（比べるピン出力）
        rest = replace_once(
            rest,
            '"{table}", bounds = [1e7, 0.625, 1e-5] }',
            f'"{output}", bounds = [1e7, 0.625, 1e-5], pins = true }}',
            case,
        )
    else:
        rest = replace_once(
            rest,
            "ufg-tolerance = ",
            f'few-group = {{ output = "{output}", bounds = [1e7, 0.625, 1e-5] }}\nufg-tolerance = ',
            case,
        )
    return (
        "# このファイルは inputs.py が書く（手で直さない）。#0286 の UFG の自己遮蔽の項目ごとの時間と差: "
        f"ケース {case}、変種 {name}\n{head}\n{rest}"
    )


def cases() -> list[str]:
    """ケースの名前（[`bases`] と [`cp_bases`] の順）"""
    return [*PIN_SOURCES, "pin-mox", "2A-60", *ASSEMBLIES, *CP_ASSEMBLIES, "pin-mox-40"]


def bases() -> dict:
    cases = {name: pin(name) for name in PIN_SOURCES}
    cases["pin-mox"] = mox()
    cases["2A-60"] = assembly_2a60()
    cases.update({name: assembly(name) for name in ASSEMBLIES})
    return {name: with_exact_points(text, name) for name, text in cases.items()}


def cp_bases() -> dict:
    """衝突確率の近似だけで比べるケースの既定の入力"""
    cases = {name: assembly(name) for name in CP_ASSEMBLIES}
    cases["pin-mox-40"] = mox_burned()
    return {name: with_exact_points(text, name) for name, text in cases.items()}


def main() -> None:
    old = bases()
    for case, base in old.items():
        for name in VARIANTS:
            (HERE / f"{case}-{name}.toml").write_text(variant(base, case, name))
    for case in SHARE_CASES:
        for name in SHARE_VARIANTS:
            (HERE / f"{case}-{name}.toml").write_text(variant(old[case], case, name))
    for case, base in {**old, **cp_bases()}.items():
        for name in CP_VARIANTS:
            (HERE / f"{case}-{name}.toml").write_text(variant(base, case, name))
    base = with_exact_points(mox_assembly(), MOX_ASSEMBLY_CASE)
    for name in MOX_ASSEMBLY_TIMING:
        (HERE / f"{MOX_ASSEMBLY_CASE}-{name}.toml").write_text(
            variant(base, MOX_ASSEMBLY_CASE, name)
        )


if __name__ == "__main__":
    main()
