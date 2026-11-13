"""#0286 の、UFG の自己遮蔽の項目ごとの時間と差の入力を書く。

燃料の種類（ケース）ごとに、自己遮蔽の条件（変種）を替えた 1 状態の入力を `<ケース>-<変種>.toml` に書く。変種 `default` は
既定（ADR-0209、ADR-0243、ADR-0246 の自由ガス、ufg-tolerance = 1e-4）で、ほかは既定に 1 つの項目を足したもの。比べるのは
compare.py。

ケース:

- ピンセル（VERA 1C の形状。../issue-0274 の #0133 のピンセルの組成と同じ。燃料は 5 リングで一様）: 新しいウラン燃料の 900 K
  （`pin-uo2`）と 600 K、1200 K（`pin-uo2-600`、`pin-uo2-1200`）、2A の 40 MWd/kg の組成（`pin-uo2-40`）、MOX-1（`pin-mox`。
  ../lattice-mox1-rings1.toml の燃焼度 0 の組成）
- 集合体（燃焼の履歴の燃焼点から、履歴と同じ条件の 1 状態の分岐。../issue-0180/fuels/inputs.py と同じ作り）: 2A の 60 MWd/kg
  （`2A-60`。../issue-0180/branch/2A-60.toml の微分を除いたもの）、Gd の 2O の 10 MWd/kg（`2O-10`）、Pyrex の 2F の 10 MWd/kg
  （`2F-10`）、制御棒を入れて燃焼させた 2A の 20 MWd/kg（`2A-rods-20`。棒は入ったまま）

衝突確率の近似（`CP_VARIANTS`）は、上のケースに 2 つを足して比べる（変種の名前は `cp-` で始まる）: Gd の 2O の 5 MWd/kg（`2O-5`。
Gd のピンのリングの組成の差が 10 MWd/kg より大きい途中）と、燃焼した MOX のピンセル（`pin-mox-40`。mox-history.toml の
40 MWd/kg からの 1 状態の分岐。燃料はリングごとに組成が違う）。

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

# 変種: 既定の入力の `library = ` の行の後に足す行（default は何も足さない）
VARIANTS = {
    "default": "",
    "static": 'ufg-scattering = "static"\n',
    "heavy-static": "ufg-free-gas = { heavy = false }\n",
    "light-400kt": "ufg-free-gas = { upper = 400 }\n",
    "sweeps-2": "ufg-free-gas = { sweeps = 2 }\n",
    "sweeps-1": "ufg-free-gas = { sweeps = 1 }\n",
    "tolerance-3e-4": "",
    "tolerance-1e-3": "",
    "table-8": "track-table = 8\n",
    "share-1e-4": "ufg-free-gas = { share = 1e-4 }\n",
    "share-1e-3": "ufg-free-gas = { share = 1e-3 }\n",
    "share-1e-2": "ufg-free-gas = { share = 1e-2 }\n",
}
# 衝突確率の近似の変種（#0286 の衝突確率の節。`default` は既定のまま。ほかは既定に足す行）: 弦の求積の点（`ufg-cp-points`）、
# 飛程の表の細かさ（`track-table`）、展開の上限（`track-table-deviation`）と、その組み合わせ、自由ガスの核種の絞り込みとの組み合わせ
CP_COMBINED = "ufg-cp-points = 8\ntrack-table = 8\ntrack-table-deviation = 0.03\n"
CP_VARIANTS = {
    "cp-default": "",
    "cp-points-8": "ufg-cp-points = 8\n",
    "cp-points-6": "ufg-cp-points = 6\n",
    "cp-points-4": "ufg-cp-points = 4\n",
    "cp-table-8": "track-table = 8\n",
    "cp-deviation-0.03": "track-table-deviation = 0.03\n",
    "cp-table-8-deviation-0.03": "track-table = 8\ntrack-table-deviation = 0.03\n",
    "cp-combined": CP_COMBINED,
    "cp-share-1e-4": "ufg-free-gas = { share = 1e-4 }\n",
    "cp-share-1e-4-combined": "ufg-free-gas = { share = 1e-4 }\n" + CP_COMBINED,
}
# 衝突確率の近似だけで比べるケース: （履歴の入力（examples/depletion から）、分岐を解く燃焼点 (MWd/kg)）
CP_ASSEMBLIES = {"2O-5": ("p2/lattice-2O.toml", 5)}
# 燃焼した MOX のピンセルの分岐（mox-history.toml の燃焼点）
MOX_BURNUP = 40

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
    """mox-history.toml の `[depletion]` より前に、少数群定数の出力と、履歴と同じ条件の 1 状態の分岐と `[history]` を足したもの"""
    source = "issue-0286/mox-history.toml"
    text = (HERE / "mox-history.toml").read_text()
    body = text[text.find("library = ") : text.find("[depletion]")]
    body = replace_once(
        body,
        'library = "../p2/depletion.msgpack"\n',
        'library = "../p2/depletion.msgpack"\n'
        'few-group = { output = "{pin-table}", bounds = [1e7, 0.625, 1e-5] }\n'
        "ufg-tolerance = 1e-4\n",
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


def variant(base: str, case: str, name: str) -> str:
    """ケースの既定の入力 `base` に変種 `name` の行を足し、少数群定数（ピンセルは few-group、集合体はテーブル）の出力の名前を付ける"""
    text = base
    if name in TOLERANCES:
        text = replace_once(
            text,
            "ufg-tolerance = 1e-4\n",
            f"ufg-tolerance = {TOLERANCES[name]}\n",
            case,
        )
    head, rest = text.split("\n", 1)
    library_end = rest.find("\n", rest.find("library = ")) + 1
    rest = rest[:library_end] + {**VARIANTS, **CP_VARIANTS}[name] + rest[library_end:]
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
    return cases


def cp_bases() -> dict:
    """衝突確率の近似だけで比べるケースの既定の入力"""
    cases = {name: assembly(name) for name in CP_ASSEMBLIES}
    cases["pin-mox-40"] = mox_burned()
    return cases


def main() -> None:
    old = bases()
    for case, base in old.items():
        for name in VARIANTS:
            (HERE / f"{case}-{name}.toml").write_text(variant(base, case, name))
    for case, base in {**old, **cp_bases()}.items():
        for name in CP_VARIANTS:
            (HERE / f"{case}-{name}.toml").write_text(variant(base, case, name))


if __name__ == "__main__":
    main()
