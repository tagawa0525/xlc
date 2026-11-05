"""VERA Problem 6（3D の HFP の単一集合体）を xlc で解く入力を書く。

Problem 6 の形状と組成は Problem 3A と同じ（仕様書 CASL-U-2012-0131-004 の 78 ページ）なので、格子、平板、炉心のメッシュは
examples/vera/p3/p3_inputs.py の 3A の関数で書き、物質の状態を、ADR-0230 の格子点の分岐のテーブル（ADR-0097）に置き換える。
置き換える文字列は p3_inputs.py が書く文字列そのもので、見つからなければ止まる（p3_inputs.py の形が変わったことを黙って
見過ごさないため）。

このファイルのある場所に次を書く。

- `lattice-6-<層>.toml`: 燃料の高さの 2D の単一集合体の格子（3A の fuel、zgrid、igrid）の分岐のテーブル。支持格子のない層は、
  平板の燃料にする媒質（`medium-6-fuel.msgpack`）も書く（ADR-0128）
- `axial-6.toml`: 軸方向の反射体の平板（ADR-0128）の分岐のテーブル。軸は格子と同じ（媒質は格子の同じ状態を使う）
- `core-6.toml` と、共通の `core-6-cycle.toml`、`core-6-plant.toml`、`core-6-assemblies.toml`（ADR-0105）: 表 P6-1 の HFP の
  熱水力のフィードバックの解（ADR-0229）
- `core-6-low-power.toml`: 同じ熱水力で出力だけを 1e-6 倍にした解（ノードの温度がほぼ入口の温度になる）
- `core-6-isothermal.toml`: 熱水力なしで、状態の減速材と燃料を入口の温度にした解（上の解と k が合うことを確かめる）
- `dvc.yaml`: ステージ

使い方: python3 p6_inputs.py
"""

import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "p3"))
import p3_inputs as p3

# 表 P6-1（仕様書 78 ページ）
INLET = 565  # K
PRESSURE_PSIA = 2250
BORON = 1300  # ppm
POWER = 17.67  # MW（3411 MW / 193 体）
FLOW_MLBS_HR = 0.6823  # バイパスの 9% を除き 193 体で割った流量
# B-10 の原子分率: Problem 3A の冷却材（表 P3-3）の B-10 / (B-10 + B-11) = 1.07070e-05 / 5.38041e-05（ADR-0106）
B10 = 0.199
# 1 lb = 0.45359237 kg（国際ポンド）、1 psi = 6.894757293168e-3 MPa（crates/core/src/thermal/water.rs の PSI と同じ）
FLOW = FLOW_MLBS_HR * 1e6 * 0.45359237 / 3600.0  # kg/s
PRESSURE = PRESSURE_PSIA * 6.894757293168e-3  # MPa
# 燃料で発生する熱の割合（表 18、19 ページ）
FUEL_HEAT_FRACTION = 0.974
# ギャップのコンダクタンス (W/(m² K)): NEA の PWR MOX/UO2 炉心過渡のベンチマーク（NEA/NSC/DOC(2006)20）の仕様の式 (8) の
# hgap。ADR-0229 は既定を持たない必須の入力とし（案 H1）、物理から導いた値ではない（案 H3）。機構のモデルは #0228 に残る
GAP_CONDUCTANCE = 1.0e4
# 燃料の実効温度: 体積平均（ADR-0229 の案 E1。Kelly ら（Nuclear Engineering and Technology 49 (2017) 1326）の MC21/CTF、
# MPACT/CTF の比較が、CTF の体積平均の温度を MC21 に渡すのと同じ）
EFFECTIVE = "volume-average"
TOLERANCE = 1.0e-5  # ADR-0229 の既定と同じ値を、記録が引けるように書く
# 燃料棒（表 1）と、案内管と計装管の外半径（表 2）、集合体のピッチ（表 P3-1）
PELLET, CLAD_IN, CLAD_OUT = 0.4096, 0.418, 0.475
RODS, GUIDE_TUBES, GUIDE_OUT, INSTRUMENT_OUT = 264, 24, 0.602, 0.605
ASSEMBLY_PITCH = 21.50
FUEL_DENSITY = 10.257  # g/cc（表 P6-1）
FLOW_AREA = ASSEMBLY_PITCH**2 - math.pi * (
    RODS * CLAD_OUT**2 + GUIDE_TUBES * GUIDE_OUT**2 + INSTRUMENT_OUT**2
)
WETTED_PERIMETER = (
    2 * math.pi * (RODS * CLAD_OUT + GUIDE_TUBES * GUIDE_OUT + INSTRUMENT_OUT)
)

# 分岐のテーブルの格子点（ADR-0230）。減速材の密度は IF97 の 15.5132 MPa の値（crates/core/src/thermal/water.rs の
# `region_one`）を小数 6 桁に丸めたもの。core は格子点の密度を IF97 と 5e-6 g/cc で確かめる（ADR-0229）
MODERATOR = [
    (565, 0.742713),
    (570, 0.732936),
    (580, 0.711901),
    (590, 0.688335),
    (600, 0.661165),
    (610, 0.628290),
]
FUEL = [565, 600, 750, 900, 1050, 1200]
NOMINAL_T, NOMINAL_DENSITY = MODERATOR[0]
if NOMINAL_T != INLET or FUEL[0] != INLET:
    raise SystemExit("the nominal state is the inlet temperature on both axes")

P3_HEADER = "# このファイルは p3_inputs.py が書く。手で直さない。\n"
HEADER = "# このファイルは p6_inputs.py が書く（形状は p3_inputs.py の 3A）。手で直さない。\n"
P3_PROBLEM = f"VERA Problem 3A（{p3.PROBLEMS['3A']['description']}"
# 状態は表 P6-1 と ADR-0230 の格子点。3A の格子の見出しは、この後に「。表 P3-1、P3-3）」（形状と組成の出典）が続く
PROBLEM = "VERA Problem 6（Problem 3A の集合体の HFP の分岐のテーブル。状態は表 P6-1 と ADR-0230 の格子点（ADR-0097）"


def replaced(text: str, old: str, new: str) -> str:
    """`old`（1 回だけ現れること）を `new` に置き換える"""
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"p3_inputs.py wrote {n} copies of {old[:60]!r}")
    return text.replace(old, new)


def renamed(text: str) -> str:
    """3A の出力のパスを 6 のものにする。残れば止まる"""
    for old, new in [
        ("lattice-3A/", "lattice-6/"),
        ("lattice-3A-", "lattice-6-"),
        ("medium-3A-", "medium-6-"),
        ("axial-3A/", "axial-6/"),
    ]:
        text = text.replace(old, new)
    if "-3A" in text or "3A/" in text:
        raise SystemExit("a Problem 3A path is left")
    return text


def materials(names: list[str]) -> str:
    """物質の TOML（`names` の順）。公称の状態（入口の温度）。冷却材はホウ酸水（ADR-0106）"""
    t = NOMINAL_T
    blocks = {
        "fuel": f"[material.fuel]\ntemperature = {FUEL[0]}\nnuclides = {{ {p3.PROBLEMS['3A']['fuel']} }}",
        "water": f"[material.water]\ntemperature = {t}\n"
        f"water = {{ density = {NOMINAL_DENSITY}, boron = {BORON}, b10 = {B10} }}",
        "zircaloy": f"[material.zircaloy]\ntemperature = {t}\n\n[material.zircaloy.nuclides]\n{p3.ZIRCALOY}",
        "ss304": f"[material.ss304]\ntemperature = {t}\n\n[material.ss304.nuclides]\n{p3.SS304}",
        "inconel": f"[material.inconel]\ntemperature = {t}\n\n[material.inconel.nuclides]\n{p3.INCONEL}",
    }
    return "\n\n".join(blocks[n] for n in names)


def quoted(names: list[str]) -> str:
    return ", ".join(f'"{n}"' for n in names)


def branches(borated: list[str], moderated: list[str], fuels: list[str]) -> str:
    temperatures = ", ".join(str(t) for t, _ in MODERATOR)
    densities = ", ".join(f"{d:.6f}" for _, d in MODERATOR)
    return (
        "[branches]\n"
        f"boron = {{ ppm = [{BORON}], materials = [{quoted(borated)}] }}\n"
        f"moderator = {{ temperature = [{temperatures}], density = [{densities}], materials = [{quoted(moderated)}] }}\n"
        f"fuel = {{ temperature = [{', '.join(str(t) for t in FUEL)}], materials = [{quoted(fuels)}] }}"
    )


def lattice(name: str, insert: str, grid: str | None, what: str) -> str:
    """3A の 2D の格子の、物質を分岐のテーブルにしたもの"""
    text = p3.lattice("3A", name, insert, grid, what)
    names = ["fuel", "water", "zircaloy"] + (["inconel"] if grid == "inconel" else [])
    moderated = names[1:]
    text = replaced(
        text,
        p3.materials("3A", names),
        materials(names) + "\n\n" + branches(["water"], moderated, ["fuel"]),
    )
    text = replaced(text, P3_HEADER, HEADER)
    text = replaced(text, P3_PROBLEM, PROBLEM)
    return renamed(text)


def slab() -> str:
    """3A の軸方向の反射体の平板の、物質を分岐のテーブルにしたもの。体積で混ぜた物質はどれも水を含むので、ホウ素と
    減速材の軸に挙げる（lattice::branches）"""
    text = p3.slab("3A")
    names = ["water", "zircaloy", "ss304", "inconel"]
    mixes = list(p3.FIXED_MIXTURES)
    for _, _, cell in p3.physical_layers("3A"):
        if cell.startswith("s") and cell not in mixes:
            mixes.append(cell)
    text = replaced(
        text,
        p3.materials("3A", names),
        materials(names)
        + "\n\n"
        + branches(["water", *mixes], [*names, *mixes], ["medium-fuel"]),
    )
    text = replaced(text, P3_HEADER, HEADER)
    text = replaced(text, P3_PROBLEM, PROBLEM)
    return renamed(text)


def thermal(power: float, what: str) -> str:
    """`what` は出力の出どころ（入力のコメントに書く）"""
    return f"""[thermal-hydraulics]
# 出力 (MW) は{what}。表 P6-1 の流量 0.6823 Mlbs/hr (kg/s)、圧力 2250 psia (MPa)。燃料の発熱の割合は表 18
power = {power:g}
flow = {FLOW:.4f}
pressure = {PRESSURE:.6f}
fuel-heat-fraction = {FUEL_HEAT_FRACTION}
# 体積平均（ADR-0229 の案 E1。Kelly ら（2017）の MC21/CTF、MPACT/CTF の比較と同じ）
effective-temperature = "{EFFECTIVE}"
tolerance = {TOLERANCE:.1e}
# 燃料棒の寸法（表 1）、ペレットの密度（表 P6-1）。ギャップのコンダクタンスは NEA/NSC/DOC(2006)20 の式 (8) の値（ADR-0229 の
# 案 H1、H3。物理から導いた値ではない）
rod = {{ pellet-radius = {PELLET}, clad-inner-radius = {CLAD_IN}, clad-outer-radius = {CLAD_OUT}, pellet-density = {FUEL_DENSITY}, gap-conductance = {GAP_CONDUCTANCE:.1e}, rods = {RODS} }}
# 集合体の断面 21.50 cm 角から燃料棒 264 本、案内管 24 本、計装管 1 本の外径の円を除いた面積 (cm²) と、その円周の和 (cm)
channel = {{ flow-area = {FLOW_AREA:.4f}, wetted-perimeter = {WETTED_PERIMETER:.4f} }}

[state.hfp]
boron = {BORON}
b10 = {B10}
inlet = {INLET}
"""


CASES = {
    "core-6": (
        "表 P6-1 の HFP の熱水力のフィードバックの解（ADR-0229）",
        thermal(POWER, "表 P6-1 の定格"),
    ),
    "core-6-low-power": (
        "core-6 の出力だけを 1e-6 倍にした熱水力の解（ノードの温度はほぼ入口の温度。core-6-isothermal と k を比べる）",
        thermal(POWER * 1e-6, "表 P6-1 の定格の 1e-6 倍"),
    ),
    "core-6-isothermal": (
        "熱水力なしで、減速材と燃料を入口の温度にした解（core-6-low-power と k を比べる）",
        f"[state.inlet]\nboron = {BORON}\nb10 = {B10}\nmoderator = {INLET}\nfuel = {INLET}\n",
    ),
}


def core_files() -> dict[str, str]:
    """炉心の入力（ADR-0105）。燃料の層は 2D の格子のテーブル、軸方向の反射体は平板のテーブル（f* の方法。ADR-0128）。
    メッシュと層は 3A の f* の炉心と同じ"""
    ns = p3.nodes("3A")
    top = p3.decimal(ns[-1][2] + p3.OFFSET)
    mesh = ", ".join(p3.decimal(lo + p3.OFFSET) for _, lo, _, _ in ns) + f", {top}"
    layers = [
        f"axial-6/{name}.msgpack"
        if reflector
        else f"lattice-6/{p3.fuel_layer('3A', (lo + hi) / 2)}.msgpack"
        for name, lo, hi, reflector in ns
    ]
    files = {}
    for stem, (what, body) in CASES.items():
        files[f"{stem}.toml"] = (
            f"# このファイルは p6_inputs.py が書く。手で直さない。VERA Problem 6 の炉心: {what}。\n"
            "# 燃料は 2D の格子のテーブル、軸方向の反射体は平板のテーブル（f* の方法。ADR-0128）。径方向は鏡面反射（集合体 1 体）、\n"
            "# z の両端は真空。z は下の水の外端から（仕様書の高さ + 55 cm）\n"
            'cycle = "core-6-cycle.toml"\n'
            f"axial-mesh = [{mesh}]\n"
            "radial-nodes = 1\n"
            'boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], '
            'z = ["vacuum", "vacuum"] }\n\n' + body
        )
    files["core-6-cycle.toml"] = (
        "# このファイルは p6_inputs.py が書く\n"
        'plant = "core-6-plant.toml"\n'
        'inventory = "core-6-assemblies.toml"\n'
        'loading = ["6"]\n'
    )
    files["core-6-plant.toml"] = (
        "# このファイルは p6_inputs.py が書く。集合体 1 体（ピッチ (cm) は集合体のピッチ 21.50 cm）\n"
        f'pitch = 21.5\nheight = {top}\nmap = ["*"]\n'
    )
    files["core-6-assemblies.toml"] = (
        "# このファイルは p6_inputs.py が書く\n"
        f"[design.6]\nbounds = [{mesh}]\nlayers = [\n"
        + "".join(f'  "{layer}",\n' for layer in layers)
        + ']\n\n[assembly.6]\ndesign = "6"\n'
    )
    return files


def python_ast(dep: str) -> str:
    """Python のファイルの deps は、リポジトリの直下の dvc.yaml の python-ast が書く正規形にする（ADR-0210）"""
    source = os.path.relpath(os.path.normpath(HERE / dep), ROOT)
    return os.path.relpath(ROOT / "python-ast" / f"{source}.ast", HERE)


def stage(name: str, cmd: str, deps: list[str], outs: list[str]) -> str:
    """`git:` の出力は、キャッシュに置かず git で追跡する（`cache: false`。ADR-0176）"""
    return (
        f"  {name}:\n    cmd: {cmd}\n    deps:\n"
        + "".join(
            f"      - {python_ast(d) if d.endswith('.py') else d}\n" for d in deps
        )
        + "    outs:\n"
        + "".join(
            f"      - {o.removeprefix('git:')}:\n          cache: false\n"
            if o.startswith("git:")
            else f"      - {o}\n"
            for o in outs
        )
    )


def pipeline() -> str:
    """ステージ: ライブラリ（../dvc.yaml の xs）→ 2D の格子のテーブル → 平板のテーブル → 炉心 → 比較"""
    bin_ = "../../../bin"
    stages = []
    tables = []
    for name, _, grid, _ in p3.LATTICES["3A"]:
        stem = f"lattice-6-{name}"
        outs = [f"lattice-6/{name}.msgpack", f"{stem}.out"]
        if grid is None:
            outs.append(f"medium-6-{name}.msgpack")
        tables.append(outs[0])
        stages.append(
            stage(
                stem,
                f"mkdir -p lattice-6 && {bin_}/lattice {stem}.toml > {stem}.out",
                [f"{stem}.toml", "../vera.msgpack", f"{bin_}/lattice"],
                outs,
            )
        )
    # 平板の少数群定数はディレクトリに書く（DVC はステージの前に出力を消すので、作ってから実行する）
    stages.append(
        stage(
            "axial-6",
            f"mkdir -p axial-6 && RUST_LOG=debug {bin_}/lattice axial-6.toml > axial-6.out 2> axial-6.log",
            [
                "axial-6.toml",
                "../vera.msgpack",
                "medium-6-fuel.msgpack",
                f"{bin_}/lattice",
            ],
            ["axial-6", "axial-6.out", "axial-6.log"],
        )
    )
    shared = ["core-6-cycle.toml", "core-6-plant.toml", "core-6-assemblies.toml"]
    for stem in CASES:
        stages.append(
            stage(
                stem,
                f"RUST_LOG=debug {bin_}/core {stem}.toml > {stem}.out 2> {stem}.log",
                [f"{stem}.toml", *shared, "axial-6", *tables, f"{bin_}/core"],
                [f"{stem}.out", f"{stem}.log"],
            )
        )
    stages.append(
        stage(
            "compare-6",
            "python3 ../compare_p6.py --json compare-6.json lattice-6-fuel.toml core-6.toml core-6-low-power.toml "
            "core-6-isothermal.toml core-6.out core-6-low-power.out core-6-isothermal.out",
            [
                "../compare_p6.py",
                "lattice-6-fuel.toml",
                "core-6.toml",
                "core-6-low-power.toml",
                "core-6-isothermal.toml",
                "core-6.out",
                "core-6-low-power.out",
                "core-6-isothermal.out",
            ],
            ["git:compare-6.json"],
        )
    )
    return (
        "# このファイルは p6_inputs.py が書く。手で直さない。\n"
        "# VERA Problem 6 を xlc で解くステージ: 多群ライブラリ（../dvc.yaml の xs）→ 燃料の高さの 2D の格子の分岐のテーブル（層ごと。\n"
        "# 平板の燃料の媒質も）→ 軸方向の反射体の平板の分岐のテーブル（axial-6）→ 炉心（熱水力のフィードバックの core-6 と、\n"
        "# 確認の core-6-low-power、core-6-isothermal）→ 比較（compare-6.json）。実行するのは `dvc repro`（ADR-0159）\n"
        "stages:\n" + "".join(stages)
    )


def main() -> None:
    if sys.argv[1:]:
        raise SystemExit("usage: python3 p6_inputs.py")
    (HERE / "dvc.yaml").write_text(pipeline())
    for name, insert, grid, what in p3.LATTICES["3A"]:
        (HERE / f"lattice-6-{name}.toml").write_text(lattice(name, insert, grid, what))
    (HERE / "axial-6.toml").write_text(slab())
    for name, text in core_files().items():
        (HERE / name).write_text(text)


if __name__ == "__main__":
    main()
