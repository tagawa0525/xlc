"""制御棒の部分挿入の例（ADR-0127）の炉心の入力を書く。

集合体1体（径方向は鏡面反射）の燃料の列（高さ 365.76 cm、VERA の燃料の長さ）を 16 ノード（22.86 cm）に分け、z の
両端は Marshak とする。制御棒は VERA の棒（仕様書の表 8: 先端から AIC 101.6 cm、その上に B4C 259.08 cm。ステップは
1.5875 cm）で、0 ステップの先端を燃料の下端に置く（VERA の全挿入の先端は燃料の下端の 5.08 cm 上だが、ここでは方式の
比較だけが目的なので、燃料の列に合わせる）。燃料の層のテーブルは rod_table.py が作る制御棒の軸（out、aic、b4c）を持つ。

ケースは棒の位置（ステップ）ごとに 3 つ:
- `flux-<ステップ>`: 部分挿入のノードを ADR-0127 の中性子束の重みと z の面の不連続因子で扱う
- `volume-<ステップ>`: 体積の重み（ADR-0127 で採らなかった案）
- `split-<ステップ>`: 軸方向のメッシュを先端と AIC/B4C の境でも分けた参照（部分挿入のノードがない。同じ定数と NEM）
と、細かいメッシュの幅の収束を見る `fine-<幅>-<ステップ>`。ステップは、ノード 8（182.88〜205.74 cm）を先端が
横切る 116〜128 の 2 ステップごとと、挿入の深さの違う 20、60、90、150、180、200。どの位置も、参照のメッシュに
1 cm より薄いノードを作らない（最も薄いのは 1.27 cm）。

ケースの名前の一覧は cases.yaml に書き（dvc.yaml の foreach が読む）、段階（ADR-0159）は dvc.yaml に書く（比較の段階の
deps にすべてのケースの出力を並べるため）。使い方: python3 rods_cases.py
"""

import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
HEIGHT = 365.76
NODES = 16
STEP = 1.5875
AIC = 101.6
B4C = 259.08
STEPS = [116, 118, 120, 122, 124, 126, 128, 20, 60, 90, 150, 180, 200]
# 細かいメッシュの幅の収束（ADR-0127 の既定の 0.5 cm と比べる）
FINE = {"2": 2.0, "1": 1.0, "0.25": 0.25, "0.125": 0.125}
FINE_STEPS = [116, 122, 180]

NOTE = "# このファイルは rods_cases.py が書く。"


def mesh(extra=()):
    bounds = [round(HEIGHT * i / NODES, 6) for i in range(NODES + 1)]
    bounds[-1] = HEIGHT
    for p in extra:
        if 0.0 < p < HEIGHT and p not in bounds:
            bounds.append(p)
    return sorted(bounds)


def numbers(values):
    return "[" + ", ".join(repr(v) for v in values) + "]"


def case(name, steps, axial, partial=None):
    lines = [
        f"{NOTE}{name}",
        'cycle = "../cycle.toml"',
        f"axial-mesh = {numbers(axial)}",
        "radial-nodes = 1",
        'boundary = { x = ["reflective", "reflective"], y = ["reflective", "reflective"], z = ["marshak", "marshak"] }',
        "",
        "[state.hzp]",
        "boron = 1300",
        "b10 = 0.199",
        "moderator = 600",
        "fuel = 600",
        f"banks = {{ D = {steps} }}",
    ]
    if partial:
        lines += ["", "[partial-rods]"] + partial
    (HERE / "cases" / f"{name}.toml").write_text("\n".join(lines) + "\n")
    return name


def main():
    (HERE / "plant.toml").write_text(
        f"""{NOTE}集合体1体の燃料の列と VERA の制御棒（仕様書の表 8）
pitch = 21.5
height = {HEIGHT!r}
map = ["*"]

[rods]
step = {STEP!r}
bottom = 0.0
steps = 230
top = {HEIGHT!r}
map = ["D"]

[rods.bank.D]
type = "rcca"

[rods.type.rcca]
segments = [{{ value = "aic", length = {AIC!r} }}, {{ value = "b4c", length = {B4C!r} }}]
"""
    )
    (HERE / "assemblies.toml").write_text(
        f"""{NOTE}燃料の層は制御棒の軸を持つテーブル（rod_table.py）
[design.fuel]
bounds = [0.0, {HEIGHT!r}]
layers = ["table-rods.msgpack"]

[assembly.F01]
design = "fuel"
"""
    )
    (HERE / "cycle.toml").write_text(
        f"""{NOTE}
plant = "plant.toml"
inventory = "assemblies.toml"
loading = ["F01"]
"""
    )
    (HERE / "cases").mkdir(exist_ok=True)
    names = []
    for s in STEPS:
        tip = 0.0 + s * STEP
        names.append(case(f"flux-{s}", s, mesh()))
        names.append(case(f"volume-{s}", s, mesh(), ['weighting = "volume"']))
        names.append(case(f"split-{s}", s, mesh([tip, tip + AIC])))
    for label, width in FINE.items():
        for s in FINE_STEPS:
            names.append(
                case(f"fine-{label}-{s}", s, mesh(), [f"fine-mesh = {width!r}"])
            )
    write_dvc(names)
    (HERE / "cases.yaml").write_text(
        "# このファイルは rods_cases.py が書く。dvc.yaml の foreach が読むケースの名前\ncases:\n"
        + "".join(f"  - {n}\n" for n in names)
    )


def code(stage):
    """Rust の段階の deps のうち、段階が呼ぶ実行ファイル（../../dvc.yaml の build が作る。ADR-0159）"""
    return [f"../../bin/{stage}"]


def listed(items, indent):
    return "".join(f"{' ' * indent}- {item}\n" for item in items)


def write_dvc(names):
    xs = (HERE / "xs.toml").read_text()
    endf = re.findall(r'^endf = "([^"]+)"', xs, flags=re.MULTILINE)
    thermal = re.findall(r'^thermal = "([^"]+)"', xs, flags=re.MULTILINE)
    endf_deps = [f"../../references/endf/ENDF-B-VIII.0_neutrons/{e}" for e in endf]
    endf_deps += [
        f"../../references/endf/ENDF-B-VIII.0_thermal_scatt/{t}"
        for t in sorted(set(thermal))
    ]
    text = f"""# このファイルは rods_cases.py が書く。制御棒の部分挿入の例（ADR-0127）の段階（ADR-0159）: ENDF → 多群ライブラリ（xs）→
# 制御棒の値ごとの格子計算（lattice-<値>）→ 制御棒の軸を持つ 1 つのテーブル（table。#0093 までのつなぎ）→ ケースごとの
# 炉心（core-<ケース>）→ 比較（compare）。実行するのは `dvc repro examples/core/dvc.yaml`
vars:
  - cases.yaml
stages:
  xs:
    cmd: ../../bin/xs xs.toml
    deps:
      - xs.toml
{listed(endf_deps + code("xs"), 6)}    outs:
      - rods.msgpack
  lattice:
    foreach: [out, aic, b4c]
    do:
      # 結果の k∞ は標準出力に出るので、ファイルに受ける
      cmd: ../../bin/lattice assembly-${{item}}.toml > lattice-${{item}}.out
      deps:
        - assembly-${{item}}.toml
        - rods.msgpack
{listed(code("lattice"), 8)}      outs:
        - fg-${{item}}.msgpack
        - lattice-${{item}}.out
  table:
    cmd: python3 rod_table.py table-rods.msgpack out=fg-out.msgpack aic=fg-aic.msgpack b4c=fg-b4c.msgpack
    deps:
      - rod_table.py
      - fg-out.msgpack
      - fg-aic.msgpack
      - fg-b4c.msgpack
    outs:
      - table-rods.msgpack
  core:
    foreach: ${{cases}}
    do:
      # DVC は出力の親のディレクトリを作らないので、作ってから実行する
      cmd: mkdir -p out && ../../bin/core cases/${{item}}.toml > out/${{item}}.out
      deps:
        - cases/${{item}}.toml
        - plant.toml
        - assemblies.toml
        - cycle.toml
        - table-rods.msgpack
{listed(code("core"), 8)}      outs:
        - out/${{item}}.out
  compare:
    cmd: python3 compare_rods.py > compare.out
    deps:
      - compare_rods.py
      - rods_cases.py
{listed([f"out/{n}.out" for n in names], 6)}    outs:
      - compare.out
"""
    (HERE / "dvc.yaml").write_text(text)


if __name__ == "__main__":
    main()
