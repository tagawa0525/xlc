"""#0129 の試行の入力を書く。

燃料のリングの組成が違うピンセル（燃焼なしの 1 状態）:

- 1C-5: examples/depletion/lattice-1C-rings4.msgpack の最後の点（5 MWd/kg）のリングごとの組成（輸送する核種で 0 でないもの）
- gd-mid: 2O の Gd のピン（8 リング）で、Gd-155/157 が外側から燃えた組成を手で作ったもの（燃焼の計算の結果ではない）
- gd-fresh: 2O の Gd のピンの新しい組成（リングは同じ組成）

燃焼（ADR-0173）:

- 1C-rings4: examples/depletion/lattice-1C-rings4.toml をそのまま（燃焼度の点も同じ）
- 2O-step: examples/depletion/p2/lattice-2O.toml の燃焼度の点を 0.1 MWd/kg だけにし、2 群の少数群定数を足したもの

それぞれ track-table = false（exact）と既定（table）の 2 つを書く。
"""

import pathlib
import re

import msgpack
import tomllib

HERE = pathlib.Path(__file__).resolve().parent
DEPLETION = HERE.parent


def materials_toml(materials):
    out = []
    for name, (temperature, nuclides) in materials.items():
        out.append(
            f"[material.{name}]\ntemperature = {temperature}\n\n[material.{name}.nuclides]"
        )
        out.extend(f"{k} = {v:.6e}" for k, v in nuclides.items())
        out.append("")
    return "\n".join(out)


def write(stem, header, materials, cells, rows):
    for kind, line in (("exact", "track-table = false\n"), ("table", "")):
        text = (
            header.replace("@STEM@", f"{stem}-{kind}")
            + line
            + f"rows = {rows}\n\n"
            + materials_toml(materials)
            + "\n"
            + cells
        )
        (HERE / f"{stem}-{kind}.toml").write_text(text)


def one_c():
    src = (DEPLETION / "lattice-1C-rings4.toml").read_text()
    spec = tomllib.loads(src)
    history = msgpack.unpackb((DEPLETION / "lattice-1C-rings4.msgpack").read_bytes())
    names, keep = history["nuclides"], set(history["transportable"])
    last = history["points"][-1]
    materials = {}
    for region, densities in zip(history["regions"], last["densities"]):
        nuclides = {
            names[i]: x for i, x in enumerate(densities) if names[i] in keep and x > 0.0
        }
        materials[region["material"]] = (
            spec["material"][region["material"]]["temperature"],
            nuclides,
        )
    for name in ("zircaloy", "water"):
        materials[name] = (
            spec["material"][name]["temperature"],
            spec["material"][name]["nuclides"],
        )
    cells = src[src.index("[cell.pin]") : src.index("[depletion]")]
    header = (
        "# #0129 の試行: 1C の 4 リングの 5 MWd/kg の組成（make_inputs.py が書く）\n"
        'library = "../depletion.msgpack"\n'
        'medium = "@STEM@.medium.msgpack"\n'
        "pitch = 1.26\nazimuthal = 64\nspacing = 0.02\n"
    )
    write("1C-5", header, materials, cells, '["pin"]')


def gadolinium():
    src = (DEPLETION / "p2" / "lattice-2O.toml").read_text()
    spec = tomllib.loads(src)
    cells = src[src.index("# Gd 燃料のピン") : src.index("# 案内管")]
    header = (
        "# #0129 の試行: 2O の Gd のピンのピンセル（make_inputs.py が書く）\n"
        'library = "../p2/depletion.msgpack"\n'
        'medium = "@STEM@.medium.msgpack"\n'
        "pitch = 1.26\nazimuthal = 30\nspacing = 0.03\n"
    )
    base = spec["material"]
    shared = {
        n: (base[n]["temperature"], base[n]["nuclides"]) for n in ("zircaloy", "water")
    }
    fresh = dict(shared)
    for m in range(8):
        fresh[f"gd{m}"] = (
            base[f"gd{m}"]["temperature"],
            dict(base[f"gd{m}"]["nuclides"]),
        )
    write("gd-fresh", header, fresh, cells, '["gd"]')
    # 外側から燃えた Gd（中心から外へ）: Gd-155 の残る割合。Gd-157 はその 2 乗。Pu-239 は外側ほど多い
    left = [0.95, 0.93, 0.9, 0.85, 0.75, 0.55, 0.3, 0.08]
    pu = [1.0, 1.0, 1.0, 1.0, 1.1, 1.3, 1.6, 2.2]
    mid = dict(shared)
    for m in range(8):
        t, n = base[f"gd{m}"]["temperature"], dict(base[f"gd{m}"]["nuclides"])
        g5, g7 = n["Gd-155"], n["Gd-157"]
        n["Gd-155"], n["Gd-157"] = g5 * left[m], g7 * left[m] ** 2
        n["Gd-156"] += g5 - n["Gd-155"]
        n["Gd-158"] += g7 - n["Gd-157"]
        n["Pu-239"] = 2.0e-5 * pu[m]
        n["U-235"] *= 0.97 - 0.01 * m
        mid[f"gd{m}"] = (t, n)
    write("gd-mid", header, mid, cells, '["gd"]')


def without_track_table(text):
    """入力の track-table の行と、その前の説明の行を除く"""
    lines = text.splitlines(keepends=True)
    out = []
    for line in lines:
        if line.startswith("track-table"):
            while (
                out
                and out[-1].startswith("#")
                and ("#0129" in out[-1] or "飛程の表" in out[-1])
            ):
                out.pop()
            continue
        out.append(line)
    return "".join(out)


def depletion_pair(stem, text):
    table = without_track_table(text)
    head = table.index("library =")
    for kind, line in (("exact", "track-table = false\n"), ("table", "")):
        body = table[:head] + line + table[head:]
        body = re.sub(r'output = "[^"]*"', f'output = "{stem}-{kind}.msgpack"', body)
        (HERE / f"{stem}-{kind}.toml").write_text(
            "# #0129 の試行（make_inputs.py が書く）\n" + body
        )


def depletion():
    one = (DEPLETION / "lattice-1C-rings4.toml").read_text()
    one = one.replace(
        'library = "depletion.msgpack"',
        'library = "../depletion.msgpack"',
    )
    one = one.replace('chain = "chain.msgpack"', 'chain = "../chain.msgpack"')
    depletion_pair("1C-rings4", one)
    two = (DEPLETION / "p2" / "lattice-2O.toml").read_text()
    two = two.replace(
        'library = "depletion.msgpack"',
        'library = "../p2/depletion.msgpack"',
    )
    # p2 の chain = "../chain.msgpack" は、このディレクトリからも同じチェーンを指す
    two = re.sub(r"burnup = \[[^\]]*\]", "burnup = [0.1]", two, flags=re.DOTALL)
    two = two.replace(
        'materials = ["fuel"',
        'few-group = { bounds = [1e7, 0.625, 1e-5] }\nmaterials = ["fuel"',
    )
    depletion_pair("2O-step", two)


one_c()
gadolinium()
depletion()
