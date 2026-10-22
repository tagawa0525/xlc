"""#0045 の試行の入力を書く（ADR-0185。UFG の格子の間引き）。

格子の間引きの許容差（`ufg-tolerance`）ごとに書く。`full` は `ufg-tolerance` を書かない（間引かない）。

- 2A-step: examples/depletion/p2/lattice-2A.toml の燃焼度の点を 0.1 MWd/kg だけにしたもの（輸送 3 回）
- 2A-60: 2A の燃焼の履歴（examples/depletion/p2/lattice-2A.msgpack。Slurm 434、lock のコミット 25ae2d40b）の 60 MWd/kg の
  組成（輸送する核種で 0 でないもの）を、燃焼の領域（対称操作の組）ごとの燃料の物質にした 1 状態。2 群の少数群定数を書く
- 1C-rings1、1C-rings4: examples/depletion/lattice-1C-rings{1,4}.toml の燃焼（5 MWd/kg まで）
- 2A-depletion: examples/depletion/p2/lattice-2A.toml の燃焼（40 点、60 MWd/kg まで）を間引いた格子で（間引かない履歴は
  examples/depletion/p2/lattice-2A.msgpack）
"""

import pathlib
import re

import msgpack

HERE = pathlib.Path(__file__).resolve().parent
DEPLETION = HERE.parent

VARIANTS = {"full": None, "t1e-3": 1.0e-3, "t1e-4": 1.0e-4}


def with_tolerance(text, tolerance):
    """`library = ...` の行のあとに `ufg-tolerance` を足す"""
    line = "" if tolerance is None else f"ufg-tolerance = {tolerance:g}\n"
    return re.sub(r"(?m)^(library = .*\n)", lambda m: m.group(1) + line, text, count=1)


def header(stem, source):
    return f"# #0045 の試行（make_inputs.py が書く。手で直さない）: {source}、{stem}\n"


def two_a_step():
    src = (DEPLETION / "p2" / "lattice-2A.toml").read_text()
    body = src[src.index("library = ") :]
    body = body.replace(
        'library = "depletion.msgpack"', 'library = "../p2/depletion.msgpack"'
    )
    body = re.sub(r"burnup = \[[^\]]*\]", "burnup = [0.1]", body, flags=re.DOTALL)
    for stem, tolerance in VARIANTS.items():
        name = f"2A-step-{stem}"
        text = with_tolerance(body, tolerance).replace(
            'output = "lattice-2A.msgpack"', f'output = "{name}.msgpack"'
        )
        (HERE / f"{name}.toml").write_text(
            header(name, "2A の燃焼の 1 ステップ（0.1 MWd/kg）") + text
        )


def octant_orbit(r, c, n):
    points = set()
    for a, b in ((r, c), (c, r)):
        for x in (a, n - 1 - a):
            for y in (b, n - 1 - b):
                points.add((x, y))
    return points


def two_a_burned():
    src = (DEPLETION / "p2" / "lattice-2A.toml").read_text()
    history = msgpack.unpackb((DEPLETION / "p2" / "lattice-2A.msgpack").read_bytes())
    names, keep = history["nuclides"], set(history["transportable"])
    last = history["points"][-1]
    assert last["burnup"] == 60.0, last["burnup"]
    rows = [r.split() for r in re.findall(r'^  "([^"]+)",$', src, flags=re.MULTILINE)]
    n = len(rows)
    # 燃焼の領域の代表のセル（余白の行と列を含む格子の番号）からピンの位置
    width = n + 2
    representative = {}
    for k, region in enumerate(history["regions"]):
        cell = region["cell"]
        representative[(cell // width - 1, cell % width - 1)] = k
    material_of = {}
    for (r, c), k in representative.items():
        for p in octant_orbit(r, c, n):
            assert rows[p[0]][p[1]] == "pin", p
            material_of[p] = k
    new_rows = []
    for r, row in enumerate(rows):
        new_rows.append(
            " ".join(
                f"p{material_of[(r, c)]}" if name == "pin" else name
                for c, name in enumerate(row)
            )
        )
    pin = src[src.index("[cell.pin]") : src.index("# 案内管")]
    others = src[src.index("# 案内管") : src.index("[depletion]")]
    materials = []
    for k, densities in enumerate(last["densities"]):
        lines = [
            f"[material.fuel-{k}]",
            "temperature = 565",
            "",
            f"[material.fuel-{k}.nuclides]",
        ]
        lines += [
            f"{names[i]} = {x:.6e}"
            for i, x in enumerate(densities)
            if names[i] in keep and x > 0.0
        ]
        materials.append("\n".join(lines) + "\n")
    rest = src[src.index("[material.zircaloy]") : src.index("# 燃料のピン")]
    cells = "\n".join(
        pin.replace("[cell.pin]", f"[cell.p{k}]").replace(
            'material = "fuel"', f'material = "fuel-{k}"'
        )
        for k in range(len(history["regions"]))
    )
    head = src[src.index("library = ") : src.index("rows = [")]
    head = head.replace(
        'library = "depletion.msgpack"', 'library = "../p2/depletion.msgpack"'
    )
    for stem, tolerance in VARIANTS.items():
        name = f"2A-60-{stem}"
        text = (
            header(name, "2A の 60 MWd/kg の組成の 1 状態")
            + with_tolerance(head, tolerance)
            + f'few-group = {{ output = "{name}.fg.msgpack", bounds = [1e7, 0.625, 1e-5] }}\n'
            + "rows = [\n"
            + "".join(f'  "{r}",\n' for r in new_rows)
            + "]\n\n"
            + "".join(m + "\n" for m in materials)
            + rest
            + cells
            + "\n"
            + others
        )
        (HERE / f"{name}.toml").write_text(text)


def one_c(rings):
    src = (DEPLETION / f"lattice-1C-{rings}.toml").read_text()
    body = src[src.index("library = ") :]
    body = body.replace(
        'library = "depletion.msgpack"', 'library = "../depletion.msgpack"'
    )
    body = body.replace('chain = "chain.msgpack"', 'chain = "../chain.msgpack"')
    for stem, tolerance in VARIANTS.items():
        name = f"1C-{rings}-{stem}"
        text = with_tolerance(body, tolerance).replace(
            f'output = "lattice-1C-{rings}.msgpack"', f'output = "{name}.msgpack"'
        )
        (HERE / f"{name}.toml").write_text(header(name, f"1C の燃焼（{rings}）") + text)


def two_a_depletion():
    src = (DEPLETION / "p2" / "lattice-2A.toml").read_text()
    body = src[src.index("library = ") :]
    body = body.replace(
        'library = "depletion.msgpack"', 'library = "../p2/depletion.msgpack"'
    )
    for stem, tolerance in VARIANTS.items():
        if tolerance is None:
            continue
        name = f"2A-depletion-{stem}"
        text = with_tolerance(body, tolerance).replace(
            'output = "lattice-2A.msgpack"', f'output = "{name}.msgpack"'
        )
        (HERE / f"{name}.toml").write_text(header(name, "2A の燃焼（40 点）") + text)


if __name__ == "__main__":
    two_a_depletion()
    two_a_step()
    two_a_burned()
    one_c("rings1")
    one_c("rings4")
