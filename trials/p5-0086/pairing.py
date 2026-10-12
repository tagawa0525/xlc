"""#0086 の試行: 5-2D の外周の、代表（3.10%、Pyrex なし）と違う種類の燃料に 1 面で面する反射体の位置に、代表の I 型の格子
の反射体のノードの定数でなく、その燃料の I 型の格子の反射体のノードの定数を置く（反射体の定数の燃料の種類への依存）。

`<問題>-paired.toml` ほか（p5 の `core-<問題>.toml` の炉心の形だけを替えたもの）を書く。

使い方（リポジトリの直下で、references の devShell の中）: python3 trials/p5-0086/pairing.py
"""

import re
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
P5 = HERE.parent.parent / "examples/vera/p5"
REL = "../../examples/vera/p5"
# 向き（行の差、列の差）。上の行が +y。回転 0 の I 型の格子の燃料は +x に反射体
STEPS = [(0, 1), (-1, 0), (0, -1), (1, 0)]


def main() -> None:
    for problem in ["5A", "5B", "5C"]:
        cycle = tomllib.loads((P5 / f"core-{problem}-cycle.toml").read_text())
        plant_text = (P5 / cycle["plant"]).read_text()
        plant = tomllib.loads(plant_text)
        inv = tomllib.loads((P5 / cycle["inventory"]).read_text())
        grid = [row.split() for row in plant["map"]]
        load = [row.split() for row in cycle["loading"]]
        kinds = {}
        for i, row in enumerate(load):
            for j, cell in enumerate(row):
                if cell == ".":
                    continue
                aid, _, turn = cell.partition("@")
                design = inv["assembly"][aid]["design"]
                m = re.fullmatch(r"(.+)-I", design)
                if not m or m.group(1) == "31":
                    continue
                di, dj = STEPS[int(turn or 0) // 90]
                site = grid[i + di][j + dj]
                name, _, rot = site.partition("@")
                if name != "I":
                    raise SystemExit(f"{aid} faces {site}")
                kind = "I" + m.group(1).replace("-", "")
                kinds[kind] = m.group(1)
                grid[i + di][j + dj] = kind + (f"@{rot}" if rot else "")
        width = max(len(c) for row in grid for c in row) + 1
        rows = "".join(
            '  "' + "".join(c.ljust(width) for c in row).rstrip() + '",\n'
            for row in grid
        )
        blocks = "".join(
            f'\n[reflector.{k}]\nbounds = {v["bounds"]}\nlayers = ["{REL}/{v["layers"][0]}"]\n'
            for k, v in plant["reflector"].items()
        )
        blocks += "".join(
            f'\n[reflector.{k}]\nbounds = [0.0, 21.5]\nlayers = ["{REL}/lattice/reflector-I-{d}.msgpack"]\n'
            for k, d in sorted(kinds.items())
        )
        name = f"{problem}-paired"
        (HERE / f"{name}-plant.toml").write_text(
            f"pitch = {plant['pitch']}\nheight = {plant['height']}\nmap = [\n{rows}]\n{blocks}"
        )
        loading = "".join(
            '  "' + "".join(c.ljust(width) for c in row).rstrip() + '",\n'
            for row in load
        )
        (HERE / f"{name}-cycle.toml").write_text(
            f'plant = "{name}-plant.toml"\ninventory = "{REL}/{cycle["inventory"]}"\nloading = [\n{loading}]\n'
        )
        case = (P5 / f"core-{problem}.toml").read_text()
        case = case.replace(f"core-{problem}-cycle.toml", f"{name}-cycle.toml")
        case = case.replace(f"core-{problem}-pins.msgpack", f"{name}-pins.msgpack")
        (HERE / f"{name}.toml").write_text(case)


if __name__ == "__main__":
    main()
