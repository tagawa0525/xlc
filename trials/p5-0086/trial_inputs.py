"""#0086 の切り分けの試行の入力を書く: 反射体に面する燃料のノードの少数群定数の項（断面積、拡散係数、ADF）を、反射体の
格子の燃料のノードのものに替えた表と、それを置いた core の入力。

- `tables/<格子>-<項>.msgpack`: 単一集合体 31（3.10%、Pyrex なし）の定数のうち、<項> だけを反射体の格子（I 型か角型）の
  燃料のノードのものに替えた表。項は xs（total、absorption、nu_fission、fission、chi、scatter）、D（diffusion）、adf と、
  その組（xsD、xsDadf）
- `I-<項>.toml`: I 型の格子と同じ体系（燃料 1 体と I 型の反射体）。燃料のノードに上の表を置く。比べる相手は
  p5/reflector-I.out の k
- `<炉心>-<名前>.toml`: p5 の炉心の入力（`core-<炉心>.toml`）の、反射体に面する燃料の位置（1 面なら I 型、2 面なら
  角型の格子の燃料のノード）の定数を替えたもの。元の設計が 31 の位置だけを替える（Pyrex のある 31-p8、31-p12 の位置は
  その組の格子がない）。名前は、両方を替える IC、I 型だけの I、角型だけの C、項だけを替える IC-<項>

使い方（リポジトリの直下で、references の devShell の中）: python3 trials/p5-0086/trial_inputs.py
"""

from pathlib import Path

import msgpack
import tomllib

HERE = Path(__file__).resolve().parent
P5 = HERE.parent.parent / "examples/vera/p5"
LATTICE = P5 / "lattice"
REL = "../../examples/vera/p5"
XS = ["total", "absorption", "nu_fission", "fission", "chi", "scatter"]
TERMS = {
    "xs": XS,
    "D": ["diffusion"],
    "xsD": XS + ["diffusion"],
    "adf": ["adf"],
    "xsDadf": XS + ["diffusion", "adf"],
}
# 向き（行の差、列の差）。上の行が +y
DIRS = {"+x": (0, 1), "+y": (-1, 0), "-x": (0, -1), "-y": (1, 0)}
TURN_OF_I = {"+x": 0, "+y": 90, "-x": 180, "-y": 270}
# 角型の格子の燃料は +x と +y に反射体がある。その 2 つを回した組 → 回転
TURN_OF_CORNER = {
    ("+x", "+y"): 0,
    ("+y", "-x"): 90,
    ("-x", "-y"): 180,
    ("-y", "+x"): 270,
}


def mix(donor: str, keys: list[str], out: Path) -> None:
    base = msgpack.unpackb((LATTICE / "31.msgpack").read_bytes())
    other = msgpack.unpackb((LATTICE / f"reflector-{donor}-fuel.msgpack").read_bytes())
    if len(base["states"]) != 1 or len(other["states"]) != 1:
        raise SystemExit("only tables with one state")
    for key in keys:
        if key not in base["states"][0]:
            raise SystemExit(f"unknown key {key}")
        base["states"][0][key] = other["states"][0][key]
    out.write_bytes(msgpack.packb(base))


def write_lattice_core(term: str) -> None:
    """I 型の格子と同じ体系の core の入力"""
    (HERE / f"I-{term}-assemblies.toml").write_text(
        f'[design.f]\nbounds = [0.0, 21.5]\nlayers = ["tables/I-{term}.msgpack"]\n\n[assembly.f]\ndesign = "f"\n'
    )
    (HERE / f"I-{term}-cycle.toml").write_text(
        f'plant = "{REL}/core-lattice-I-plant.toml"\ninventory = "I-{term}-assemblies.toml"\nloading = ["f ."]\n'
    )
    (HERE / f"I-{term}.toml").write_text(
        f'cycle = "I-{term}-cycle.toml"\naxial-mesh = [0.0, 21.5]\nradial-nodes = 1\n'
        'boundary = { x = ["reflective", "vacuum"], y = ["reflective", "reflective"], z = ["reflective", "reflective"] }\n'
    )


def write_periphery(problem: str, name: str, one: str | None, two: str | None) -> None:
    """炉心の入力の、反射体に 1 面で面する設計 31 の燃料の位置に `one`、2 面で面する位置に `two` の表（反射体が +x、
    または +x と +y の向きの定数）を、反射体の向きに回して置く。None なら替えない"""
    cycle = tomllib.loads((P5 / f"core-{problem}-cycle.toml").read_text())
    plant = tomllib.loads((P5 / cycle["plant"]).read_text())
    inv = tomllib.loads((P5 / cycle["inventory"]).read_text())
    grid = [row.split() for row in plant["map"]]
    load = [row.split() for row in cycle["loading"]]
    assemblies = dict(inv["assembly"])
    added = {}
    for i, row in enumerate(grid):
        for j, mark in enumerate(row):
            if mark != "*":
                continue
            faces = [d for d, (di, dj) in DIRS.items() if grid[i + di][j + dj] != "*"]
            if any(grid[i + DIRS[d][0]][j + DIRS[d][1]] == "-" for d in faces):
                raise SystemExit(f"fuel at {(i, j)} faces `-`")
            if not faces or assemblies[load[i][j]]["design"] != "31":
                continue
            if len(faces) == 1:
                table, turn, kind = one, TURN_OF_I[faces[0]], "I"
            elif len(faces) == 2:
                pair = next((p for p in TURN_OF_CORNER if set(p) == set(faces)), None)
                if pair is None:
                    raise SystemExit(f"fuel at {(i, j)} faces {faces}")
                table, turn, kind = two, TURN_OF_CORNER[pair], "C"
            else:
                raise SystemExit(f"fuel at {(i, j)} faces {faces}")
            if table is None:
                continue
            design = f"{kind}-{Path(table).stem}"
            added[design] = table
            aid = load[i][j]
            assemblies[aid] = {"design": design}
            load[i][j] = aid + (f"@{turn}" if turn else "")
    blocks = [
        f"[design.{k}]\nbounds = {v['bounds']}\nlayers = [{', '.join(f'{chr(34)}{REL}/{x}{chr(34)}' for x in v['layers'])}]\n"
        for k, v in inv["design"].items()
    ]
    blocks += [
        f'[design.{k}]\nbounds = [0.0, 21.5]\nlayers = ["{t}"]\n'
        for k, t in added.items()
    ]
    blocks += [
        f'[assembly.{a}]\ndesign = "{v["design"]}"\n' for a, v in assemblies.items()
    ]
    (HERE / f"{name}-assemblies.toml").write_text("\n".join(blocks))
    width = max(len(c) for row in load for c in row) + 1
    rows = "".join(
        '  "' + "".join(c.ljust(width) for c in row).rstrip() + '",\n' for row in load
    )
    (HERE / f"{name}-cycle.toml").write_text(
        f'plant = "{REL}/{cycle["plant"]}"\ninventory = "{name}-assemblies.toml"\nloading = [\n{rows}]\n'
    )
    case = (P5 / f"core-{problem}.toml").read_text()
    case = case.replace(f"core-{problem}-cycle.toml", f"{name}-cycle.toml")
    case = case.replace(f"core-{problem}-pins.msgpack", f"{name}-pins.msgpack")
    (HERE / f"{name}.toml").write_text(case)


def main() -> None:
    (HERE / "tables").mkdir(exist_ok=True)
    for donor in ["I", "corner"]:
        for term, keys in TERMS.items():
            mix(donor, keys, HERE / "tables" / f"{donor}-{term}.msgpack")
    for term in ["xs", "D", "xsD", "adf"]:
        write_lattice_core(term)
    full_i = f"{REL}/lattice/reflector-I-fuel.msgpack"
    full_c = f"{REL}/lattice/reflector-corner-fuel.msgpack"
    for problem in ["5A", "5B", "5C", "5A-Inext", "5A-Ly", "small-E", "small-I"]:
        write_periphery(problem, f"{problem}-IC", full_i, full_c)
    for problem in ["5A", "small-E", "small-I"]:
        write_periphery(problem, f"{problem}-I", full_i, None)
        write_periphery(problem, f"{problem}-C", None, full_c)
        for term in TERMS:
            write_periphery(
                problem,
                f"{problem}-IC-{term}",
                f"tables/I-{term}.msgpack",
                f"tables/corner-{term}.msgpack",
            )


if __name__ == "__main__":
    main()
