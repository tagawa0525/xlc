"""格子計算の燃焼の検証（VERA Problem 2 の集合体。ADR-0178）の入力を書く: OpenMC の燃焼（p2/openmc-<ジョブ>.toml）と、xlc の
格子計算の燃焼（p2/lattice-<ジョブ>.toml）。同じ組成、形状、燃焼度の点を 1 か所（openmc_depletion_assembly.py の表）から作り、
2 つの入力がずれないようにする。このファイルが書くものは手で直さない。

使い方（references の devShell で）:

    python3 p2_inputs.py          # 入力を書く
    python3 p2_inputs.py --check  # 書かれている入力が、この表から作ったものと同じか確かめる
"""

import sys
from pathlib import Path

import openmc_depletion_assembly as asm

P2 = Path(__file__).resolve().parent / "p2"

# 燃焼の出力密度 (W/g)。VERA Problem 9 のサイクル 1 の平均（16.939 GWd/MT / 441 EFPD）。ピンセルの検証（openmc-1C-rings1.toml）と同じ
POWER = 38.41
# 1 世代の粒子数、世代数、捨てる世代数。数える 60 世代 × 100 万 = 6×10⁷ 粒子で、対称操作の組で平均したピン出力の統計誤差が約 0.065%
# （ADR-0032 の目標 rms 0.2% の約 1/3）。2A の 2 万粒子・数える 20 世代（4×10⁵）の短い走行で 0.80%（組で平均した rms。新しい燃料と
# 燃焼した燃料で 0.80% と 0.84%）で、誤差は粒子数の平方根に反比例する（ADR-0178）
PARTICLES = 1_000_000
BATCHES = 100
INACTIVE = 40
SEED = 1
THREADS = 16
# Gd 燃料のリングの数（等しい面積。数は Gd で測って決める。docs/plans/m4-depletion.md の「格子計算の燃焼」）
GD_RINGS = 8

# ジョブ名 → (問題, 燃焼度の点, 説明)。ジョブ名は AWS のジョブ名（examples/depletion/p2/aws.toml の jobs）の openmc-<ジョブ>
JOBS = {
    "2A": (
        "2A",
        asm.BURNUP_40,
        "毒物のない 3.1% の集合体（燃料、被覆管、減速材 565 K、減速材 0.743 g/cc、ホウ素 1300 ppm。仕様書の表 P2-1）",
    ),
    "2F": (
        "2F",
        asm.BURNUP_40,
        "2A の案内管 24 本すべてに Pyrex の棒を入れたもの（600 K、0.743 g/cc。表 P2-1。Pyrex も燃焼する）",
    ),
    "2M": (
        "2M",
        asm.BURNUP_70,
        "IFBA 128 本（図 P2-2。ZrB2 の被覆の B-10 が燃焼する。600 K、0.743 g/cc）",
    ),
    "2N": (
        "2N",
        asm.BURNUP_40,
        "IFBA 104 本と WABA 20 本（図 P2-2。どちらも燃焼する。600 K、0.743 g/cc）",
    ),
    "2O": (
        "2O",
        asm.BURNUP_70,
        "Gd 燃料（1.8% 濃縮、5% Gd2O3）12 本（図 P2-2。リングごとに燃焼する。600 K、0.743 g/cc）",
    ),
    "2P": (
        "2P",
        asm.BURNUP_70,
        "Gd 燃料 24 本（図 P2-2。リングごとに燃焼する。600 K、0.743 g/cc）",
    ),
    "2A-density": (
        "2A-density",
        asm.BURNUP_40,
        "履歴を変えた燃焼 1: 2A の減速材の密度を 0.661 g/cc に（仕様書の 2B〜2D の密度。温度は 565 K のまま）",
    ),
    "2A-rods": (
        "2A-rods",
        asm.BURNUP_40,
        "履歴を変えた燃焼 2: 2A の案内管 24 本すべてに AIC の制御棒を入れたまま燃やす（2G と同じ配置。制御棒は燃焼しない）",
    ),
}

# 時間の倍率 s を測る短い走行（2A。粒子 2 万・世代 40、手元の Slurm の `-c 8`）。初めの組成は新しい燃料と、ピンセルの燃焼
# （../openmc-1C-rings1.msgpack。燃焼度の点 [0, 0.1, 0.5, 1, 2, 3, 5] の 5 MWd/kg、番号 6）の燃料の組成
S_PARTICLES = 20_000
S_BATCHES = 40
S_INACTIVE = 20
S_THREADS = 8
S_BURNUP = [0.1]
S_RUNS = {"s-2A-fresh": None, "s-2A-burned": ("../openmc-1C-rings1.msgpack", 6)}

# セルの種類 → lattice の入力のセルの名前
CELL_NAMES = {
    "fuel": "pin",
    "ifba": "ib",
    "gd": "gd",
    "gt": "gt",
    "it": "it",
    "pyrex": "py",
    "waba": "wb",
    "aic": "cr",
}
CELL_COMMENTS = {
    "pin": "燃料のピン（3.1%）",
    "ib": "IFBA の燃料のピン（ペレットの外側に ZrB2 の被覆 10 μm。ペレットと被覆管の間は真空）",
    "gd": "Gd 燃料のピン（リングごとに別の物質。ギャップは真空）",
    "gt": "案内管（中は減速材）",
    "it": "計装管（中は減速材）",
    "py": "Pyrex の棒を入れた案内管（内管の中と、Pyrex の内外の隙間は真空）",
    "wb": "WABA の棒を入れた案内管（内管の中は減速材、内管と吸収体の間、吸収体と被覆管の間は真空）",
    "cr": "AIC の制御棒を入れた案内管（吸収体と被覆管の間は真空）",
}


def xlc_name(nuclide: str) -> str:
    """OpenMC の核種名（`Zr90`）を xlc の名前（`Zr-90`）にする。"""
    head = nuclide.rstrip("0123456789")
    return f"{head}-{nuclide[len(head) :]}"


def number(value: float) -> str:
    return f"{value:.5e}"


def material_block(name: str, densities: dict[str, float], temperature: float) -> str:
    entries = [f"{xlc_name(n)} = {number(d)}" for n, d in densities.items()]
    if len(entries) <= 6:
        return (
            f"[material.{name}]\ntemperature = {temperature:g}\n"
            f"nuclides = {{ {', '.join(entries)} }}\n"
        )
    body = "\n".join(entries)
    return (
        f"[material.{name}]\ntemperature = {temperature:g}\n\n"
        f"[material.{name}.nuclides]\n{body}\n"
    )


def radius(value: float | None) -> str:
    return "" if value is None else f"radius = {value:.10g}, "


def cell_block(kind: str, gd_rings: int | None) -> str:
    name = CELL_NAMES[kind]
    zones = ",\n".join(
        f'  {{ material = "{z.material}", {radius(z.radius)}rings = {z.rings} }}'
        for z in asm.cell_zones(kind, gd_rings)
    )
    return (
        f"# {CELL_COMMENTS[name]}\n[cell.{name}]\nzones = [\n{zones},\n]\nsectors = 8\n"
    )


def burnup_lines(burnup: list[float]) -> str:
    chunks = [burnup[i : i + 10] for i in range(0, len(burnup), 10)]
    lines = ",\n".join("  " + ", ".join(f"{b:g}" for b in chunk) for chunk in chunks)
    return f"[\n{lines},\n]"


def lattice_input(job: str) -> str:
    problem, burnup, description = JOBS[job]
    c = asm.CASES[problem]
    run = {"problem": problem, **({"gd-rings": GD_RINGS} if c.gd else {})}
    gd_rings = run.get("gd-rings")
    kinds = [[asm.cell_kind(problem, i, j) for j in range(asm.N)] for i in range(asm.N)]
    used = [k for k in CELL_NAMES if any(k in row for row in kinds)]
    width = max(len(CELL_NAMES[k]) for k in used)
    rows = ",\n".join(
        '  "' + " ".join(CELL_NAMES[k].ljust(width) for k in row) + '"' for row in kinds
    )
    region_kinds = list(dict.fromkeys(r.kind for r in asm.regions(run)))
    names = ["fuel", "zircaloy", "water"]
    fixed = {n for k in used for z in asm.cell_zones(k, gd_rings) for n in [z.material]}
    names += [n for n in ("ss304", "pyrex", "waba", "ifba", "aic") if n in fixed]
    names += sorted(
        (n for n in fixed if n.startswith("gd") and n[2:].isdigit()),
        key=lambda n: int(n[2:]),
    )
    blocks = [
        material_block(n, asm.composition(n, c), c.temperature)
        for n in dict.fromkeys(names)
    ]
    cells = "\n".join(cell_block(k, gd_rings) for k in used)
    return (
        f"# 格子計算の燃焼の検証（ADR-0178）: VERA Problem 2 の {description}を、分岐なしの 1 つの状態で、\n"
        f"# {len(burnup)} 点（最後は {burnup[-1]:g} MWd/kg）まで燃やす。OpenMC の参照解は openmc-{job}.toml。\n"
        "# 出力密度 38.41 W/g は VERA Problem 9 のサイクル 1 の平均（16.939 GWd/MT / 441 EFPD）。形状、組成と燃焼度の点は\n"
        "# OpenMC の入力と同じ表（openmc_depletion_assembly.py）から p2_inputs.py が書く。手で直さない。集合体の周りに水のギャップ 0.04 cm。\n"
        f"# 燃やす物質は {', '.join(region_kinds)}（対称操作の組ごとに 1 つの領域。ADR-0173）。\n"
        'library = "depletion.msgpack"\n'
        "pitch = 1.26\n"
        "azimuthal = 30\n"
        "spacing = 0.03\n"
        'gap = { width = 0.04, material = "water" }\n'
        f"rows = [\n{rows},\n]\n\n"
        + "\n".join(blocks)
        + "\n"
        + cells
        + "\n[depletion]\n"
        'chain = "../chain.msgpack"\n'
        f"power = {POWER:g}\n"
        f"burnup = {burnup_lines(burnup)}\n"
        f"materials = [{', '.join(chr(34) + k + chr(34) for k in region_kinds)}]\n"
        f'output = "lattice-{job}.msgpack"\n'
    )


def openmc_input(job: str) -> str:
    problem, burnup, description = JOBS[job]
    gd = f"gd-rings = {GD_RINGS}\n" if asm.CASES[problem].gd else ""
    return (
        f"# 格子計算の燃焼の検証（ADR-0178）の参照解: VERA Problem 2 の {description}を、{len(burnup)} 点（最後は {burnup[-1]:g} MWd/kg）\n"
        "# まで OpenMC の燃焼で燃やす（openmc_depletion_assembly.py。AWS で流す。aws.toml）。xlc は lattice-"
        + job
        + ".toml。\n"
        "# 出力密度 38.41 W/g は VERA Problem 9 のサイクル 1 の平均。粒子は 1 世代 "
        f"{PARTICLES:,}・{BATCHES} 世代（最初の {INACTIVE} 世代を捨てる）: ピン出力の\n"
        "# 統計誤差を ADR-0032 の目標（rms 0.2%）の約 1/3（約 0.065%）にする数（ADR-0178）。p2_inputs.py が書く。手で直さない。\n"
        f'problem = "{problem}"\n'
        'cross-sections = "openmc-library/cross_sections.xml"\n'
        'chain = "../openmc-chain.xml"\n'
        f"power = {POWER:g}\n"
        f"burnup = {burnup_lines(burnup)}\n"
        f"particles = {PARTICLES}\n"
        f"batches = {BATCHES}\n"
        f"inactive = {INACTIVE}\n"
        f"seed = {SEED}\n"
        f"threads = {THREADS}\n" + gd
    )


def s_input(name: str) -> str:
    start = S_RUNS[name]
    source = (
        "燃料の初めの組成は新しい燃料"
        if start is None
        else f"燃料の初めの組成はピンセルの燃焼（{start[0]}）の点 {start[1]}（5 MWd/kg）"
    )
    extra = (
        ""
        if start is None
        else f'start-from = "{start[0]}"\nstart-point = {start[1]}\n'
    )
    return (
        "# 時間の倍率 s（燃焼した燃料と燃焼の tally による輸送の時間の倍率。docs/plans/m4-depletion.md の「OpenMC の燃焼の参照解の\n"
        f"# 物理コア時間」）を測る短い走行（ADR-0178）: 2A の 1 ステップ（{S_BURNUP[0]:g} MWd/kg、輸送 3 回）、粒子 {S_PARTICLES:,}・{S_BATCHES} 世代\n"
        f"# （最初の {S_INACTIVE} 世代を捨てる）、{S_THREADS} スレッド。{source}。openmc_depletion_assembly.py の run と transport で解く\n"
        "# （手元の Slurm で `-c 8`。時間は記録に書かない）。p2_inputs.py が書く。手で直さない。\n"
        'problem = "2A"\n'
        'cross-sections = "openmc-library/cross_sections.xml"\n'
        'chain = "../openmc-chain.xml"\n'
        f"power = {POWER:g}\n"
        f"burnup = {burnup_lines(S_BURNUP)}\n"
        f"particles = {S_PARTICLES}\n"
        f"batches = {S_BATCHES}\n"
        f"inactive = {S_INACTIVE}\n"
        f"seed = {SEED}\n"
        f"threads = {S_THREADS}\n" + extra
    )


def files() -> dict[Path, str]:
    out = {P2 / f"{name}.toml": s_input(name) for name in S_RUNS}
    for job in JOBS:
        out[P2 / f"openmc-{job}.toml"] = openmc_input(job)
        out[P2 / f"lattice-{job}.toml"] = lattice_input(job)
    return out


def main(argv: list[str]) -> None:
    match argv:
        case []:
            for path, text in files().items():
                path.write_text(text)
        case ["--check"]:
            stale = [
                str(p)
                for p, t in files().items()
                if not p.exists() or p.read_text() != t
            ]
            if stale:
                raise SystemExit(f"stale: {stale}; run p2_inputs.py")
        case _:
            raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
