"""#0279 の段 2（ADR-0247）の k の感度の比較の入力を、2A の燃焼の履歴の 60 MWd/kg の組成から書く。

../p9-product/lattice-2A.msgpack（ADR-0239 の公称の状態で燃焼させた 2A。燃料は 5 つのリングの物質 fuel-0〜fuel-4。ADR-0245）の
最後の燃焼点（60 MWd/kg）の、最初の燃焼の領域の組（角のセル 20 の組）の 5 つのリングの数密度を、輸送する核種のうちどれかの
リングで正のものだけ燃料の物質に書く。負の数密度（燃焼の丸め。60 MWd/kg の C-12 の −1.4e-21 など）は、履歴からの分岐が解く物質
（crates/lattice/src/depletion.rs の `region_material`）と同じく 0 にする。被覆管、減速材、離散化、温度は ../issue-0180/branch/2A-60.toml と同じ。

- pincell.toml: 燃料のピン 1 本（鏡面反射）
- lattice-3x3.toml: 2A の案内管 1 本を燃料のピン 8 本で囲んだ 3×3 の格子（2A の小さい格子。鏡面反射）

- pincell-generalized.toml: pincell.toml の燃料のピンで、段 3 の一般化随伴の収束の方法を比べる（`generalized` の説明）

pincell.toml と lattice-3x3.toml は、31 核種の感度を、直接の摂動（ε = 0.1）と一緒に解く。MOC の判定は k 1e-12、中性子束 1e-11 で、#0279 の計画の
k 1e-10、中性子束 1e-9 より 2 桁細かい: 1e-10 では、感度の小さい核種（Gd-160 の S = −5.3e-7、Pm-148m の −8.7e-9）の直接の
摂動の k の変化が判定の誤差に近く、ピンセルで Gd-160 が 2.6e-2、Pm-148m が 5.5e-2 の相対の差になった（1e-12 では 8.1e-5、
4.6e-4。一般化摂動の値は両方の判定で 1.5e-10 しか変わらない）。

使い方: python3 inputs.py（このディレクトリで）
"""

import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
HISTORY = HERE / "../p9-product/lattice-2A.msgpack"
BURNUP = 60.0

# ADR-0035 の元素の輸送する核種 26 と、その輸送する先行核 5（ADR-0216。../issue-0180/branch/2A-60.toml の micro）
NUCLIDES = [
    "U-234", "U-235", "U-236", "U-238", "Np-237", "Np-239", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242",
    "Am-241", "Am-243", "Xe-135", "Sm-147", "Sm-149", "Sm-150", "Sm-151", "Sm-152",
    "Gd-152", "Gd-154", "Gd-155", "Gd-156", "Gd-157", "Gd-158", "Gd-160",
    "Cm-242", "Cm-244", "Eu-153", "I-135", "Pm-148m",
]  # fmt: skip

HEADER = """library = "../p2/depletion.msgpack"
ufg-tolerance = 1e-4
pitch = 1.26
azimuthal = 30
spacing = 0.03
"""

ZIRCALOY = """[material.zircaloy]
temperature = 580

[material.zircaloy.nuclides]
Cr-50 = 3.30121e-06
Cr-52 = 6.36606e-05
Cr-53 = 7.21860e-06
Cr-54 = 1.79686e-06
Fe-54 = 8.68307e-06
Fe-56 = 1.36306e-04
Fe-57 = 3.14789e-06
Fe-58 = 4.18926e-07
Zr-90 = 2.18865e-02
Zr-91 = 4.77292e-03
Zr-92 = 7.29551e-03
Zr-94 = 7.39335e-03
Zr-96 = 1.19110e-03
Sn-112 = 4.68066e-06
Sn-114 = 3.18478e-06
Sn-115 = 1.64064e-06
Sn-116 = 7.01616e-05
Sn-117 = 3.70592e-05
Sn-118 = 1.16872e-04
Sn-119 = 4.14504e-05
Sn-120 = 1.57212e-04
Sn-122 = 2.23417e-05
Sn-124 = 2.79392e-05

[material.water]
temperature = 580
water = { density = 0.711901, boron = 650, b10 = 0.199 }
"""

PIN = """[cell.pin]
zones = [
  { material = "fuel-0", radius = 0.18317868871678278, rings = 1 },
  { material = "fuel-1", radius = 0.25905378592099365, rings = 1 },
  { material = "fuel-2", radius = 0.3172747957213116, rings = 1 },
  { material = "fuel-3", radius = 0.36635737743356556, rings = 1 },
  { material = "fuel-4", radius = 0.4096, rings = 1 },
  { material = "void", radius = 0.418, rings = 1 },
  { material = "zircaloy", radius = 0.475, rings = 1 },
  { material = "water", radius = 0.63, rings = 3 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8
"""

GUIDE_TUBE = """[cell.gt]
zones = [
  { material = "water", radius = 0.561, rings = 3 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8
"""


def sensitivity(output: str) -> str:
    names = ", ".join(f'"{n}"' for n in NUCLIDES)
    return (
        "[sensitivity]\n"
        f"nuclides = [{names}]\n"
        f'output = "{output}"\n'
        "direct = 0.1\n"
        "k-tolerance = 1e-12\n"
        "flux-tolerance = 1e-11\n"
    )


RESPONSES = [
    "total-1", "total-2", "nu-fission-1", "nu-fission-2", "fission-1", "fission-2",
    "scatter-1-1", "scatter-1-2", "scatter-2-1", "scatter-2-2",
]  # fmt: skip


def generalized(output: str) -> str:
    """段 3 の一般化随伴の比較（ADR-0247）。少数群定数の微分の 10 項目の応答（2 群）を、3 つの収束の方法の 2 つの判定（GMRES は
    残差の 2 ノルムで止めるので、同じ判定では誤差が大きい）と、細かい判定の GMRES の参照で解く。CMFD と GMRES は 10 本を一緒に
    スイープしても解く。GMRES の再始動は 100 回（ピンセルの判定 1e-6 で、30 回の 1,474 スイープ（手元の試行。xlc d32a6fd8e をコミットする前の同じ内容の作業ツリー）に対し 927（Slurm 773））。MOC の判定は既定（本番の判定）"""
    names = ", ".join(f'"{n}"' for n in NUCLIDES)
    responses = ", ".join(f'"{r}"' for r in RESPONSES)
    return (
        "[sensitivity]\n"
        f"nuclides = [{names}]\n"
        f'output = "{output}"\n\n'
        "[sensitivity.generalized]\n"
        "bounds = [1e7, 0.625, 1e-5]\n"
        f"responses = [{responses}]\n"
        'methods = ["source-iteration", "cmfd", "gmres"]\n'
        "restart = 100\n"
        "tolerances = [1e-6, 1e-7]\n"
        'batched = ["cmfd", "gmres"]\n'
        'reference = "gmres"\n'
        "reference-tolerance = 1e-10\n"
    )


def fuels(history: dict) -> str:
    """最後の燃焼点の、最初の燃焼の領域の組の 5 つのリングの物質（自己遮蔽はリングの核種が同じ順であることを求めるので、どれかの
    リングで正の核種をすべてのリングに同じ順で書く）"""
    point = history["points"][-1]
    if point["burnup"] != BURNUP:
        sys.exit(f"the last point is {point['burnup']} MWd/kg, not {BURNUP}")
    transportable = set(history["transportable"])
    first = history["regions"][0]["cell"]
    rings = [
        (region["material"], densities)
        for region, densities in zip(history["regions"], point["densities"])
        if region["cell"] == first
    ]
    if len(rings) != 5:
        sys.exit("the first region group must have the five fuel rings")
    kept = [
        i
        for i, name in enumerate(history["nuclides"])
        if name in transportable and any(d[i] > 0.0 for _, d in rings)
    ]
    lines = []
    for material, densities in rings:
        lines.append(
            f"[material.{material}]\ntemperature = 900\n\n[material.{material}.nuclides]"
        )
        for i in kept:
            lines.append(f"{history['nuclides'][i]} = {max(densities[i], 0.0)!r}")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    with HISTORY.open("rb") as f:
        history = msgpack.unpackb(f.read(), strict_map_key=False)
    materials = fuels(history) + "\n" + ZIRCALOY
    pincell = (
        HEADER
        + 'rows = ["pin"]\n\n'
        + materials
        + "\n"
        + PIN
        + "\n"
        + sensitivity("pincell-sensitivity.toml")
    )
    lattice = (
        HEADER
        + 'rows = ["pin pin pin", "pin gt  pin", "pin pin pin"]\n\n'
        + materials
        + "\n"
        + PIN
        + "\n"
        + GUIDE_TUBE
        + "\n"
        + sensitivity("lattice-3x3-sensitivity.toml")
    )
    pincell_generalized = (
        HEADER
        + 'rows = ["pin"]\n\n'
        + materials
        + "\n"
        + PIN
        + "\n"
        + generalized("pincell-generalized-sensitivity.toml")
    )
    comment = "# このファイルは inputs.py が書く（手で直さない）。#0279 の段 2 の k の感度の比較（ADR-0247）\n"
    (HERE / "pincell.toml").write_text(comment + pincell)
    (HERE / "lattice-3x3.toml").write_text(comment + lattice)
    (HERE / "pincell-generalized.toml").write_text(
        "# このファイルは inputs.py が書く（手で直さない）。#0279 の段 3 の一般化随伴の比較（ADR-0247）\n"
        + pincell_generalized
    )


if __name__ == "__main__":
    main()
