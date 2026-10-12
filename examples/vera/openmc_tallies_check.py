"""標準のタリーの組（openmc_tallies.py、ADR-0152）の定義を固定する確認。

使い方: python3 examples/vera/openmc_tallies_check.py（リポジトリの直下で、references の devShell の中で）

問題の種類ごとに、タリーの名前、フィルタの種類と並び、スコア、核種が下の表と一致することを確かめる。組を変えて
この確認が落ちたら、`VERSION` を上げてから表を直す（版が同じなら中身も同じ、を保つため）。
"""

import re
import sys
from pathlib import Path

import openmc

sys.path.insert(0, str(Path(__file__).resolve().parent))
import openmc_tallies as std

ROOT = Path(__file__).resolve().parents[2]

MATERIALS = ("material groups", "material transfer", "material fine", "nuclides")
PIN = ("pin groups", "pin nuclides", "pin nuclides flux", "boundary")
LAYER = ("layer groups", "layer transfer", "layer current")
ASSEMBLY = ("assembly groups", "assembly transfer", "assembly current")

M, E, EO, L, F, MS = (
    "MaterialFilter",
    "EnergyFilter",
    "EnergyoutFilter",
    "LegendreFilter",
    "MeshFilter",
    "MeshSurfaceFilter",
)
ZONE_SCORES = ["flux", "total", "scatter", "absorption", "fission", "nu-fission"]
MULTIPLYING = [
    "(n,2nd)",
    "(n,2n)",
    "(n,3n)",
    "(n,2na)",
    "(n,3na)",
    "(n,2n2a)",
    "(n,4n)",
    "(n,2np)",
    "(n,3np)",
]
NUCLIDE_SCORES = [
    "total",
    "elastic",
    "absorption",
    "fission",
    "nu-fission",
    "(n,gamma)",
    "(n,2n)",
]
# 名前 → (フィルタの並び, スコア, 核種)。核種は下の小さなモデル（燃料、減速材、被覆管、AIC）にあるもの
TRACKED = ["U235", "U238", "O16", "H1", "B10", "Ag107", "In115"]
FUEL_NUCLIDES = ["U235", "U238", "O16"]
EXPECTED = {
    "material groups": ([M, E], ZONE_SCORES + MULTIPLYING, []),
    "material transfer": ([M, E, EO, L], ["nu-scatter"], []),
    "material fine": ([M, E], ["flux", "absorption", "nu-fission"], []),
    "nuclides": ([M, E], NUCLIDE_SCORES, TRACKED),
    "pin groups": ([F, E], ZONE_SCORES, []),
    "pin nuclides": ([F, M, E], NUCLIDE_SCORES, FUEL_NUCLIDES),
    "pin nuclides flux": ([F, M, E], ["flux"], []),
    "boundary": ([F, E], ["flux"], []),
    "layer groups": ([F, E], ZONE_SCORES, []),
    "layer transfer": ([F, E, EO], ["nu-scatter"], []),
    "layer current": ([MS, E], ["current"], []),
    "assembly groups": ([F, E], ZONE_SCORES, []),
    "assembly transfer": ([F, E, EO], ["nu-scatter"], []),
    "assembly current": ([MS, E], ["current"], []),
}


def model_materials() -> tuple[list[openmc.Material], list[openmc.Material]]:
    fuel = openmc.Material(name="fuel")
    for n, d in (("U235", 7e-4), ("U238", 2.2e-2), ("O16", 4.6e-2)):
        fuel.add_nuclide(n, d)
    water = openmc.Material(name="moderator")
    for n, d in (("H1", 5e-2), ("O16", 2.5e-2), ("B10", 1e-5)):
        water.add_nuclide(n, d)
    clad = openmc.Material(name="zircaloy-4")
    clad.add_nuclide("Zr90", 2e-2)
    aic = openmc.Material(name="aic")
    for n, d in (("Ag107", 2e-2), ("In115", 8e-3)):
        aic.add_nuclide(n, d)
    return [fuel, water, clad, aic], [fuel]


def pins(n: int, pitch: float) -> openmc.RegularMesh:
    mesh = openmc.RegularMesh()
    mesh.dimension = (n, n)
    mesh.lower_left = (-n * pitch / 2, -n * pitch / 2)
    mesh.upper_right = (n * pitch / 2, n * pitch / 2)
    return mesh


def check_set(label: str, tallies: list[openmc.Tally], kinds: tuple[str, ...]) -> None:
    names = [t.name for t in tallies]
    expected = [f"standard-{std.VERSION}/{k}" for k in kinds]
    assert names == expected, (label, names, expected)
    for t, kind in zip(tallies, kinds):
        filters, scores, nuclides = EXPECTED[kind]
        got = [type(f).__name__ for f in t.filters]
        assert got == filters, (label, kind, got, filters)
        assert list(t.scores) == scores, (label, kind, list(t.scores), scores)
        got_nuclides = [n for n in t.nuclides if n != "total"]
        assert got_nuclides == nuclides, (label, kind, got_nuclides, nuclides)
        for f in t.filters:
            if isinstance(f, (openmc.EnergyFilter, openmc.EnergyoutFilter)):
                assert len(f.values) in (len(std.XLC_BINS), len(std.FINE_BINS))
            if isinstance(f, openmc.LegendreFilter):
                assert f.order == 3


def check_wims69() -> None:
    """群の境界が xlc の WIMS 69 群（crates/format/src/group.rs）と同じこと。"""
    source = (ROOT / "crates/format/src/group.rs").read_text()
    body = re.search(r"WIMS69_ASCENDING: \[f64; 70\] = \[(.*?)\];", source, re.DOTALL)
    assert body is not None
    rust = [float(v) for v in body.group(1).replace("\n", " ").split(",") if v.strip()]
    assert list(std.WIMS69) == rust
    assert std.XLC_BINS == (0.0, *rust, 1.0e10)
    assert len(std.FINE_BINS) == 1969 and std.FINE_BINS[0] > 0


def check_strips() -> None:
    """集合体の外周の帯と、1/4 炉心の集合体の境の帯の境界。"""
    half, gap = 10.75, 0.04
    assert std.strip_edges(-half, [-half, half], half, gap) == [
        -half,
        -half + gap,
        half - gap,
        half,
    ]
    edges = std.strip_edges(0.0, [10.75, 32.25], 32.25 + gap, gap)
    expected = [0.0, 10.75 - gap, 10.75 + gap, 32.25 - gap, 32.25 + gap]
    assert all(abs(a - b) < 1e-12 for a, b in zip(edges, expected, strict=True))


def check_enabled() -> None:
    assert std.enabled({}) is True
    assert std.enabled({"standard-tallies": False}) is False
    try:
        std.enabled({"standard-tallies": "no"})
    except SystemExit:
        pass
    else:
        raise AssertionError("a non-boolean standard-tallies is accepted")


def main() -> None:
    check_wims69()
    check_strips()
    check_enabled()
    materials, fuels = model_materials()
    check_set("pincell", std.pincell(materials), MATERIALS)
    check_set(
        "assembly",
        std.assembly(materials, fuels, pins(17, 1.26), 10.75, 0.04),
        (*MATERIALS, *PIN),
    )
    check_set(
        "assembly 3d",
        std.assembly_3d(materials, fuels, pins(17, 1.26), 10.75, 0.04, [0.0, 1.0, 2.0]),
        (*MATERIALS, *PIN[:3], *LAYER, PIN[3]),
    )
    core = openmc.RegularMesh()
    core.dimension = (2, 2)
    core.lower_left = (-10.75, -32.25)
    core.upper_right = (32.25, 10.75)
    check_set(
        "core 2d", std.core_2d(materials, core, 0.04), (*MATERIALS, *ASSEMBLY, PIN[3])
    )
    print(f"standard-{std.VERSION}: ok")


if __name__ == "__main__":
    main()
