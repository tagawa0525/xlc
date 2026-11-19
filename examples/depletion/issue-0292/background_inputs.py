"""#0292、ADR-0261: 無限希釈で数える形の残る差を分ける（compare_background.py）ための、同じ組成で輸送しない核種の扱いだけを変えた
ブランチ計算の入力を書く。

`python3 background_inputs.py` は、燃料ごと、状態（`CASES`）ごとに、background/<燃料>/<状態>.toml と、その履歴
background/<燃料>/out/<状態>-history.msgpack を書く。履歴は、組成の履歴（燃焼度ごとの数密度）に、扱いの履歴の輸送する核種
（`transportable`）と無限希釈の吸収体（`dilute_absorbers`）を当てたもの。lattice は、この 2 つがライブラリから求め直した組と
同じことを確かめるので、扱いの履歴はそのライブラリで燃焼させた履歴にする（どれも 565 K の燃焼で、組は燃料によらない）。

ブランチ計算は `BURNUP` の 1 点の、燃焼と同じ基準条件（減速材 0.743 g/cc、565 K、ホウ素 1300 ppm）の 1 状態で、少数群の境界は
ライブラリの群の境界そのもの、`micro` は輸送するすべての核種。Gd の 3×3 は燃焼の領域が複数なので、compare_background.py が
分けられず、入れない。
"""

import copy
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from inputs import NOMINAL_DENSITY, branch_input, read

BURNUP = 60.0
# 燃料ごとの状態: (名前, 組成の履歴, 扱いの履歴, ライブラリ)。履歴は issue-0292/<置き場>/nominal.msgpack、ライブラリは
# issue-0292/<置き場>/depletion.msgpack。並びは compare_background.py が隣り合う 2 つを比べる順:
# - ref-65: 65 核種を輸送して残りを数えた燃焼そのもの（MOX の比べる相手）
# - ref-69: 65 核種と外の核分裂の 4 核種を輸送して残りを数える扱い（U は比べる相手の燃焼そのもの）
# - form: 13 核種を輸送して残りを数える扱い（調べる形）を、ref-69 と同じ組成で
# - own-13: 調べる形の燃焼そのもの
CASES = {
    "u": [
        ("ref-69", "dilute/t69/u", "dilute/t69/u", "dilute/library-69"),
        ("form", "dilute/t69/u", "dilute/t13/u", "dilute/library-13"),
        ("own-13", "dilute/t13/u", "dilute/t13/u", "dilute/library-13"),
    ],
    "mox": [
        ("ref-65", "dilute/t65/mox", "dilute/t65/mox", "trial/library"),
        ("ref-69", "dilute/t65/mox", "dilute/t69/u", "dilute/library-69"),
        ("form", "dilute/t65/mox", "dilute/t13/mox", "dilute/library-13"),
        ("own-13", "dilute/t13/mox", "dilute/t13/mox", "dilute/library-13"),
    ],
}


def with_lists(composition: dict, lists: dict) -> dict:
    """`composition` の数密度に、`lists` の輸送する核種と無限希釈の吸収体を当てた履歴"""
    for h in (composition, lists):
        if h.get("format") != "xlc-depletion-history":
            raise SystemExit("not a lattice depletion history")
    if composition["nuclides"] != lists["nuclides"]:
        raise SystemExit("the histories list different nuclides")
    mixed = copy.deepcopy(composition)
    mixed["transportable"] = list(lists["transportable"])
    mixed["dilute_absorbers"] = list(lists.get("dilute_absorbers", []))
    return mixed


def case_input(
    fuel: str, name: str, library: str, bounds: list[float], micro: list[str]
) -> str:
    """background/<燃料>/ に置く状態 `name` の入力。`library` は examples/depletion からのパス"""
    note = (
        f"# #0292、ADR-0261: {fuel} の {name}（out/{name}-history.msgpack の組成と扱い）の {BURNUP:g} MWd/kg の基準条件。"
        "background_inputs.py が書く。手で直さない"
    )
    return branch_input(
        fuel,
        f"out/{name}",
        f"out/{name}-history.msgpack",
        [BURNUP],
        micro,
        NOMINAL_DENSITY,
        library=library,
        bounds=bounds,
        note=note,
    )


def main() -> None:
    for fuel, cases in CASES.items():
        directory = HERE / "background" / fuel
        (directory / "out").mkdir(parents=True, exist_ok=True)
        for name, composition, lists, library in cases:
            mixed = with_lists(
                read(HERE / composition / "nominal.msgpack"),
                read(HERE / lists / "nominal.msgpack"),
            )
            bounds = read(HERE / library / "depletion.msgpack")["group_bounds"]
            (directory / "out" / f"{name}-history.msgpack").write_bytes(
                msgpack.packb(mixed)
            )
            (directory / f"{name}.toml").write_text(
                case_input(
                    fuel,
                    name,
                    f"issue-0292/{library}/depletion.msgpack",
                    bounds,
                    mixed["transportable"],
                )
            )


if __name__ == "__main__":
    main()
