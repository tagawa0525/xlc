"""公称の燃焼の履歴の組成のうち、ミクロ燃焼で追う核種だけを、条件を変えた燃焼の履歴の値に置き換えた履歴を書く。

    python3 make_hybrid.py <公称の履歴> <条件を変えた履歴> <出力>
    python3 make_hybrid.py <公称の履歴> <条件を変えた履歴> <出力> --micro <分岐の入力>

履歴は lattice の燃焼の履歴（`xlc-depletion-history`。ADR-0173）。2 つの履歴は、同じ核種の並び、同じ燃焼の領域、同じ燃焼点を
持つこと（違えば止まる。燃焼の領域は物質、セル、体積で比べる）。置き換えるのは各燃焼点の各領域の数密度だけで、k∞ などのほかの項目は公称の履歴のまま
（`[history]` の分岐は数密度だけを読む。ADR-0187）。置き換える核種は、`--micro` が無ければ ADR-0035 の候補の元素（U、Np、Pu、Am、
Xe、Sm、Gd）の核種のすべて（「履歴の効き（2A）」）、あれば分岐の入力（TOML）の `[history]` の `micro` の核種（履歴の残り R の第 3 項と
同じ追う核種。ADR-0215、ADR-0216）。
"""

import argparse

import msgpack
import tomllib

TRACKED = ("U", "Np", "Pu", "Am", "Xe", "Sm", "Gd")


def element(nuclide):
    """`Pu-239`、`Am-242m` の元素の記号"""
    return nuclide.split("-")[0]


def identity(region):
    """燃焼の領域の物質、セル、体積。MOC の領域の番号は格子（案内管の中身）で変わり、lattice は入力から求め直すので比べない（ADR-0187）"""
    return (region["material"], region["cell"], region["volume"])


def tracked_indices(nuclides):
    return [i for i, n in enumerate(nuclides) if element(n) in TRACKED]


def listed_indices(nuclides, listed):
    missing = [n for n in listed if n not in nuclides]
    if missing:
        raise SystemExit(f"the histories have no {', '.join(missing)}")
    return [nuclides.index(n) for n in listed]


def micro_nuclides(path):
    """分岐の入力の `[history]` の `micro`（追う核種）"""
    with open(path, "rb") as f:
        micro = tomllib.load(f).get("history", {}).get("micro")
    if not micro:
        raise SystemExit(f"{path}: the input has no [history] micro")
    return micro


def hybrid(nominal, varied, listed=None):
    """`listed` が None なら追う元素の核種、そうでなければ `listed` の核種を置き換える"""
    if (
        nominal["format"] != "xlc-depletion-history"
        or varied["format"] != "xlc-depletion-history"
    ):
        raise SystemExit("both files must be xlc-depletion-history")
    if nominal["nuclides"] != varied["nuclides"]:
        raise SystemExit("the histories list different nuclides")
    if [identity(r) for r in nominal["regions"]] != [
        identity(r) for r in varied["regions"]
    ]:
        raise SystemExit("the histories have different depletion regions")
    if [p["burnup"] for p in nominal["points"]] != [
        p["burnup"] for p in varied["points"]
    ]:
        raise SystemExit("the histories have different burnup points")
    indices = (
        tracked_indices(nominal["nuclides"])
        if listed is None
        else listed_indices(nominal["nuclides"], listed)
    )
    for point, other in zip(nominal["points"], varied["points"]):
        for densities, source in zip(point["densities"], other["densities"]):
            for i in indices:
                densities[i] = source[i]
    return nominal


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("nominal")
    parser.add_argument("varied")
    parser.add_argument("output")
    parser.add_argument(
        "--micro", help="replace the [history] micro nuclides of this branch input"
    )
    args = parser.parse_args()
    listed = None if args.micro is None else micro_nuclides(args.micro)
    with open(args.output, "wb") as f:
        msgpack.pack(hybrid(read(args.nominal), read(args.varied), listed), f)


if __name__ == "__main__":
    main()
