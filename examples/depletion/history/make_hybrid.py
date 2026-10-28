"""公称の燃焼の履歴の組成のうち、ミクロ燃焼で追う元素の核種だけを、条件を変えた燃焼の履歴の値に置き換えた履歴を書く。

    python3 make_hybrid.py <公称の履歴> <条件を変えた履歴> <出力>

履歴は lattice の燃焼の履歴（`xlc-depletion-history`。ADR-0173）。2 つの履歴は、同じ核種の並び、同じ燃焼の領域、同じ燃焼点を
持つこと（違えば止まる。燃焼の領域は物質、セル、体積で比べる）。置き換えるのは各燃焼点の各領域の数密度だけで、k∞ などのほかの項目は公称の履歴のまま
（`[history]` の分岐は数密度だけを読む。ADR-0187）。追う元素は ADR-0035 の候補（U、Np、Pu、Am、Xe、Sm、Gd）。
"""

import sys

import msgpack

TRACKED = ("U", "Np", "Pu", "Am", "Xe", "Sm", "Gd")


def element(nuclide):
    """`Pu-239`、`Am-242m` の元素の記号"""
    return nuclide.split("-")[0]


def identity(region):
    """燃焼の領域の物質、セル、体積。MOC の領域の番号は格子（案内管の中身）で変わり、lattice は入力から求め直すので比べない（ADR-0187）"""
    return (region["material"], region["cell"], region["volume"])


def tracked_indices(nuclides):
    return [i for i, n in enumerate(nuclides) if element(n) in TRACKED]


def hybrid(nominal, varied):
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
    indices = tracked_indices(nominal["nuclides"])
    for point, other in zip(nominal["points"], varied["points"]):
        for densities, source in zip(point["densities"], other["densities"]):
            for i in indices:
                densities[i] = source[i]
    return nominal


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    nominal, varied, output = sys.argv[1:]
    with open(output, "wb") as f:
        msgpack.pack(hybrid(read(nominal), read(varied)), f)


if __name__ == "__main__":
    main()
