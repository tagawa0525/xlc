"""燃焼度の軸を持つ分岐のテーブル（`xlc-table`、ADR-0187）から、燃焼度の点を選んだテーブルを書く（ADR-0192 の、点を間引いた
テーブルで core を実行する確認）。ほかの軸、状態の値、`heavy_metal` などの項目はそのまま。

使い方: python3 thin_burnup_table.py <テーブル（msgpack）> <出力（msgpack）> <燃焼度 (MWd/kg)>...

燃焼度はテーブルの燃焼度の点と厳密に等しく、増える順であること（違えばエラー）。燃焼度は最初の軸で、解いた点（`points`）の
うち選んだ燃焼度の点のものを、燃焼度の番号を付け直して残す（点ごとに部分ノードの数（split²）の状態が続く）。点ごとに解いた燃焼度の
点が違うテーブル（ADR-0200）も同じ。
"""

import sys
from itertools import pairwise
from pathlib import Path

import msgpack


def thin(table: dict, burnups: list[float]) -> dict:
    """燃焼度の点 `burnups` の状態だけを残したテーブル"""
    axes = table["axes"]
    if not axes or axes[0]["name"] != "burnup":
        raise SystemExit("the first axis of the table must be the burnup")
    points = axes[0]["points"]
    if not burnups or any(b not in points for b in burnups):
        raise SystemExit(f"the burnups {burnups} must be points of the table {points}")
    if any(a >= b for a, b in pairwise(burnups)):
        raise SystemExit(f"the burnups {burnups} must increase")
    sub = table["split"] ** 2
    if len(table["states"]) != len(table["points"]) * sub:
        raise SystemExit(
            f"the table has {len(table['states'])} states for {len(table['points'])} points"
        )
    index = {points.index(b): n for n, b in enumerate(burnups)}
    kept_points = []
    states = []
    for n, point in enumerate(table["points"]):
        if point[0] in index:
            kept_points.append([index[point[0]], *point[1:]])
            states += table["states"][n * sub : (n + 1) * sub]
    return {
        **table,
        "axes": [{**axes[0], "points": list(burnups)}, *axes[1:]],
        "points": kept_points,
        "states": states,
    }


def main() -> None:
    args = sys.argv[1:]
    if len(args) < 3:
        raise SystemExit(__doc__)
    table = msgpack.unpackb(Path(args[0]).read_bytes())
    out = thin(table, [float(a) for a in args[2:]])
    Path(args[1]).write_bytes(msgpack.packb(out))


if __name__ == "__main__":
    main()
