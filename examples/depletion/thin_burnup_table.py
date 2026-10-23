"""燃焼度の軸を持つ分岐のテーブル（`xlc-table`、ADR-0187）から、燃焼度の点を選んだテーブルを書く（ADR-0192 の、点を間引いた
テーブルで core を実行する確認）。ほかの軸、状態の値、`heavy_metal` などの項目はそのまま。

使い方: python3 thin_burnup_table.py <テーブル（msgpack）> <出力（msgpack）> <燃焼度 (MWd/kg)>...

燃焼度はテーブルの燃焼度の点と厳密に等しく、増える順であること（違えばエラー）。燃焼度は最初の軸（いちばん遅く変わる）で、
状態は燃焼度の点ごとに、ほかの軸の点の数の積 × 部分ノードの数（split²）だけ続く。
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
    block, rest = divmod(len(table["states"]), len(points))
    if rest or block == 0:
        raise SystemExit(
            f"the table has {len(table['states'])} states for {len(points)} burnups"
        )
    states = []
    for b in burnups:
        i = points.index(b)
        states += table["states"][i * block : (i + 1) * block]
    return {
        **table,
        "axes": [{**axes[0], "points": list(burnups)}, *axes[1:]],
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
