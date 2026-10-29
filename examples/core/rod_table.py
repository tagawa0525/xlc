"""制御棒の値ごとの格子計算のテーブルを、制御棒の離散の軸（ADR-0097、ADR-0127）を持つ 1 つのテーブルにまとめる。

格子計算が制御棒の分岐を作れるようになるまで（#0093）のつなぎ。格子計算は制御棒を入れた集合体を別の入力で解き、
制御棒の軸のないテーブルを出す。値ごとのテーブルは、公称の状態を制御棒を入れた状態にした別の組（ADR-0200）で、解いた点
（`points`）が値ごとに違ってよい。ここでは、群境界、燃焼度でない軸（名前、格子点、減速材の密度、B-10 の割合、内挿の変数、公称の
値）、部分ノードの数、初めの重金属の質量が同じで、燃焼度の点が最初のテーブル（棒のない状態）の燃焼度の点から選んだものである
ことを確かめ（違えばエラー。黙って結合しない）、値ごとの点に制御棒の軸の番号を足して（最後の軸）、点の辞書式の順に並べる。
燃焼度の軸は最初のテーブルのもので、ほかのテーブルの点の燃焼度の番号を付け直す。状態の項目はそのままコピーする（数値の
ビットを変えない）。

ミクロ燃焼の微視的断面積（`micro`。ADR-0216）は、どのテーブルも持つか、どれも持たないこと。持てば、追う核種が同じで、点ごとの
微視的断面積を状態と同じく並べ替える。燃焼度の点ごとの均質化した数密度は燃焼の履歴だけで決まり、制御棒の値によらないので、
最初のテーブルのものを使い、ほかのテーブルの同じ燃焼度の値が丸めの範囲で同じことを確かめる（ADR-0217）。

使い方: python3 rod_table.py <出力> <値>=<テーブル> ...（最初の値は棒のない状態 out）
"""

import sys

import msgpack


def read(path):
    with open(path, "rb") as f:
        table = msgpack.unpackb(f.read(), raw=False, strict_map_key=False)
    if table.get("format") != "xlc-table":
        sys.exit(f"{path}: not an xlc-table")
    if any(axis["name"] == "rod" for axis in table["axes"]):
        sys.exit(f"{path}: already has a rod axis")
    if "points" not in table:
        sys.exit(f"{path}: has no points")
    return table


def burnup_axis(table):
    """燃焼度の軸の番号（なければ None）"""
    names = [axis["name"] for axis in table["axes"]]
    return names.index("burnup") if "burnup" in names else None


def merge(pairs):
    """`pairs` は (値, パス, テーブル) の列。最初が棒のない状態"""
    first = pairs[0][2]
    b = burnup_axis(first)
    others = [a for i, a in enumerate(first["axes"]) if i != b]
    burnups = first["axes"][b]["points"] if b is not None else None
    sub = first["split"] ** 2
    first_micro = first.get("micro")
    rows = []
    for r, (value, path, table) in enumerate(pairs):
        for key in ("bounds", "split"):
            if table[key] != first[key]:
                sys.exit(
                    f"{path} ({value}): {key} differs from {pairs[0][1]}: {table[key]} vs {first[key]}"
                )
        # 初めの重金属の質量（履歴から作ったテーブルだけが持つ。ADR-0190）は、どのテーブルも同じ値か、どれも持たない
        if table.get("heavy_metal") != first.get("heavy_metal"):
            sys.exit(
                f"{path} ({value}): heavy_metal differs from {pairs[0][1]}: "
                f"{table.get('heavy_metal')} vs {first.get('heavy_metal')}"
            )
        tb = burnup_axis(table)
        if (tb is None) != (b is None):
            sys.exit(f"{path} ({value}): the burnup axis differs from {pairs[0][1]}")
        axes = [a for i, a in enumerate(table["axes"]) if i != tb]
        if axes != others or (tb is not None and tb != b):
            sys.exit(
                f"{path} ({value}): axes differ from {pairs[0][1]}: {table['axes']} vs {first['axes']}"
            )
        own = table["axes"][tb]["points"] if tb is not None else None
        if own is not None and any(p not in burnups for p in own):
            sys.exit(
                f"{path} ({value}): the burnups {own} are not points of {pairs[0][1]} {burnups}"
            )
        if len(table["states"]) != len(table["points"]) * sub:
            sys.exit(
                f"{path} ({value}): {len(table['states'])} states for {len(table['points'])} points"
            )
        micro = table.get("micro")
        check_micro(path, value, micro, first_micro, own, burnups, table["points"])
        for n, point in enumerate(table["points"]):
            point = list(point)
            if tb is not None:
                point[tb] = burnups.index(own[point[tb]])
            rows.append(
                (
                    point + [r],
                    table["states"][n * sub : (n + 1) * sub],
                    micro["states"][n] if micro is not None else None,
                )
            )
    rows.sort(key=lambda row: row[0])
    merged = dict(first)
    merged["axes"] = first["axes"] + [
        {"name": "rod", "values": [v for v, _, _ in pairs]}
    ]
    merged["points"] = [point for point, _, _ in rows]
    merged["states"] = [s for _, states, _ in rows for s in states]
    if first_micro is not None:
        merged["micro"] = dict(first_micro)
        merged["micro"]["states"] = [m for _, _, m in rows]
    return merged


# 同じ履歴の均質化した数密度の、制御棒の値ごとの格子の違いで許す相対差（格子の体積の和の丸め。数個の ulp の大きさ）
DENSITY_ROUNDING = 1e-12


def check_micro(path, value, micro, first_micro, own, burnups, points):
    """`micro`（燃焼度の点 `own`、点 `points` のテーブルのもの）が、最初のテーブルの `first_micro`（燃焼度の点 `burnups`）と合うこと"""
    if (micro is None) != (first_micro is None):
        sys.exit(f"{path} ({value}): micro data must be in every table or in none")
    if micro is None:
        return
    if micro["nuclides"] != first_micro["nuclides"]:
        sys.exit(
            f"{path} ({value}): the micro nuclides {micro['nuclides']} differ from {first_micro['nuclides']}"
        )
    if own is None or len(micro["densities"]) != len(own):
        sys.exit(f"{path} ({value}): micro densities need one row per burnup point")
    if len(micro["states"]) != len(points):
        sys.exit(
            f"{path} ({value}): {len(micro['states'])} micro states for {len(points)} points"
        )
    for burnup, row in zip(own, micro["densities"]):
        expected = first_micro["densities"][burnups.index(burnup)]
        # 比較を否定して、NaN も合わないものとして止める
        if len(row) != len(expected) or not all(
            abs(a - b) <= DENSITY_ROUNDING * max(abs(a), abs(b))
            for a, b in zip(row, expected)
        ):
            sys.exit(
                f"{path} ({value}): the micro densities at burnup {burnup} differ from the first table"
            )


def main():
    if len(sys.argv) < 4:
        sys.exit("usage: rod_table.py <output> <value>=<table> <value>=<table> ...")
    output = sys.argv[1]
    pairs = [arg.split("=", 1) for arg in sys.argv[2:]]
    if any(len(p) != 2 for p in pairs):
        sys.exit("each table is <value>=<path>")
    values = [v for v, _ in pairs]
    if values[0] != "out" or len(set(values)) != len(values):
        sys.exit(f"the values must start with out and be distinct, got {values}")
    merged = merge([(value, path, read(path)) for value, path in pairs])
    with open(output, "wb") as f:
        f.write(msgpack.packb(merged, use_bin_type=True))


if __name__ == "__main__":
    main()
