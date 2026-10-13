"""制御棒の値ごとの格子計算のテーブルを、制御棒の離散の軸（ADR-0097、ADR-0127）を持つ 1 つのテーブルにまとめる。

格子計算が制御棒の分岐を作れるようになるまで（#0093）のつなぎ。格子計算は制御棒を入れた集合体を別の入力で解き、
軸のない（連続の軸だけの）テーブルを出す。ここでは、群境界、連続の軸（名前、格子点、減速材の密度、B-10 の割合）、
部分ノードの数、状態の数が同じであることを確かめ（違えばエラー。黙って結合しない）、値ごとの状態を制御棒の軸
（最後の軸。いちばん速く変わる）に並べる。状態の項目はそのまま写す（数値のビットを変えない）。

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
    return table


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
    tables = [read(path) for _, path in pairs]
    first = tables[0]
    for (value, path), table in zip(pairs, tables):
        for key in ("bounds", "axes", "split"):
            if table[key] != first[key]:
                sys.exit(
                    f"{path} ({value}): {key} differs from {pairs[0][1]}: {table[key]} vs {first[key]}"
                )
        if len(table["states"]) != len(first["states"]):
            sys.exit(
                f"{path} ({value}): {len(table['states'])} states, {pairs[0][1]} has {len(first['states'])}"
            )
    sub = first["split"] ** 2
    points = len(first["states"]) // sub
    states = []
    for i in range(points):
        for table in tables:
            states.extend(table["states"][i * sub : (i + 1) * sub])
    merged = dict(first)
    merged["axes"] = first["axes"] + [{"name": "rod", "values": values}]
    merged["states"] = states
    with open(output, "wb") as f:
        f.write(msgpack.packb(merged, use_bin_type=True))


if __name__ == "__main__":
    main()
