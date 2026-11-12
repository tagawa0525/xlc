"""テーブルの状態の k∞ と、同じ燃焼度の格子計算の燃焼の履歴の k∞ の反応度の差を書く（#0283）。

使い方: python3 table_history_k.py <格子計算の履歴（msgpack）> <テーブル（msgpack）>...

テーブルごとに、状態の燃焼度、テーブルの k∞、履歴の k∞、(1/k_履歴 − 1/k_テーブル) × 1e5 pcm を 1 行ずつ書く。履歴は 1 つの状態の
燃焼なので、燃焼度のほかの軸が 1 点だけ（燃焼度ごとに状態が 1 つ）で、集合体を分けないテーブルだけを受け付ける。テーブルの燃焼度の
点が履歴の燃焼点に無ければ止める。
"""

import sys
from pathlib import Path

import msgpack


def load(path: str) -> dict:
    return msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)


def rows(history: dict, table: dict) -> list[tuple[float, float, float, float]]:
    """テーブルの状態ごとの (燃焼度, k_テーブル, k_履歴, pcm)"""
    if table["split"] != 1:
        raise SystemExit("the table must not split the assembly")
    names = [axis["name"] for axis in table["axes"]]
    for axis in table["axes"]:
        if axis["name"] != "burnup" and len(axis["points"]) != 1:
            raise SystemExit(
                f"the table has {len(axis['points'])} {axis['name']} points; "
                "the history holds one state per burnup"
            )
    # 履歴の燃焼点は燃焼度の値で引く（テーブルの燃焼度の軸は履歴の点の値をそのまま持つ）
    history_k = {point["burnup"]: point["k_inf"] for point in history["points"]}
    axis = names.index("burnup")
    burnups = table["axes"][axis]["points"]
    out = []
    for point, state in zip(table["points"], table["states"], strict=True):
        burnup = burnups[point[axis]]
        if burnup not in history_k:
            raise SystemExit(f"burnup {burnup} is not a history point")
        k_history = history_k[burnup]
        k_table = state["k_inf"]
        out.append(
            (burnup, k_table, k_history, (1.0 / k_history - 1.0 / k_table) * 1e5)
        )
    return out


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    history = load(sys.argv[1])
    for path in sys.argv[2:]:
        print(path)
        for burnup, k_table, k_history, pcm in rows(history, load(path)):
            print(f"  {burnup:g} {k_table:.9f} {k_history:.9f} {pcm:+.3f}")


if __name__ == "__main__":
    main()
