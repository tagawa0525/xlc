"""少数群定数のテーブルの不連続因子（辺の ADF と角の CDF）をすべて 1 にしたコピーを書く。2 段階の手法から不連続因子を
外したときの誤差を見る（#0037 の分け方の 3）。

使い方: python3 without_df.py <入力の fg.msgpack> <出力の fg.msgpack>
"""

import sys
from pathlib import Path

import msgpack

table = msgpack.unpackb(Path(sys.argv[1]).read_bytes(), raw=False)
for state in table["states"]:
    state["adf"] = [[1.0 for _ in side] for side in state["adf"]]
    if state["pins"] is not None:
        state["pins"]["corners"] = [
            [1.0 for _ in corner] for corner in state["pins"]["corners"]
        ]
Path(sys.argv[2]).write_bytes(msgpack.packb(table))
