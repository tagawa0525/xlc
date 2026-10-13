"""試行: 1 状態だけ解いたテーブルの状態が、8 状態のテーブルの最初の状態（ホウ素、減速材、燃料が最小の点）とビット単位で
同じかを確かめる（実行ファイルを作り直したあとに、前の出力を使い続けてよいかの判断）。

使い方: python3 same_state.py <1 状態のテーブル> <8 状態のテーブル>...（対で並べる）
"""

import sys

import msgpack


def load(path: str) -> dict:
    return msgpack.unpackb(open(path, "rb").read(), raw=False, strict_map_key=False)


def main() -> None:
    args = sys.argv[1:]
    if len(args) % 2:
        raise SystemExit(__doc__)
    bad = 0
    for one, eight in zip(args[::2], args[1::2]):
        a, b = load(one), load(eight)
        if len(a["states"]) != 1 or len(b["states"]) != 8:
            raise SystemExit(f"{one}, {eight}: expected 1 and 8 states")
        same = a["states"][0] == b["states"][0] and a["bounds"] == b["bounds"]
        bad += not same
        print(f"{'same' if same else 'DIFFERENT'}: {one} vs state 0 of {eight}")
    if bad:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
