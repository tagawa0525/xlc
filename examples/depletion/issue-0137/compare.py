"""#0137 の比較: 短い燃焼の履歴（xlc の msgpack）を、長い履歴の先頭の同じ数の点と比べる。

Usage: python3 compare.py <短い履歴> <長い履歴>

核種、燃焼の領域、輸送できる核種の一致と、各点（時間、燃焼度、k∞、数密度、ピン出力、少数群定数）がビットで同じか（msgpack に詰め直したバイト列の一致）を出し、
違えば k∞ と数密度の最大の相対差を出す。実行時間（`seconds`）と輸送の回数（`solves`）は比べない。
"""

import sys

import msgpack


def load(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def main():
    short, long = load(sys.argv[1]), load(sys.argv[2])
    n = len(short["points"])
    print(
        f"points {n} of {len(long['points'])}, solves {short['solves']} and {long['solves']}"
    )
    for key in ["format", "nuclides", "regions", "transportable"]:
        print(key, "equal" if short[key] == long[key] else "DIFFERS")
    ours, theirs = short["points"], long["points"][:n]
    # 値の == は -0.0 と 0.0、整数と浮動小数を等しいとするので、ビットの一致は詰め直したバイト列で比べる
    print("points bit-identical:", msgpack.packb(ours) == msgpack.packb(theirs))
    dk = dn = 0.0
    for a, b in zip(ours, theirs):
        dk = max(dk, abs(a["k_inf"] - b["k_inf"]) / abs(b["k_inf"]))
        for ra, rb in zip(a["densities"], b["densities"]):
            for x, y in zip(ra, rb):
                if x != y:
                    dn = max(dn, abs(x - y) / max(abs(x), abs(y)))
    print("max relative difference: k_inf", dk, "densities", dn)
    print("seconds", short["seconds"])


if __name__ == "__main__":
    main()
