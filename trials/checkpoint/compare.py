"""ADR-0196 の試行: 再開した履歴と一気に流した履歴、本番の履歴の最初の点を比べる。

使い方: python3 compare.py <再開した履歴> <一気に流した履歴> <本番の履歴>
"""

import hashlib
import sys

import msgpack


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def without_times(history):
    h = dict(history)
    h.pop("seconds")
    return msgpack.packb(h, use_bin_type=True)


resumed, straight, production = (read(p) for p in sys.argv[1:4])
for name, path in zip(["resumed", "straight", "production"], sys.argv[1:4]):
    with open(path, "rb") as f:
        print(name, "md5", hashlib.md5(f.read()).hexdigest())
print(
    "resumed == straight (bytes without seconds):",
    without_times(resumed) == without_times(straight),
)
print("solves", resumed["solves"], straight["solves"])
print("seconds resumed", resumed["seconds"])
print("seconds straight", straight["seconds"])
n = len(straight["points"])
same = all(
    a["k_inf"] == b["k_inf"]
    and a["densities"] == b["densities"]
    and a["burnup"] == b["burnup"]
    for a, b in zip(straight["points"], production["points"][:n])
)
print(
    f"straight vs production first {n} points (k_inf, densities, burnup bit-equal):",
    same,
)
print(
    "pin_powers equal:",
    all(
        a["pin_powers"] == b["pin_powers"]
        for a, b in zip(straight["points"], production["points"][:n])
    ),
)
print(
    "nuclides/regions/transportable equal:",
    all(straight[k] == production[k] for k in ["nuclides", "regions", "transportable"]),
)
