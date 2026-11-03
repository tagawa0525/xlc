"""1 つの入力の制御棒の分岐のテーブルと、値ごとの格子計算を rod_table.py でまとめたテーブルを比べる（#0093、ADR-0224）。

読んだ値の厳密な一致（`==`）と、違う項目ごとの最大の相対差、状態ごとの k∞ を出す。

使い方（このディレクトリで）:
    for v in out aic b4c; do ../../../bin/lattice assembly-$v.toml > lattice-$v.out; done
    python3 ../rod_table.py table-merged.msgpack out=fg-out.msgpack aic=fg-aic.msgpack b4c=fg-b4c.msgpack
    ../../../bin/lattice assembly.toml > lattice.out
    python3 compare_rod_branch.py table-branch.msgpack table-merged.msgpack
"""

import sys

import msgpack


def read(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def leaves(x, path=""):
    if isinstance(x, dict):
        for k in sorted(x, key=str):
            yield from leaves(x[k], f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, f"{path}[{i}]")
    else:
        yield path, x


a, b = read(sys.argv[1]), read(sys.argv[2])
print("equal:", a == b)
la, lb = dict(leaves(a)), dict(leaves(b))
print("keys only in branch:", sorted(set(la) - set(lb))[:20])
print("keys only in merged:", sorted(set(lb) - set(la))[:20])
worst = {}
for k in set(la) & set(lb):
    x, y = la[k], lb[k]
    if x == y:
        continue
    if isinstance(x, float) and isinstance(y, float):
        d = abs(x - y) / max(abs(x), abs(y))
    else:
        d = float("inf")
    field = k.split("/")[-1].split("[")[0] or k
    if d > worst.get(field, (0, ""))[0]:
        worst[field] = (d, k, x, y)
for field, (d, k, x, y) in sorted(worst.items()):
    print(f"{field}: max relative difference {d:.3e} at {k}: {x!r} vs {y!r}")
for n, s in enumerate(a["states"]):
    print(
        f"state {n} point {a['points'][n]}: k_inf {s['k_inf']!r} vs {b['states'][n]['k_inf']!r}"
    )
