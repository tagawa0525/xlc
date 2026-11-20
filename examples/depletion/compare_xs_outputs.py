"""2 つの `xs nuclide` の出力（msgpack）を、核種と温度ごとに値で比べる（#0325）。

使い方: python3 compare_xs_outputs.py <出力 a> <出力 b>

ビルドの違う 2 つの xs で同じ入力を処理した出力を比べ、変わった量と、その大きさと群を 1 行ずつ書く。核種は名前で、温度は並びで
対応させ、形（核種の名前、温度、群の数、散乱の次数、σ0 の数）が違えば止まる。

- `chi`: χ の差の絶対値の最大と、その群（0 始まり）、差の絶対値の和
- `scatter P<l>`: 散乱行列の l 次のモーメントの差の絶対値を、a の P0 のその行（入射の群）の和で割った値の最大と、その（入射の群、
  出射の群）
- `shielded <量>`: 自己遮蔽した群定数（全断面積、弾性散乱、核分裂、νσf、捕獲）の σ0 と群ごとの相対差 |a − b| / max(|a|, |b|)
  の最大と、その（σ0 の番号、群）
- 差が無い量は書かず、核種と温度のどの量にも差が無ければ `same` と書く
"""

import sys

import msgpack

SHIELDED = ["total", "elastic", "fission", "nu_fission", "capture"]


def load(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read())


def relative(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def compare_temperature(a, b):
    """1 つの温度の差の行（文字列）の一覧"""
    lines = []
    if len(a["chi"]) != len(b["chi"]):
        raise ValueError("chi lengths differ")
    diffs = [abs(x - y) for x, y in zip(a["chi"], b["chi"])]
    if any(diffs):
        g = max(range(len(diffs)), key=diffs.__getitem__)
        lines.append(
            f"chi: max |diff| {diffs[g]:.3e} in group {g}, sum {sum(diffs):.3e}"
        )
    if len(a["scatter"]) != len(b["scatter"]):
        raise ValueError("scattering orders differ")
    rows = [sum(row) for row in a["scatter"][0]]
    for order, (ma, mb) in enumerate(zip(a["scatter"], b["scatter"])):
        worst = (0.0, 0, 0)
        for g, (ra, rb) in enumerate(zip(ma, mb)):
            if len(ra) != len(rb):
                raise ValueError("scattering rows differ")
            for h, (x, y) in enumerate(zip(ra, rb)):
                d = abs(x - y) / rows[g] if rows[g] > 0 else abs(x - y)
                if d > worst[0]:
                    worst = (d, g, h)
        if worst[0] > 0:
            lines.append(
                f"scatter P{order}: max |diff| / P0 row sum {worst[0]:.3e} in {worst[1]} -> {worst[2]}"
            )
    if len(a["shielded"]) != len(b["shielded"]):
        raise ValueError("sigma0 counts differ")
    for name in SHIELDED:
        worst = (0.0, 0, 0)
        for s, (sa, sb) in enumerate(zip(a["shielded"], b["shielded"])):
            for g, (x, y) in enumerate(zip(sa[name], sb[name])):
                d = relative(x, y)
                if d > worst[0]:
                    worst = (d, s, g)
        if worst[0] > 0:
            lines.append(
                f"shielded {name}: max rel diff {worst[0]:.3e} at sigma0 {worst[1]}, group {worst[2]}"
            )
    return lines


def compare(a, b):
    """核種と温度ごとの差の行の一覧"""
    by_name = {n["name"]: n for n in b["nuclides"]}
    if sorted(by_name) != sorted(n["name"] for n in a["nuclides"]):
        raise ValueError("nuclides differ")
    lines = []
    for na in a["nuclides"]:
        nb = by_name[na["name"]]
        if [t["temperature"] for t in na["temperatures"]] != [
            t["temperature"] for t in nb["temperatures"]
        ]:
            raise ValueError(f"{na['name']}: temperatures differ")
        for ta, tb in zip(na["temperatures"], nb["temperatures"]):
            prefix = f"{na['name']} {ta['temperature']:g} K"
            diffs = compare_temperature(ta, tb)
            lines.extend(f"{prefix} {d}" for d in diffs) if diffs else lines.append(
                f"{prefix} same"
            )
    return lines


def main(argv):
    if len(argv) != 3:
        sys.exit("usage: compare_xs_outputs.py <output a> <output b>")
    for line in compare(load(argv[1]), load(argv[2])):
        print(line)


if __name__ == "__main__":
    main(sys.argv)
