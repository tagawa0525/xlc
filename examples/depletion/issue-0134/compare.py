"""#0134 の比較: 自己遮蔽をステップの始まりだけで解いた燃焼の履歴（案 F2）を、輸送のたびに解いた履歴（案 F1）の同じ点と比べる。

Usage: python3 compare.py <F2 の履歴> <F1 の履歴>

燃焼度、核種、燃焼の領域が同じことを確かめ、F2 の点の数だけ（F1 が長ければ先頭）、点ごとに k∞ の差 (k_F2 − k_F1) × 1e5 pcm、
燃焼の領域の物質の種類（燃料 `fuel`、Gd 燃料のリング `gd*`）ごとの主要な核種の数密度の相対差（種類の領域の体積の重みの平均の差と、
領域ごとの差の絶対値の最大）、ピン出力の相対差の絶対値の最大（F1 のピン出力が正のセル）を出す。F1 が 0 の量で F2 が 0 でなければ止める。最後に、それぞれの絶対値の最大と、
時間の内訳（`seconds`）、輸送の回数（`solves`）、自己遮蔽の回数（`shieldings`。F2 だけ。F1 は輸送の回数に等しい）を出す。F2 が輸送ごとに自己遮蔽を解いていれば止める。
"""

import sys

import msgpack

NUCLIDES = {
    "fuel": ["U-235", "Pu-239", "Xe-135", "Sm-149"],
    "gd": ["U-235", "Pu-239", "Xe-135", "Sm-149", "Gd-155", "Gd-157"],
}


def load(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def kind(material):
    return "gd" if material.startswith("gd") else material


def main():
    f2, f1 = load(sys.argv[1]), load(sys.argv[2])
    n = len(f2["points"])
    if n > len(f1["points"]):
        raise SystemExit(f"{n} points of F2 for {len(f1['points'])} of F1")
    for key in ["nuclides", "regions"]:
        if f2[key] != f1[key]:
            raise SystemExit(f"{key} differ")
    a, b = f2["points"], f1["points"][:n]
    if [p["burnup"] for p in a] != [p["burnup"] for p in b]:
        raise SystemExit("burnups differ")
    index = {name: i for i, name in enumerate(f2["nuclides"])}
    regions = f2["regions"]
    kinds = {}
    for r, region in enumerate(regions):
        kinds.setdefault(kind(region["material"]), []).append(r)
    worst = {"k": 0.0, "pins": 0.0}
    print(
        "burnup | k F1 | dk (pcm) | "
        + " | ".join(f"{k} {name} mean/max (%)" for k in kinds for name in NUCLIDES[k])
        + " | pins max (%)"
    )
    for p2, p1 in zip(a, b):
        dk = (p2["k_inf"] - p1["k_inf"]) * 1e5
        worst["k"] = max(worst["k"], abs(dk))
        cells = []
        for k, members in kinds.items():
            total = sum(regions[r]["volume"] for r in members)
            for name in NUCLIDES[k]:
                i = index[name]
                m2 = (
                    sum(p2["densities"][r][i] * regions[r]["volume"] for r in members)
                    / total
                )
                m1 = (
                    sum(p1["densities"][r][i] * regions[r]["volume"] for r in members)
                    / total
                )
                # F1 が 0 の量は相対差を取れないので、F2 も 0 であることだけを確かめる（初めに 0 の Xe-135 など）
                if m1 <= 0 and m2 != m1:
                    raise SystemExit(
                        f"{p2['burnup']} MWd/kg {k} {name}: F1 {m1}, F2 {m2}"
                    )
                mean = (m2 - m1) / m1 if m1 > 0 else 0.0
                largest = 0.0
                for r in members:
                    x, y = p2["densities"][r][i], p1["densities"][r][i]
                    if y > 0:
                        largest = max(largest, abs(x - y) / y)
                    elif x != y:
                        raise SystemExit(
                            f"{p2['burnup']} MWd/kg region {r} {name}: F1 {y}, F2 {x}"
                        )
                key = (k, name)
                worst[key] = (
                    max(worst.get(key, (0.0, 0.0))[0], abs(mean)),
                    max(worst.get(key, (0.0, 0.0))[1], largest),
                )
                cells.append(f"{mean * 100:+.4f}/{largest * 100:.4f}")
        pins = 0.0
        for r2, r1 in zip(p2["pin_powers"], p1["pin_powers"]):
            for x, y in zip(r2, r1):
                if y > 0:
                    pins = max(pins, abs(x / y - 1))
                elif x != y:
                    raise SystemExit(f"{p2['burnup']} MWd/kg pin power: F1 {y}, F2 {x}")
        worst["pins"] = max(worst["pins"], pins)
        print(
            f"{p2['burnup']:g} | {p1['k_inf']:.6f} | {dk:+.2f} | "
            + " | ".join(cells)
            + f" | {pins * 100:.4f}"
        )
    print(f"max |dk| {worst['k']:.2f} pcm, max pins {worst['pins'] * 100:.4f} %")
    for key, (mean, largest) in (
        (k, v) for k, v in worst.items() if isinstance(k, tuple)
    ):
        print(
            f"max {key[0]} {key[1]}: mean {mean * 100:.4f} %, region {largest * 100:.4f} %"
        )
    # F2 の履歴は自己遮蔽の回数を持つ。比べる本番の F1 の履歴は回数を書く前の形式で、F1 は輸送ごとに解くので回数は輸送の回数
    if f2["shieldings"] >= f2["solves"]:
        raise SystemExit(
            f"F2 solved the shielding {f2['shieldings']} times in {f2['solves']} transports"
        )
    if "shieldings" in f1 and f1["shieldings"] != f1["solves"]:
        raise SystemExit(
            f"F1 solved the shielding {f1['shieldings']} times in {f1['solves']} transports"
        )
    print(
        "F2 solves",
        f2["solves"],
        "shieldings",
        f2["shieldings"],
        "seconds",
        f2["seconds"],
    )
    print("F1 solves", f1["solves"], "seconds", f1["seconds"])


if __name__ == "__main__":
    main()
