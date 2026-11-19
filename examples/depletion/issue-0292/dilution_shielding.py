"""足す核種の、燃料の組成の実際の希釈での自己遮蔽の目安を、輸送計算なしで出す（#0292）。

使い方: python3 dilution_shielding.py --json <出力> <履歴（msgpack）> <テーブル（msgpack）> <吸収の割合の JSON>
          <輸送する核種の xs のファイルの置き場> <足す核種の xs のファイルの置き場> <足した核種の一覧> [--temperature 565]
        python3 dilution_shielding.py --report <出力の JSON>

核種ごとの xs のファイル（`xs nuclide` の出力。ADR-0138）の Bondarenko の表（`shielded`。σ0 ごとの群の捕獲と核分裂）を、燃焼の領域の
組成の実際の希釈 σ0,g = Σ_{j≠i} Nⱼ σt,j,g(∞) / Nᵢ（燃焼の領域の体積で平均した数密度。ほかの核種の全断面積は無限希釈）で、ln σ0 の
線形で内挿する。燃料の外への逃げの項（等価原理の Σe）は足さないので、σ0 は低めで、自己遮蔽は大きめに出る（保守的な側）。

1 群にまとめる中性子束の重みは、輸送計算の細かい群の中性子束が無いので、0.625 eV より上は単位のレサジーあたり一定（1/E）、下は燃料の
温度の Maxwell 分布とし、2 つの範囲の比はテーブルの同じ燃焼度の基準条件の状態の 2 群の無限媒質の中性子束（untransported_absorption.py の
`two_group_flux`）にする。比 r = Σ σₐ,g(σ0) w_g / Σ σₐ,g(∞) w_g。目安の差は 1e5 fᵢ cᵢ (1 − r) / k∞（compare_shielding.py と同じ）。
重みの確かさは、無限希釈の 1 群の吸収断面積を同じ重みで作った値の、テーブルの冷却計算の値（MOC の中性子束で重みを付けた値）に対する
比（`weight_check`）で見る。
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import msgpack
from untransported_absorption import (
    absorption,
    nominal_states,
    region_sum,
    two_group_flux,
)

# 熱群と高速群の境界 (eV)。テーブルの 2 群の境界と同じ
THERMAL = 0.625
BOLTZMANN = 8.617333262e-5


def read(path):
    return msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)


def temperature_data(nuclide_file, temperature):
    data = nuclide_file["nuclides"][0]
    matched = [t for t in data["temperatures"] if t["temperature"] == temperature]
    if len(matched) != 1:
        raise SystemExit(f"{data['name']} has no data at {temperature} K")
    return matched[0]


def weights(bounds, temperature, fast, thermal):
    """群ごとの中性子束の重み（`bounds` は降順の群の境界 (eV)。熱外と熱の和はそれぞれ `fast` と `thermal`）"""
    kt = BOLTZMANN * temperature
    out = []
    for hi, lo in itertools.pairwise(bounds):
        if lo >= THERMAL:
            out.append(("epithermal", math.log(hi / lo)))
        else:
            top = min(hi, THERMAL)

            def cdf(e):
                x = e / kt
                return 1.0 - (1.0 + x) * math.exp(-x)

            out.append(("thermal", cdf(top) - cdf(lo)))
    sums = {"epithermal": 0.0, "thermal": 0.0}
    for kind, w in out:
        sums[kind] += w
    scale = {"epithermal": fast, "thermal": thermal}
    return [scale[kind] * w / sums[kind] for kind, w in out]


def interpolate(sigma0, table, value):
    """ln σ0 の線形の内挿（`sigma0` は降順。範囲の外は端の値）"""
    if value >= sigma0[0]:
        return table[0]
    for (s0, t0), (s1, t1) in zip(zip(sigma0, table), zip(sigma0[1:], table[1:])):
        if value >= s1:
            x = (math.log(value) - math.log(s1)) / (math.log(s0) - math.log(s1))
            return t1 + (t0 - t1) * x
    return table[-1]


def share(point, name, pruned):
    """吸収の割合の JSON の点の、核種の割合と数密度の比。JSON から外した核種（最大の目安が床以下。untransported_absorption.py の
    PRUNE_PCM）は 0。どちらでもない名前はエラー（足す核種の一覧と JSON の燃料の取り違え）"""
    if name in point["untransported_all"]:
        return point["untransported_all"][name]
    if name in pruned:
        return {"fraction": 0.0, "captured_ratio": 1.0}
    raise SystemExit(f"{name} is not an untransported nuclide of the absorption JSON")


def check_grid(reference, other, name):
    """xs のファイルの σ0 の点と群の境界が、最初のファイルと同じこと（違えば表の行と σ0 の対応がずれる）"""
    if (
        other["sigma0"] != reference["sigma0"]
        or other["group_bounds"] != reference["group_bounds"]
    ):
        raise SystemExit(
            f"{name}: the xs file has another sigma0 grid or group structure"
        )


def directory_md5(directory, names):
    """置き場の核種の xs のファイルの md5 を名前の順につないだものの md5"""
    digest = hashlib.md5()
    for name in sorted(names):
        digest.update(
            hashlib.md5((Path(directory) / f"{name}.msgpack").read_bytes()).digest()
        )
    return digest.hexdigest()


def absorption_groups(shielded_entry):
    return [c + f for c, f in zip(shielded_entry["capture"], shielded_entry["fission"])]


def ratio(added, background, density, sigma0, w, temperature):
    """足す核種の、組成の希釈の 1 群の吸収断面積と無限希釈の比、と無限希釈の 1 群の値（重み `w` の和で割る）"""
    data = temperature_data(added, temperature)
    tables = [absorption_groups(s) for s in data["shielded"]]
    groups = len(w)
    shielded = []
    for g in range(groups):
        value = background[g] / density if density > 0.0 else math.inf
        shielded.append(interpolate(sigma0, [t[g] for t in tables], value))
    infinite = tables[0]
    one = sum(a * x for a, x in zip(infinite, w))
    return sum(a * x for a, x in zip(shielded, w)) / one, one / sum(w)


def summarize(history, table, shares, transported_dir, added_dir, names, temperature):
    index = {n: i for i, n in enumerate(history["nuclides"])}
    volume = sum(r["volume"] for r in history["regions"])
    points = {p["burnup"]: p for p in history["points"]}
    cache = {}
    first = read(Path(added_dir) / f"{names[0]}.msgpack")

    def nuclide(directory, name):
        key = (directory, name)
        if key not in cache:
            cache[key] = read(Path(directory) / f"{name}.msgpack")
            check_grid(first, cache[key], name)
        return cache[key]

    burnups, states = nominal_states(table)
    cooling = table["cooling"]
    pruned = set(shares.get("pruned", {}).get("nuclides", []))
    bounds, sigma0 = first["group_bounds"], first["sigma0"]
    rows = {n: [] for n in names}
    for p in shares["points"]:
        b = p["burnup"]
        hp = points.get(b)
        if hp is None:
            raise SystemExit(f"burnup {b} is not a point of the history")

        def mean(name, hp=hp):
            return region_sum(history, hp, index[name]) / volume

        k = burnups.index(b)
        fast, thermal = two_group_flux(table["states"][states[k]])
        w = weights(bounds, temperature, fast, thermal)
        total = [0.0] * (len(bounds) - 1)
        for name in history["transportable"]:
            n = mean(name)
            if n <= 0.0:
                continue
            sigma = temperature_data(nuclide(transported_dir, name), temperature)[
                "shielded"
            ][0]["total"]
            for g, s in enumerate(sigma):
                total[g] += n * s
        for name in names:
            n = mean(name)
            r, one = ratio(nuclide(added_dir, name), total, n, sigma0, w, temperature)
            cooled = absorption(cooling["points"][k][cooling["nuclides"].index(name)])
            entry = share(p, name, pruned)
            rows[name].append(
                {
                    "burnup": b,
                    "density": n,
                    "ratio": r,
                    "infinite_one_group": one,
                    "weight_check": one / cooled,
                    "shielding_pcm": 1e5
                    * entry["fraction"]
                    * entry["captured_ratio"]
                    * (1.0 - r)
                    / p["k_inf"],
                }
            )
    out = []
    for name, pts in rows.items():
        largest = max(pts, key=lambda x: abs(x["shielding_pcm"]))
        out.append(
            {
                "nuclide": name,
                "max_pcm": largest["shielding_pcm"],
                "burnup": largest["burnup"],
                "points": pts,
            }
        )
    return sorted(out, key=lambda r: -abs(r["max_pcm"]))


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        "核種 | 自己遮蔽による目安の差の最大 (pcm) | 燃焼度 | その点の σₐ の比 | 重みの確かさ（1 群の無限希釈 / 冷却計算）"
    )
    for r in summary["nuclides"]:
        at = next(p for p in r["points"] if p["burnup"] == r["burnup"])
        lines.append(
            f"{r['nuclide']} | {r['max_pcm']:+.2f} | {r['burnup']:g} | {at['ratio']:.4f} | {at['weight_check']:.3f}"
        )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument("--temperature", type=float, default=565.0)
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not args.json or len(args.files) != 6:
        raise SystemExit(__doc__)
    history_path, table_path, shares_path, transported_dir, added_dir, list_path = (
        args.files
    )
    names = Path(list_path).read_text().split()
    summary = {
        "format": "xlc-dilution-shielding",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in (history_path, table_path, shares_path, list_path)
        ],
        # 比を決める xs のファイルの置き場ごとの md5（足す核種と、履歴の輸送する核種）
        "xs_files": {
            "added": directory_md5(added_dir, names),
            "transported": directory_md5(
                transported_dir, read(history_path)["transportable"]
            ),
        },
        "temperature": args.temperature,
        "nuclides": summarize(
            read(history_path),
            read(table_path),
            json.loads(Path(shares_path).read_text()),
            transported_dir,
            added_dir,
            names,
            args.temperature,
        ),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
