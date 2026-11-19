"""燃焼の領域の核種の中性子の吸収を、ミクロ燃焼で追う 31 核種、それ以外の輸送する核種、輸送しない核種（ライブラリの放射化の核種）に
分け、燃焼度の点ごとに格子全体の吸収に占める割合を出す（#0292）。

使い方: python3 untransported_absorption.py --json <出力> --chain <チェーン（msgpack）> <燃焼の履歴（msgpack）> <テーブル（msgpack）>
        python3 untransported_absorption.py --report <出力の JSON>

テーブルは燃焼の履歴からのブランチ計算（`[history]`）で、燃焼度の点ごとの冷却計算の 1 群の微視的断面積（`cooling`。ADR-0238）を持ち、燃焼度の
ほかの軸がどれも 1 点のもの。1 群の微視的断面積 σⱼ は、燃焼の領域の反応率の和を、格子の体積で均質化した数密度 N̄ⱼ と、体積で
積分した中性子束で割った値（crates/lattice/src/history.rs の `cooling`）なので、核種 j の吸収率は格子で均質化すると N̄ⱼ σₐⱼ φ̄
（φ̄ は全群の中性子束）。σₐ は、チェーンの反応の中性子を失う分: (n,γ)、(n,p)、(n,α) と核分裂の和から、(n,2n) を 1 倍、(n,3n) を 2 倍、
(n,4n) を 3 倍引く（出ていく中性子の数 − 1）。N̄ⱼ は、履歴の燃焼の領域の数密度を体積で足して格子の体積で割る。格子の体積は
テーブルの追う核種の均質化した数密度と同じになるように、追う核種の U-238 から求め、ほかの追う核種で相対 1e-9 に合うことを確かめる。

格子全体の吸収は、テーブルの状態の 2 群の吸収断面積と、同じ状態の 2 群の無限媒質の中性子束（compare_linearization.py の k_inf と
同じ式）で Σ_G Σₐ,G φ̄_G。中性子束の規格は比に入らない（σⱼ の 1 群の φ̄ は 2 群の和）。割合 f = N̄ⱼ σₐⱼ φ̄ / Σ_G Σₐ,G φ̄_G。
反応度の目安は、その吸収を除いたときの 1 次の変化 Δρ ≈ −f / k∞（k∞ = P / A の無限媒質で、吸収 δA を足すと Δρ = −δA / P）。
スペクトルと自己遮蔽の変化を含まない（輸送しない核種は格子計算の燃焼で崩壊だけで追うので、その数密度は捕獲で減らない。ADR-0238）。
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import msgpack

# ADR-0035 の元素の輸送する核種 26 と、その輸送する先行核 5（ADR-0216。inputs.py の TRACKED と同じ）
TRACKED = [
    "U-234", "U-235", "U-236", "U-238", "Np-237", "Np-239", "Pu-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242",
    "Am-241", "Am-243", "Xe-135", "Sm-147", "Sm-149", "Sm-150", "Sm-151", "Sm-152",
    "Gd-152", "Gd-154", "Gd-155", "Gd-156", "Gd-157", "Gd-158", "Gd-160",
    "Cm-242", "Cm-244", "Eu-153", "I-135", "Pm-148m",
]  # fmt: skip
# チェーンの反応の、出ていく中性子の数 − 1 の符号を変えたもの（吸収に数える重み）
LOSS = {102: 1.0, 103: 1.0, 107: 1.0, 16: -1.0, 17: -2.0, 37: -3.0}
VOLUME_TOLERANCE = 1e-9
# 一覧に出す核種の数
TOP = 15
# 順位の累積で、足す核種の数を数える残りの閾値 (pcm)
THRESHOLDS_PCM = (50.0, 5.0)
# xs が輸送する核種として処理しない核種（examples/depletion/p2/xs.toml の説明。#0127）。放射化の断面積は持つので一覧に出る
# JSON から外す核種の最大の目安の床 (pcm)。外した核種の数と目安の和は JSON に残す
PRUNE_PCM = 0.01
XS_EXCLUDED = ("Pm-147", "Pm-148", "Pm-149", "Am-242m", "Eu-151")


def read(path):
    return msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)


def absorption(nuclide):
    """1 群の微視的断面積の、中性子を失う分 (barn)"""
    total = nuclide["fission"]
    for r in nuclide["reactions"]:
        if r["mt"] not in LOSS:
            raise SystemExit(f"no neutron balance for MT {r['mt']}")
        total += LOSS[r["mt"]] * sum(r["channels"])
    return total


def solve2(matrix, vector):
    (a, b), (c, d) = matrix
    det = a * d - b * c
    return [
        (d * vector[0] - b * vector[1]) / det,
        (a * vector[1] - c * vector[0]) / det,
    ]


def two_group_flux(state):
    """(diag(Σt) − Sᵀ) φ = χ の解（k = νΣf · φ）"""
    total, scatter = state["total"], state["scatter"]
    if len(total) != 2:
        raise SystemExit("the table must have two groups")
    removal = [
        [(total[g] if g == h else 0.0) - scatter[h][g] for h in range(2)]
        for g in range(2)
    ]
    return solve2(removal, state["chi"])


def nominal_states(table):
    """燃焼度の点ごとの、ほかの軸がどれも基準点の状態の番号（燃焼度の軸は最初）"""
    if table.get("format") != "xlc-table" or "cooling" not in table:
        raise SystemExit("the table must be an xlc-table with cooling cross sections")
    burnup, *rest = table["axes"]
    if burnup["name"] != "burnup":
        raise SystemExit("the first axis must be the burnup axis")
    nominal = []
    for axis in rest:
        if len(axis["points"]) == 1:
            nominal.append(0)
        elif "nominal" in axis and axis["nominal"] in axis["points"]:
            nominal.append(axis["points"].index(axis["nominal"]))
        else:
            raise SystemExit(f"the {axis['name']} axis has no nominal point")
    by_burnup = {}
    for s, point in enumerate(table["points"]):
        if list(point[1:]) == nominal:
            by_burnup[point[0]] = s
    if sorted(by_burnup) != list(range(len(burnup["points"]))):
        raise SystemExit("the table lacks the nominal state at some burnup")
    return burnup["points"], [by_burnup[k] for k in range(len(burnup["points"]))]


def lattice_volume(history, table, point, burnup):
    """追う核種の均質化した数密度が履歴の燃焼の領域の体積の和 / 格子の体積になる格子の体積"""
    names = table["micro"]["nuclides"]
    densities = table["micro"]["densities"][point]
    index = {n: i for i, n in enumerate(history["nuclides"])}
    hp = history_point(history, burnup)

    def summed(name):
        return sum(
            r["volume"] * max(n[index[name]], 0.0)
            for r, n in zip(history["regions"], hp["densities"])
        )

    if "U-238" not in names:
        raise SystemExit("the table must track U-238 to fix the lattice volume")
    volume = summed("U-238") / densities[names.index("U-238")]
    for name, d in zip(names, densities):
        if d > 0.0 and abs(summed(name) / volume / d - 1.0) > VOLUME_TOLERANCE:
            raise SystemExit(
                f"{name} at {burnup} MWd/kg: the history and the table disagree in the homogenized density"
            )
    return volume


def history_point(history, burnup):
    matched = [p for p in history["points"] if p["burnup"] == burnup]
    if len(matched) != 1:
        raise SystemExit(f"burnup {burnup} is not a point of the history")
    return matched[0]


def category(name, transportable):
    if name in TRACKED:
        return "tracked"
    return "transported" if name in transportable else "untransported"


def region_sum(history, point, i):
    """燃焼の領域の数密度 × 体積の和（負の丸めは 0）"""
    return sum(
        r["volume"] * max(n[i], 0.0)
        for r, n in zip(history["regions"], point["densities"])
    )


def interpolate(xs, ys, x):
    """区分線形（範囲の外は端の値）"""
    if x <= xs[0]:
        return ys[0]
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return ys[-1]


def captured_ratios(history, burnups, sigmas, names, half_lives):
    """輸送しない核種が自分の吸収でも減るとしたときの数密度の、崩壊だけの数密度（履歴の値）への比を、テーブルの燃焼度の点ごとに。

    履歴の燃焼点の区間ごとに、崩壊だけの数密度 N_d の変化から一定の生成率 P を求め（dN_d/dt = P − λ N_d）、同じ P で
    dN/dt = P − (λ + σₐ φ̄) N を解く。σₐ φ̄ は 1 群の吸収断面積（区間の中点の燃焼度でテーブルの点の値を線形に内挿）と、U-235 の
    数密度の減り方から求めた φ̄ = −d ln N_U235 / dt / σₐ,U235（U-234 の捕獲による U-235 の生成を無視する）の積。ほかの核種の
    吸収で生まれる分（Eu-154 の捕獲で増える Eu-155 など）は足さないので、捕獲で生まれる核種の比は低めに出る
    """
    index = {n: i for i, n in enumerate(history["nuclides"])}
    u235 = index["U-235"]
    points = history["points"]
    phis = []
    for a, b in itertools.pairwise(points):
        dt = b["time"] - a["time"]
        middle = 0.5 * (a["burnup"] + b["burnup"])
        rate = (
            -math.log(region_sum(history, b, u235) / region_sum(history, a, u235)) / dt
        )
        phis.append((dt, middle, rate / interpolate(burnups, sigmas["U-235"], middle)))
    out = {}
    for name in names:
        i = index[name]
        half_life = half_lives.get(name)
        decay = math.log(2.0) / half_life if half_life else 0.0
        captured = region_sum(history, points[0], i)
        ratios = {}
        for (dt, middle, phi), a, b in zip(phis, points, points[1:]):
            before, after = region_sum(history, a, i), region_sum(history, b, i)
            if decay * dt > 1e-12:
                production = (
                    (after - before * math.exp(-decay * dt))
                    * decay
                    / -math.expm1(-decay * dt)
                )
            else:
                production = (after - before) / dt
            removal = decay + interpolate(burnups, sigmas[name], middle) * phi
            if removal * dt > 1e-12:
                captured = (
                    captured * math.exp(-removal * dt)
                    + production * -math.expm1(-removal * dt) / removal
                )
            else:
                captured += production * dt
            if b["burnup"] in burnups:
                ratios[b["burnup"]] = captured / after if after > 0.0 else 1.0
        out[name] = ratios
    return out


def summarize(history, table, half_lives):
    if history.get("format") != "xlc-depletion-history":
        raise SystemExit("the first file must be a lattice depletion history")
    burnups, states = nominal_states(table)
    cooling = table["cooling"]
    transportable = set(history["transportable"])
    index = {n: i for i, n in enumerate(history["nuclides"])}
    missing = [n for n in cooling["nuclides"] if n not in index]
    if missing:
        raise SystemExit(f"the history has no {missing[:5]}")
    untransported = [
        n for n in cooling["nuclides"] if category(n, transportable) == "untransported"
    ]
    sigmas = {
        name: [absorption(point[j]) for point in cooling["points"]]
        for j, name in enumerate(cooling["nuclides"])
        if name in untransported or name == "U-235"
    }
    ratios = captured_ratios(history, burnups, sigmas, untransported, half_lives)
    points = []
    for k, burnup in enumerate(burnups):
        if burnup == 0.0:
            continue
        state = table["states"][states[k]]
        flux = two_group_flux(state)
        lattice_absorption = sum(a * f for a, f in zip(state["absorption"], flux))
        volume = lattice_volume(history, table, k, burnup)
        hp = history_point(history, burnup)
        shares = {}
        for name, sigma in zip(cooling["nuclides"], cooling["points"][k]):
            density = (
                sum(
                    r["volume"] * max(n[index[name]], 0.0)
                    for r, n in zip(history["regions"], hp["densities"])
                )
                / volume
            )
            shares[name] = density * absorption(sigma) * sum(flux) / lattice_absorption
        totals = {c: 0.0 for c in ("tracked", "transported", "untransported")}
        for name, f in shares.items():
            totals[category(name, transportable)] += f
        captured = sum(shares[n] * ratios[n][burnup] for n in untransported)
        k_inf = state["k_inf"]
        ranked = sorted(
            (n for n in shares if category(n, transportable) == "untransported"),
            key=lambda n: -abs(shares[n]),
        )
        points.append(
            {
                "burnup": burnup,
                "untransported_all": {
                    n: {"fraction": shares[n], "captured_ratio": ratios[n][burnup]}
                    for n in untransported
                },
                "k_inf": k_inf,
                "fractions": totals,
                "rho_pcm": {c: -1e5 * f / k_inf for c, f in totals.items()},
                "untransported_captured": captured,
                "untransported_captured_rho_pcm": -1e5 * captured / k_inf,
                "untransported": [
                    {
                        "nuclide": n,
                        "fraction": shares[n],
                        "rho_pcm": -1e5 * shares[n] / k_inf,
                        "captured_ratio": ratios[n][burnup],
                    }
                    for n in ranked[:TOP]
                ],
                "transported": [
                    {
                        "nuclide": n,
                        "fraction": shares[n],
                        "rho_pcm": -1e5 * shares[n] / k_inf,
                    }
                    for n in sorted(
                        (
                            n
                            for n in shares
                            if category(n, transportable) == "transported"
                        ),
                        key=lambda n: -abs(shares[n]),
                    )[:TOP]
                ],
            }
        )
    return prune(points, ranking(points), PRUNE_PCM)


def prune(points, ranked, floor):
    """最大の目安が `floor` pcm 以下の核種を、点ごとの核種の一覧と順位から外す（順位の残りの和は外す前の値のまま）"""
    kept = {r["nuclide"] for r in ranked["ranking"] if r["max_pcm"] > floor}
    dropped = [r for r in ranked["ranking"] if r["nuclide"] not in kept]
    for p in points:
        p["untransported_all"] = {
            n: v for n, v in p["untransported_all"].items() if n in kept
        }
    return {
        "points": points,
        "ranking": [r for r in ranked["ranking"] if r["nuclide"] in kept],
        "counts": ranked["counts"],
        "pruned": {
            "count": len(dropped),
            "sum_pcm": sum(r["max_pcm"] for r in dropped),
        },
    }


def ranking(points):
    """ADR-0165 の案 S3 の形の順位: 核種ごとに、自分の吸収で減らした割合の燃焼度の点を通じた最大と、その反応度の目安 1e5 f / k∞ の
    最大。最大の大きい順に並べ、上から足したときの残りの核種の目安の和（核種ごとの最大の和で、保守的な側）を添える"""
    names = list(points[0]["untransported_all"])
    best = {}
    for name in names:
        values = [
            (
                1e5
                * p["untransported_all"][name]["fraction"]
                * p["untransported_all"][name]["captured_ratio"]
                / p["k_inf"],
                p["untransported_all"][name]["fraction"]
                * p["untransported_all"][name]["captured_ratio"],
                p["burnup"],
            )
            for p in points
        ]
        pcm, _, burnup = max(values)
        best[name] = {
            "pcm": pcm,
            "fraction": max(v[1] for v in values),
            "burnup": burnup,
        }
    order = sorted(names, key=lambda n: -best[n]["pcm"])
    rest = sum(best[n]["pcm"] for n in order)
    rows = []
    for name in order:
        rest -= best[name]["pcm"]
        rows.append(
            {
                "nuclide": name,
                "max_pcm": best[name]["pcm"],
                "max_fraction": best[name]["fraction"],
                "burnup": best[name]["burnup"],
                "rest_pcm": max(rest, 0.0),
                "xs_excluded": name in XS_EXCLUDED,
            }
        )
    counts = {
        str(t): next(i + 1 for i, r in enumerate(rows) if r["rest_pcm"] <= t)
        for t in THRESHOLDS_PCM
    }
    return {"ranking": rows, "counts": counts}


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        "燃焼度 (MWd/kg) | k∞ | 追う 31 核種 (%) | ほかの輸送する核種 (%) | 輸送しない核種 (%) | その Δρ の目安 (pcm) | "
        "自分の吸収で減らした輸送しない核種 (%) | その Δρ の目安 (pcm)"
    )
    for p in summary["points"]:
        f = p["fractions"]
        lines.append(
            f"{p['burnup']:g} | {p['k_inf']:.6f} | {100 * f['tracked']:.3f} | {100 * f['transported']:.3f} | "
            f"{100 * f['untransported']:.4f} | {p['rho_pcm']['untransported']:+.1f} | "
            f"{100 * p['untransported_captured']:.4f} | {p['untransported_captured_rho_pcm']:+.1f}"
        )
    lines.append(
        "足す核種の数（残りの目安の和が閾値以下になる順位）: "
        + ", ".join(f"{t} pcm で {n}" for t, n in summary["counts"].items())
    )
    lines.append(
        "順位 | 核種 | 最大の目安 (pcm) | その割合 (%) | 燃焼度 | 残りの和 (pcm) | xs が輸送する核種にしない"
    )
    for i, r in enumerate(summary["ranking"][:40]):
        lines.append(
            f"{i + 1} | {r['nuclide']} | {r['max_pcm']:.1f} | {100 * r['max_fraction']:.4f} | {r['burnup']:g} | "
            f"{r['rest_pcm']:.1f} | {'#0127' if r['xs_excluded'] else ''}"
        )
    last = summary["points"][-1]
    for key, label in [
        ("untransported", "輸送しない核種"),
        ("transported", "ほかの輸送する核種"),
    ]:
        lines.append(f"{label}の上位（{last['burnup']:g} MWd/kg）:")
        for n in last[key]:
            ratio = (
                f"（自分の吸収で減らすと × {n['captured_ratio']:.3f}）"
                if "captured_ratio" in n
                else ""
            )
            lines.append(
                f"  {n['nuclide']} {100 * n['fraction']:.4f}% {n['rho_pcm']:+.1f} pcm{ratio}"
            )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument(
        "--chain", help="the depletion chain (msgpack) for the half-lives"
    )
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not args.json or not args.chain or len(args.files) != 2:
        raise SystemExit(__doc__)
    history, table = (read(f) for f in args.files)
    half_lives = {
        n["name"]: n["decay"]["half_life"]
        for n in read(args.chain)["nuclides"]
        if n["decay"] is not None and n["decay"].get("half_life")
    }
    summary = {
        "format": "xlc-untransported-absorption",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in [*args.files, args.chain]
        ],
        **summarize(history, table, half_lives),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
