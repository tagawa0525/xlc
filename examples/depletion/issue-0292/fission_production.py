"""無限希釈の吸収体として数えた核分裂する核種（ADR-0261）の、落とした核分裂の 1 次の見積もりを出し、輸送する核種に足す核種を選ぶ（#0292）。

使い方: python3 fission_production.py --json <出力> --xs <ライブラリの xs.toml> --added <足した核種の一覧> --shielded <自己遮蔽する核種の一覧>
          --fuel <名前> <前> <テーブル> <希釈> <遮蔽の希釈> <比較の JSON> [--fuel ...]
        python3 fission_production.py --report <出力の JSON>

ENDF の ν̄（MF1 MT452）を読むので、`references` の devShell の python3（endf パッケージ）で実行する。

燃料ごとの入力（compare_dilute.py と同じ燃焼の履歴）:

- 前: p2 のライブラリで数えない燃焼の履歴。テーブルの格子の体積を決める（untransported_absorption.py の `lattice_volume`）
- テーブル: 前の履歴からのブランチ計算。燃焼度の点ごとの冷却計算の 1 群の核分裂の断面積と、2 群の中性子束と吸収
- 希釈: p2 のライブラリで全核種を無限希釈で数えた燃焼の履歴。足す核種の順位に使う
- 遮蔽の希釈: 自己遮蔽する核種だけを輸送する核種にして残りを数えた燃焼の履歴。足した核種のうち自己遮蔽しない核種の見積もりに使い、
  比較の JSON の `dilute_unshielded_pcm`（遮蔽の希釈 − 足した希釈）と並べる

ADR-0261 の形は、無限希釈の吸収体の捕獲だけを輸送に入れ、核分裂は吸収にも生成にも入れない。無限媒質の k = P / A で、核分裂の率 F を
吸収と生成 ν F の両方から外すと、1 次で Δρ = 1/k_数える − 1/k_外す = −(F/A)(ν − k)/k²（P = k A）。F/A は
untransported_absorption.py の割合と同じ形で、履歴の数密度を格子で均質化した N̄ と、テーブルの 1 群の σf、2 群の中性子束の和 φ̄ から
N̄ σf φ̄ / Σ_G Σₐ,G φ̄_G。k は履歴の同じ燃焼度の k∞。スペクトルと自己遮蔽の変化、核分裂の生成物と数密度の変化は含まない。

ν は 2 つ持つ: 熱中性子（0.0253 eV）の ν̄ と、評価の全エネルギーでの ν̄ の最大。熱の値は見積もりに、最大は足す核種の順位に使う
（核分裂の起きるエネルギーによらず (ν − k) の上限になり、順位の見積もりは保守的な側になる）。

足す核種の選び方は untransported_absorption.py の順位（ADR-0165 の案 S3）と同じ形: 核分裂する吸収体の核種ごとに、燃料と燃焼度を
通じた見積もりの絶対値の最大で並べ、上から足して、残りの最大の和が項目ごとの目安 5 pcm（ADR-0200、ADR-0215）以下になるまでを選ぶ。
"""

import argparse
import hashlib
import json
from pathlib import Path

import tomllib
from untransported_absorption import (
    history_point,
    lattice_volume,
    nominal_states,
    read,
    region_sum,
    two_group_flux,
)

# 項目ごとの目安 (pcm)。ADR-0200、ADR-0215
THRESHOLD_PCM = 5.0
THERMAL_EV = 0.0253
# 多項式の ν̄（LNU = 1）の最大を探すエネルギーの範囲 (eV) と点の数（対数で等間隔）
POLYNOMIAL_RANGE = (1e-5, 2e7)
POLYNOMIAL_POINTS = 200


def dropped_production_pcm(fraction, nu, k):
    """核分裂の率が格子の吸収の `fraction` の核種の、核分裂を吸収と生成の両方から外した Δρ の 1 次の見積もり (pcm)"""
    return -1e5 * fraction * (nu - k) / k**2


def nubar(section):
    """MF1 MT452 の節（endf パッケージの辞書）から、(熱中性子の ν̄, 全エネルギーでの ν̄ の最大)"""
    if section["LNU"] == 2:
        nu = section["nu"]
        return float(nu(THERMAL_EV)), float(max(nu.y))
    if section["LNU"] == 1:
        c = section["C"]

        def poly(e):
            return sum(ci * e**i for i, ci in enumerate(c))

        lo, hi = POLYNOMIAL_RANGE
        grid = [
            lo * (hi / lo) ** (i / (POLYNOMIAL_POINTS - 1))
            for i in range(POLYNOMIAL_POINTS)
        ]
        return poly(THERMAL_EV), max(poly(e) for e in grid)
    raise SystemExit(f"unknown LNU {section['LNU']} in MF1 MT452")


def fission_fractions(before, table, history):
    """燃焼度の点（0 を除く）ごとに、履歴の無限希釈の吸収体のうち核分裂の断面積を持つ核種の F/A と、履歴の k∞"""
    if not history.get("dilute_absorbers"):
        raise SystemExit("the history counts no dilute absorbers")
    burnups, states = nominal_states(table)
    cooling = table["cooling"]
    column = {n: j for j, n in enumerate(cooling["nuclides"])}
    missing = [n for n in history["dilute_absorbers"] if n not in column]
    if missing:
        raise SystemExit(f"the table has no cooling cross sections of {missing[:5]}")
    index = {n: i for i, n in enumerate(history["nuclides"])}
    fissionable = [
        n
        for n in history["dilute_absorbers"]
        if any(p[column[n]]["fission"] > 0.0 for p in cooling["points"])
    ]
    points = []
    for k, burnup in enumerate(burnups):
        if burnup == 0.0:
            continue
        state = table["states"][states[k]]
        flux = two_group_flux(state)
        lattice_absorption = sum(a * f for a, f in zip(state["absorption"], flux))
        volume = lattice_volume(before, table, k, burnup)
        hp = history_point(history, burnup)
        fractions = {
            n: region_sum(history, hp, index[n])
            / volume
            * cooling["points"][k][column[n]]["fission"]
            * sum(flux)
            / lattice_absorption
            for n in fissionable
        }
        points.append({"burnup": burnup, "k_inf": hp["k_inf"], "fractions": fractions})
    return points


def select(estimates, threshold):
    """`estimates` は核種ごとの (燃料, 燃焼度, 見積もり pcm) の列。絶対値の最大で並べ、残りの最大の和が `threshold` 以下になるまでを選ぶ"""
    best = {n: max(rows, key=lambda r: abs(r[2])) for n, rows in estimates.items()}
    order = sorted(best, key=lambda n: -abs(best[n][2]))
    rest = sum(abs(best[n][2]) for n in order)
    selected = []
    ranking = []
    for name in order:
        if rest > threshold:
            selected.append(name)
        rest -= abs(best[name][2])
        fuel, burnup, pcm = best[name]
        ranking.append(
            {
                "nuclide": name,
                "max_pcm": pcm,
                "fuel": fuel,
                "burnup": burnup,
                "rest_pcm": max(rest, 0.0),
            }
        )
    return {"ranking": ranking, "selected": selected}


def estimate(points, nus, names, measured):
    """燃焼度の点ごとに、`names` の核種の見積もり（熱中性子の ν̄）とその和を、比較の JSON の測った差と並べる"""
    by_burnup = {p["burnup"]: p for p in measured}
    rows = []
    for p in points:
        if p["burnup"] not in by_burnup:
            raise SystemExit(f"the comparison has no burnup {p['burnup']}")
        nuclides = {
            n: dropped_production_pcm(p["fractions"][n], nus[n][0], p["k_inf"])
            for n in names
            if n in p["fractions"]
        }
        rows.append(
            {
                "burnup": p["burnup"],
                "nuclides": nuclides,
                "estimate_pcm": sum(nuclides.values()),
                "measured_pcm": by_burnup[p["burnup"]]["dilute_unshielded_pcm"],
            }
        )
    return rows


def read_nus(xs_toml, names):
    import endf

    config = tomllib.loads(Path(xs_toml).read_text())
    directory = Path(xs_toml).parent / config["endf-directory"]
    out = {}
    for name in names:
        material = endf.Material(str(directory / config["nuclide"][name]["endf"]))
        out[name] = nubar(material.section_data[1, 452])
    return out


def summarize(fuels, xs_toml, added, shielded):
    ranked = {}
    unshielded = {}
    for name, (before, table, dilute, shielded_dilute, compared) in fuels.items():
        if set(compared["shielded"]) != set(shielded) or set(compared["added"]) != set(
            added
        ):
            raise SystemExit(f"{name}: the comparison used other nuclide lists")
        ranked[name] = fission_fractions(before, table, dilute)
        unshielded[name] = fission_fractions(before, table, shielded_dilute)
    names = sorted({n for pts in ranked.values() for p in pts for n in p["fractions"]})
    nus = read_nus(xs_toml, names)
    estimates = {
        n: [
            (
                fuel,
                p["burnup"],
                dropped_production_pcm(p["fractions"][n], nus[n][1], p["k_inf"]),
            )
            for fuel, pts in ranked.items()
            for p in pts
        ]
        for n in names
    }
    chosen = select(estimates, THRESHOLD_PCM)
    for row in chosen["ranking"]:
        row["nu_thermal"], row["nu_max"] = nus[row["nuclide"]]
    unshielded_names = sorted(set(added) - set(shielded))
    return {
        "threshold_pcm": THRESHOLD_PCM,
        **chosen,
        "unshielded": {
            fuel: estimate(
                unshielded[fuel], nus, unshielded_names, fuels[fuel][4]["points"]
            )
            for fuel in fuels
        },
    }


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    lines.append(
        f"足す核種（残りの和が {summary['threshold_pcm']:g} pcm 以下）: {' '.join(summary['selected'])}"
    )
    lines.append(
        "核種 | ν̄（熱） | ν̄（最大） | 最大の見積もり (pcm) | 燃料 | 燃焼度 | 残りの和 (pcm)"
    )
    for r in summary["ranking"]:
        if abs(r["max_pcm"]) < 0.01:
            continue
        lines.append(
            f"{r['nuclide']} | {r['nu_thermal']:.3f} | {r['nu_max']:.3f} | {r['max_pcm']:.2f} | {r['fuel']} | "
            f"{r['burnup']:g} | {r['rest_pcm']:.2f}"
        )
    for fuel, rows in summary["unshielded"].items():
        lines.append(f"{fuel}: 燃焼度 | 見積もり (pcm) | 測った差 (pcm) | 上位の核種")
        for r in rows:
            top = sorted(r["nuclides"].items(), key=lambda x: x[1])[:3]
            lines.append(
                f"  {r['burnup']:g} | {r['estimate_pcm']:+.1f} | {r['measured_pcm']:+.1f} | "
                + ", ".join(f"{n} {v:+.1f}" for n, v in top)
            )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument("--xs")
    parser.add_argument("--added")
    parser.add_argument("--shielded")
    parser.add_argument("--fuel", nargs=6, action="append", metavar="FILE")
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not (args.json and args.xs and args.added and args.shielded and args.fuel):
        raise SystemExit(__doc__)
    fuels = {}
    paths = [args.xs, args.added, args.shielded]
    for name, *files in args.fuel:
        *histories, compared = files
        fuels[name] = (
            *(read(f) for f in histories),
            json.loads(Path(compared).read_text()),
        )
        paths.extend(files)
    added = Path(args.added).read_text().split()
    shielded = Path(args.shielded).read_text().split()
    summary = {
        "format": "xlc-dropped-fission-production",
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in paths
        ],
        **summarize(fuels, args.xs, added, shielded),
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
