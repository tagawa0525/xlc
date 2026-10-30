"""ミクロ燃焼の線形の補正（ADR-0216）の誤差を、#0191 の候補（自己遮蔽の数密度への依存、数密度の差の空間の分布、少数群の重み）に分ける。

    python3 compare_split.py --json split.json --single table-split-pu240.msgpack   # このディレクトリの出力から JSON を書く
    python3 compare_split.py --report split.json  # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは 5 つの分岐のテーブル（`xlc-table`。燃焼度の軸と 1 点の分岐の軸、点ごとに状態 1 つ、どれも `micro` を持つ）:

- nominal（`table-nominal.msgpack`、nominal.toml）: 2A の履歴の点
- shifted（`table-split-shifted.msgpack`）: 点 k の組成を、点 k - 1 の組成のうち追う核種だけを点 k の値にしたもの（make_shifted.py）
- fixed（`table-split-fixed.msgpack`）: shifted の組成で、燃焼の領域の自己遮蔽を点 k - 1 の組成で解いたもの（`shielding-file`）
- spread（`table-split-spread.msgpack`）: 追う核種の点 k の平均を保ち、領域の間の分布を点 k - 1 のものにした組成を、自己遮蔽を
  点 k - 1 の組成で解いたもの（make_shifted.py --kind spread）
- spread_live（`table-split-spread-live.msgpack`）: spread の組成を、状態の組成で自己遮蔽したもの

点の組 (k - 1, k) ごとに、輸送の効き t_X = ρ(X の点 k の k∞) − ρ(nominal の点 k - 1 の k∞) と、線形の効き（compare_linearization.py
と同じ）から、誤差（線形 − t_shifted）を次の和に分ける（隣り合う計算の違いは置き換えた量 1 つだけ）:

- 自己遮蔽: t_fixed − t_shifted。燃焼の領域の実効断面積と SPH 因子が組成の差で変わる分（自分と他の核種の数密度への依存）
- 空間の分布: t_spread − t_fixed。追う核種の数密度の差の領域の間の分布が、σᵢ の重みにした数密度の分布と違う分
- 少数群の重み: 線形 − t_spread。残りの、組成の差による中性子束（少数群の中のスペクトルと領域の間の分布）の変化の分

`--single` のテーブル（追う核種 1 つだけを `micro` に持つ。split-pu240.toml）は、その核種の数密度の差だけの組成を、自己遮蔽と空間の
分布を固定して解いたもので、その核種の線形の補正（点 k - 1 の σᵢ × ΔN̄ᵢ）の誤差を各行の `single` に出す。少数群の重みの分のうち、
その核種の数密度の差が起こす分（ほかの核種の数密度の差との交差を含まない）。その誤差も、下と同じ内訳（その核種の項と残り）に分ける。

順を入れ替えた分け方（空間の分布を状態の組成の自己遮蔽で t_spread_live − t_shifted、自己遮蔽を spread の組成で t_spread − t_spread_live）
も出し、2 つの候補の交差の大きさを示す。

それぞれの分は、2 群の少数群定数の差を核種ごとに分けて内訳にする。比べる 2 つの状態 A、B の追う核種の均質化した数密度 N̄ᵢ は同じ
（自己遮蔽は組成が同じ。空間の分布は領域ごとの数密度が違い、平均だけが同じ。少数群の重みは A が線形の予測で、追う核種の分は点 k の
N̄ᵢ × 点 k - 1 の σᵢ）。σᵢ は均質化した反応率を N̄ᵢ と中性子束で割ったものなので、B − A は追う核種 i の N̄ᵢ (σᵢ,B − σᵢ,A) と、残り（追わない
核種と燃焼させない物質の巨視的断面積の差）と、χ の差の和にちょうどなる。各項の効きは ρ(A) − ρ(A + その項)（2 群の k）で、行の値（表の
k∞ から求めた分）と各項の和の差を交差の項として出す（2 群の k とテーブルの k∞ の差 `self_check_pcm` もここに入る）。

ρ は反応度 1 − 1/k で、差は pcm（Δρ = (1/k_a − 1/k_b) × 1e5、a から b へ）。
"""

import argparse
import json
from pathlib import Path

from compare_linearization import burnups, corrected, k_inf, md5, read, rho

FILES = {
    "nominal": "table-nominal.msgpack",
    "shifted": "table-split-shifted.msgpack",
    "fixed": "table-split-fixed.msgpack",
    "spread": "table-split-spread.msgpack",
    "spread_live": "table-split-spread-live.msgpack",
}
# 燃焼度 0 の組成は 5 つの計算で同じなので、k∞ は同じ問題の値（自己遮蔽を前の点で解く計算も、燃焼度 0 の前の点は同じ組成）
FRESH_TOLERANCE_PCM = 1e-6
# 追う核種の均質化した数密度は、置き換えた 4 つの組成で点 k の値（spread は平均を保つ書き換えなので丸めの差だけ）
DENSITY_TOLERANCE = 1e-12
# 内訳に核種ごとに出す閾値（どれかの組のどれかの分で効きの絶対値がこれ以上の追う核種。ほかは「ほかの追う核種」にまとめる）
NUCLIDE_THRESHOLD_PCM = 5.0


def add(state, delta, sign=1.0):
    """少数群定数 `state` に差 `delta`（total、nu_fission、scatter）を `sign` 倍して足したもの"""
    out = dict(state)
    out["total"] = [a + sign * d for a, d in zip(state["total"], delta["total"])]
    out["nu_fission"] = [
        a + sign * d for a, d in zip(state["nu_fission"], delta["nu_fission"])
    ]
    out["scatter"] = [
        [a + sign * d for a, d in zip(row, drow)]
        for row, drow in zip(state["scatter"], delta["scatter"])
    ]
    return out


def difference(b, a):
    """少数群定数の差 b − a（total、nu_fission、scatter）"""
    return {
        "total": [x - y for x, y in zip(b["total"], a["total"])],
        "nu_fission": [x - y for x, y in zip(b["nu_fission"], a["nu_fission"])],
        "scatter": [
            [x - y for x, y in zip(rb, ra)]
            for rb, ra in zip(b["scatter"], a["scatter"])
        ],
    }


def nuclide_delta(density, sigma_b, sigma_a):
    """追う核種の巨視的断面積の差 N (σ_b − σ_a)"""
    return {
        "total": [
            density * (x - y) for x, y in zip(sigma_b["total"], sigma_a["total"])
        ],
        "nu_fission": [
            density * (x - y)
            for x, y in zip(sigma_b["nu_fission"], sigma_a["nu_fission"])
        ],
        "scatter": [
            [density * (x - y) for x, y in zip(rb, ra)]
            for rb, ra in zip(sigma_b["scatter"], sigma_a["scatter"])
        ],
    }


def split(a, b, sigmas_a, sigmas_b, densities, names, whole):
    """ρ(A) − ρ(B) の分 `whole`（行の値）を、追う核種ごと、残り（others）、χ に分ける。B − A のうち追う核種 i の分は
    N̄ᵢ (σᵢ,B − σᵢ,A)、残りは total、nu_fission、scatter の差のうちその和でない分、χ は A の χ を B の χ にした分（それぞれ 2 群の k）。
    交差の項（cross）は `whole` と各項の和の差"""
    rho_a = rho(k_inf(a))
    rest = difference(b, a)
    nuclides = {}
    for name, n, sa, sb in zip(names, densities, sigmas_a, sigmas_b):
        delta = nuclide_delta(n, sb, sa)
        rest = add(rest, delta, -1.0)
        nuclides[name] = rho_a - rho(k_inf(add(a, delta)))
    others = rho_a - rho(k_inf(add(a, rest)))
    chi = rho_a - rho(k_inf(dict(a, chi=b["chi"])))
    total = sum(nuclides.values()) + others + chi
    return {
        "whole_pcm": whole,
        "nuclides": nuclides,
        "others_pcm": others,
        "chi_pcm": chi,
        "cross_pcm": whole - total,
    }


def check_single(single, nominal, names):
    """1 核種のテーブルが nominal と同じ燃焼度の点で、nominal の追う核種の 1 つだけを追い、燃焼度 0 の k∞ と、その核種の均質化した
    数密度が nominal と同じこと。その核種の nominal の中の番号を返す"""
    points = burnups(nominal)
    if single.get("micro") is None or burnups(single) != points:
        raise SystemExit(
            "the single table needs micro data at the nominal burnup points"
        )
    tracked = single["micro"]["nuclides"]
    if len(tracked) != 1 or tracked[0] not in names:
        raise SystemExit(
            f"the single table must track one nuclide of the nominal table, got {tracked}"
        )
    i = names.index(tracked[0])
    fresh = rho(single["states"][0]["k_inf"]) - rho(nominal["states"][0]["k_inf"])
    if abs(fresh) > FRESH_TOLERANCE_PCM:
        raise SystemExit(f"the fresh point of the single table differs by {fresh} pcm")
    for k, (got, expected) in enumerate(
        zip(single["micro"]["densities"], nominal["micro"]["densities"])
    ):
        if abs(got[0] - expected[i]) > DENSITY_TOLERANCE * abs(expected[i]):
            raise SystemExit(
                f"the single table has {tracked[0]} at {got[0]} for {expected[i]} at burnup {points[k]}"
            )
    return i


def check(tables):
    """5 つのテーブルが同じ燃焼度の点、同じ追う核種を持ち、燃焼度 0 の k∞ が同じで、置き換えた組成の追う核種の均質化した数密度が
    nominal の点 k の値と同じこと"""
    nominal = tables["nominal"]
    points = burnups(nominal)
    names = None
    for label, table in tables.items():
        if table.get("micro") is None:
            raise SystemExit(f"the {label} table has no micro data")
        if burnups(table) != points:
            raise SystemExit(f"the {label} table has other burnup points")
        if names is None:
            names = table["micro"]["nuclides"]
        elif table["micro"]["nuclides"] != names:
            raise SystemExit(f"the {label} table tracks other nuclides")
        fresh = rho(table["states"][0]["k_inf"]) - rho(nominal["states"][0]["k_inf"])
        if abs(fresh) > FRESH_TOLERANCE_PCM:
            raise SystemExit(
                f"the fresh point of the {label} table differs by {fresh} pcm"
            )
        for k, (got, expected) in enumerate(
            zip(table["micro"]["densities"], nominal["micro"]["densities"])
        ):
            for name, g, e in zip(names, got, expected):
                if abs(g - e) > DENSITY_TOLERANCE * abs(e):
                    raise SystemExit(
                        f"the {label} table has {name} at {g} for {e} at burnup {points[k]}"
                    )
    return points, names


def compare(tables):
    tables = dict(tables)
    single = tables.pop("single", None)
    points, names = check(tables)
    one = None if single is None else check_single(single, tables["nominal"], names)
    nominal = tables["nominal"]
    micro = nominal["micro"]
    every = [*tables.values(), *([] if single is None else [single])]
    self_check = max(
        abs(rho(k_inf(s)) - rho(s["k_inf"])) for t in every for s in t["states"]
    )
    rows = []
    for k in range(1, len(points)):
        before, after = micro["densities"][k - 1], micro["densities"][k]
        deltas = [b - a for a, b in zip(before, after)]
        base = nominal["states"][k - 1]
        sigma0 = micro["states"][k - 1]
        predicted = corrected(base, sigma0, deltas)
        linear = rho(k_inf(predicted)) - rho(k_inf(base))
        t = {
            label: rho(tables[label]["states"][k]["k_inf"]) - rho(base["k_inf"])
            for label in ("shifted", "fixed", "spread", "spread_live")
        }
        extra = {}
        if single is not None:
            alone = corrected(base, [sigma0[one]], [deltas[one]])
            linear_alone = rho(k_inf(alone)) - rho(k_inf(base))
            transport_alone = rho(single["states"][k]["k_inf"]) - rho(base["k_inf"])
            extra["single"] = {
                "nuclide": names[one],
                "linear_pcm": linear_alone,
                "transport_pcm": transport_alone,
                "error_pcm": linear_alone - transport_alone,
                # A = その核種だけの線形の予測、B = 1 核種の組成
                "split": split(
                    alone,
                    single["states"][k],
                    [sigma0[one]],
                    single["micro"]["states"][k],
                    [after[one]],
                    [names[one]],
                    linear_alone - transport_alone,
                ),
            }
        state = {label: tables[label]["states"][k] for label in t}
        sigma = {label: tables[label]["micro"]["states"][k] for label in t}
        rows.append(
            {
                "from": points[k - 1],
                "to": points[k],
                "tracked_pcm": t["shifted"],
                "linear_pcm": linear,
                "error_pcm": linear - t["shifted"],
                "shielding_pcm": t["fixed"] - t["shifted"],
                "spatial_pcm": t["spread"] - t["fixed"],
                "weight_pcm": linear - t["spread"],
                "spatial_live_pcm": t["spread_live"] - t["shifted"],
                "shielding_spread_pcm": t["spread"] - t["spread_live"],
                "splits": {
                    # 自己遮蔽: A = fixed、B = shifted（ρ(fixed) − ρ(shifted)）
                    "shielding": split(
                        state["fixed"],
                        state["shifted"],
                        sigma["fixed"],
                        sigma["shifted"],
                        after,
                        names,
                        t["fixed"] - t["shifted"],
                    ),
                    # 空間の分布: A = spread、B = fixed
                    "spatial": split(
                        state["spread"],
                        state["fixed"],
                        sigma["spread"],
                        sigma["fixed"],
                        after,
                        names,
                        t["spread"] - t["fixed"],
                    ),
                    # 少数群の重み: A = 線形の予測（σ は点 k - 1 の値）、B = spread
                    "weight": split(
                        predicted,
                        state["spread"],
                        sigma0,
                        sigma["spread"],
                        after,
                        names,
                        linear - t["spread"],
                    ),
                },
                **extra,
            }
        )
    shown = [
        n
        for n in names
        if any(
            abs(r["splits"][part]["nuclides"][n]) >= NUCLIDE_THRESHOLD_PCM
            for r in rows
            for part in r["splits"]
        )
    ]
    for r in rows:
        for part in r["splits"].values():
            every = part.pop("nuclides")
            part["nuclides"] = {n: every[n] for n in shown}
            part["other_tracked_pcm"] = sum(
                v for n, v in every.items() if n not in shown
            )
    return {
        "format": "xlc-compare-micro-linearization-split",
        "nuclides": len(names),
        "shown_nuclides": shown,
        "self_check_pcm": self_check,
        "rows": rows,
    }


def report(summary):
    rows = summary["rows"]

    def span(key):
        values = [r[key] for r in rows]
        return f"{min(values):+.2f}〜{max(values):+.2f}"

    def share(key, after):
        """5 MWd/kg より後（after）か前の組の、誤差に対する割合の範囲（%）"""
        values = [
            100 * r[key] / r["error_pcm"]
            for r in rows
            if (r["from"] >= 5.0) == after and abs(r["error_pcm"]) > 1.0
        ]
        return f"{min(values):.1f}〜{max(values):.1f}%"

    def of_effect(key):
        """5 MWd/kg より後の組の、輸送の効きに対する割合の範囲（%）"""
        values = [100 * r[key] / r["tracked_pcm"] for r in rows if r["from"] >= 5.0]
        return f"{min(values):.2f}〜{max(values):.2f}%"

    cross = max(abs(p["cross_pcm"]) for r in rows for p in r["splits"].values())
    order = max(abs(r["spatial_live_pcm"] - r["spatial_pcm"]) for r in rows)
    lines = [
        (
            f"追う核種 {summary['nuclides']}、点の組 {len(rows)}。誤差 {span('error_pcm')} pcm のうち、"
            f"自己遮蔽 {span('shielding_pcm')} pcm、空間の分布 {span('spatial_pcm')} pcm、少数群の重み {span('weight_pcm')} pcm"
        ),
    ]
    for key, label in (
        ("shielding_pcm", "自己遮蔽"),
        ("spatial_pcm", "空間の分布"),
        ("weight_pcm", "少数群の重み"),
    ):
        lines.append(
            f"{label}: 誤差に対して 5 MWd/kg より後の組で {share(key, True)}、前の組で {share(key, False)}。"
            f"5 MWd/kg より後の組の輸送の効きに対して {of_effect(key)}"
        )
    lines.append(
        f"順を入れ替えた空間の分布との差は最大 {order:.2f} pcm（2 つの候補の交差）。内訳の交差の項は最大 {cross:.2f} pcm。"
        f"2 群の k とテーブルの k∞ の差は最大 {summary['self_check_pcm']:.2e} pcm"
    )
    if "single" in rows[0]:
        name = rows[0]["single"]["nuclide"]
        late = [r for r in rows if r["from"] >= 5.0]
        gap = max(abs(r["single"]["error_pcm"] - r["weight_pcm"]) for r in late)
        errors = [r["single"]["error_pcm"] for r in rows]
        lines.append(
            f"{name} だけの組成の誤差は {min(errors):+.2f}〜{max(errors):+.2f} pcm で、5 MWd/kg より後の組の少数群の重みの分との差は"
            f"最大 {gap:.2f} pcm"
        )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", help="write the summary to this file")
    parser.add_argument("--report", help="print the report of this summary")
    parser.add_argument(
        "--single", help="the table of a single-nuclide composition (one micro nuclide)"
    )
    args = parser.parse_args()
    if args.report:
        with open(args.report, encoding="utf-8") as f:
            print(report(json.load(f)))
        return
    here = Path(__file__).resolve().parent
    files = {name: here / file for name, file in FILES.items()}
    if args.single:
        files["single"] = here / args.single
    tables = {name: read(path) for name, path in files.items()}
    summary = compare(tables)
    summary["inputs"] = {name: md5(path) for name, path in files.items()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
