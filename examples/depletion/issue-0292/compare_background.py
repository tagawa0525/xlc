"""同じ組成で、輸送しない核種の扱いだけを変えた状態の k∞ の差を、核種、効果（吸収、散乱、生成）、エネルギーの範囲に分ける（#0292、
ADR-0261）。

使い方: python3 compare_background.py --json <出力> --burnup <燃焼度> --case <名前> <ブランチ計算の出力> <履歴> <ライブラリ> ...
        python3 compare_background.py --report <出力の JSON>

状態（`--case`）は、履歴からのブランチ計算（`[history]`）の 1 状態で、少数群の境界をライブラリの群の境界そのものにし、`micro` に
輸送するすべての核種を書いたもの（background_inputs.py が書く）。燃焼の領域は 1 つに限る（無限希釈で数える核種の群ごとの吸収を、
下の中性子束の比で作るため）。並べた順の隣り合う 2 つ（A → B）ごとに、差 Δρ(B − A) = (1/k_A − 1/k_B) × 1e5 (pcm。正は B の
k が高い向き)を分ける。

格子で均質化した群定数の無限媒質（外周が鏡面反射の格子の k と同じ問題）で、M = diag(Σt) − Sᵀ（S は [入射の群][出射の群]）、
F = χ (νΣf)ᵀ、λ = 1/k は M φ = λ F φ の固有値。A の随伴 φ*_A と B の φ_B で

    λ_B − λ_A = φ*_Aᵀ (δM − λ_B δF) φ_B / (φ*_Aᵀ F_A φ_B)

は 1 次の近似でなく恒等式（δM = M_B − M_A、δF = F_B − F_A）で、δM と δF を項に分けた和に分かれる。δF は χ_A (δνΣf)ᵀ と
δχ (νΣf_B)ᵀ に分ける。均質化した群定数を、核種ごとの寄与 Nᵢ σᵢ と残り（燃焼させない物質。SPH 因子と中性子束の重みの差を含む）に
分け、核種ごとの差を 3 つの効果に分ける:

- 吸収: diag(δ(Σt − Σ_h S_gh))。正味の吸収で、(n,2n) は負に入る
- 散乱: diag(δ Σ_h S_gh) − δSᵀ（群の中の散乱は打ち消し合う）
- 生成: δνΣf

項は入射の群（生成は核分裂を起こす群）で、ライブラリの UFG の範囲（`pointwise_range`）より上（fast）、範囲の中（ufg）、範囲より下
（thermal）に分ける。核種の寄与は、輸送する核種は `micro` の N̄ σ（全断面積、P0 の散乱行列、νσf）、無限希釈で数える核種は、
燃焼の領域の数密度を格子で均質化した N̄ と、ライブラリの放射化の断面積の (n,γ)、(n,p)、(n,α)（lattice の輸送が数える反応。
crates/lattice/src/dilute.rs の `ABSORPTION`）の群の値に、燃料の領域の中性子束の格子の中性子束に対する比 R_g を掛けた吸収だけ。
R_g は燃料の O-16 の `micro` の全断面積 / ライブラリの無限希釈の全断面積（O-16 は燃料の中で共鳴の自己遮蔽がほぼ無い）。O-16 は
MeV の群で共鳴を持ち R_g がずれうるので、無限希釈で数える核種の吸収の和を、冷却計算の 1 群の断面積（`cooling`。lattice が
燃料の領域の中性子束で作る）からの和と比べ、相対 `DILUTE_TOLERANCE` を超えればエラーにする。

核種の種類: 両方で輸送する（transported）、A で輸送し B で数える（to_dilute）、その逆（to_transported）、両方で数える（dilute）。
`parts` は差を次の和に分ける（和は差に一致する）:

- u238_resonance: 両方で輸送する U-238 の ufg の吸収。同じ組成では、U-238 の実効断面積の差（共鳴の自己遮蔽の背景の差）
- resonance_other: 両方で輸送するほかの核種の ufg の吸収
- absorption_other_ranges: 両方で輸送する核種の fast と thermal の吸収
- transported_scattering_production: 両方で輸送する核種の散乱と生成
- moved_absorption: 扱いを変えた核種（to_dilute、to_transported）の吸収
- moved_scattering_fast、moved_scattering_below: 同じ核種の散乱の、fast と、ufg と thermal
- moved_production: 同じ核種の生成
- rest: 両方で数える核種、残り、χ の差
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from untransported_absorption import read

FORMAT = "xlc-compare-background"
# lattice の輸送が無限希釈の吸収体の吸収に数える反応（crates/lattice/src/dilute.rs の ABSORPTION）
ABSORPTION = (102, 103, 107)
# 無限希釈の吸収の和の、冷却計算の 1 群の断面積からの和に対する相対の許容。R_g のずれ（O-16 の MeV の共鳴）が出る大きさで、U と
# MOX のピンでは 5e-4 以下（examples/depletion/verification.md の「無限希釈で数える核種の残る差の分け方」）
DILUTE_TOLERANCE = 1e-3
RANGES = ("fast", "ufg", "thermal")
EFFECTS = ("absorption", "scattering", "production")
PARTS = (
    "u238_resonance",
    "resonance_other",
    "absorption_other_ranges",
    "transported_scattering_production",
    "moved_absorption",
    "moved_scattering_fast",
    "moved_scattering_below",
    "moved_production",
    "rest",
)
# JSON の核種ごとの項から外す、すべての項の絶対値の最大の床 (pcm)。外した核種の数と和は JSON に残す
PRUNE_PCM = 0.005


def at_temperature(entry, temperature):
    matched = [t for t in entry["temperatures"] if t["temperature"] == temperature]
    if len(matched) != 1:
        raise SystemExit(f"{entry['name']}: no library data at {temperature} K")
    return matched[0]


def range_masks(bounds, pointwise):
    """群ごとの範囲（fast、ufg、thermal）の真偽。群の境界は高いエネルギーから"""
    low, high = pointwise
    if low not in bounds or high not in bounds:
        raise SystemExit(f"the UFG range {pointwise} is not on the group bounds")
    upper, lower = np.array(bounds[:-1]), np.array(bounds[1:])
    return {
        "fast": lower >= high,
        "ufg": (lower >= low) & (upper <= high),
        "thermal": upper <= low,
    }


def case(branch, history, library, burnup):
    """1 状態の、格子で均質化した群定数と、核種ごとの寄与 (Σt, S, νΣf)"""
    if len(branch["states"]) != 1:
        raise SystemExit("the branch output must have one state")
    if len(history["regions"]) != 1:
        raise SystemExit("the history must have one depletion region")
    if branch["bounds"] != library["group_bounds"]:
        raise SystemExit("the few groups must be the library groups")
    (state,) = branch["states"]
    total = np.array(state["total"])
    scatter = np.array(state["scatter"])
    chi = np.array(state["chi"])
    nu_fission = np.array(state["nu_fission"])
    groups = len(total)
    (temperature,) = next(a for a in branch["axes"] if a["name"] == "fuel")["points"]
    micro = branch["micro"]
    names = micro["nuclides"]
    if names != history["transportable"]:
        raise SystemExit(
            "the micro cross sections must list every transported nuclide of the history in its order"
        )
    densities = dict(zip(names, micro["densities"][0]))
    matched = [p for p in history["points"] if p["burnup"] == burnup]
    if len(matched) != 1:
        raise SystemExit(f"burnup {burnup} is not a point of the history")
    (region,) = matched[0]["densities"]
    index = {n: i for i, n in enumerate(history["nuclides"])}
    # 燃焼の領域の体積 / 格子の体積
    ratio = densities["U-238"] / region[index["U-238"]]
    contributions = {}
    for n, x in zip(names, micro["states"][0]):
        d = densities[n]
        contributions[n] = (
            d * np.array(x["total"]),
            d * np.array(x["scatter"]),
            d * np.array(x["nu_fission"]),
        )
    (oxygen,) = [x for x in library["nuclides"] if x["name"] == "O-16"]
    infinite = np.array(at_temperature(oxygen, temperature)["shielded"][0]["total"])
    flux_ratio = np.array(micro["states"][0][names.index("O-16")]["total"]) / infinite
    activation = {a["name"]: a for a in library["activation"]}
    cooling = dict(zip(branch["cooling"]["nuclides"], branch["cooling"]["points"][0]))
    zero_scatter = np.zeros((groups, groups))
    dilute_one_group = 0.0
    for n in history.get("dilute_absorbers", []):
        d = max(region[index[n]], 0.0) * ratio
        if d == 0.0:
            continue
        sigma = np.zeros(groups)
        for reaction in at_temperature(activation[n], temperature)["reactions"]:
            if reaction["mt"] in ABSORPTION:
                for channel in reaction["channels"]:
                    sigma += np.array(channel["cross_section"])
        contributions[n] = (d * sigma * flux_ratio, zero_scatter, np.zeros(groups))
        one = sum(
            sum(r["channels"]) for r in cooling[n]["reactions"] if r["mt"] in ABSORPTION
        )
        dilute_one_group += d * one
    removal = np.diag(total) - scatter.T
    phi = np.linalg.solve(removal, chi)
    dilute_names = [
        n for n in history.get("dilute_absorbers", []) if n in contributions
    ]
    dilute_check = 0.0
    if dilute_names:
        groupwise = sum(contributions[n][0] for n in dilute_names) @ phi
        dilute_check = groupwise / (dilute_one_group * phi.sum()) - 1.0
        if abs(dilute_check) > DILUTE_TOLERANCE:
            raise SystemExit(
                f"the dilute absorption from the library differs from the cooling cross sections by {dilute_check:.2e}"
            )
    rest = (
        total - sum(c[0] for c in contributions.values()),
        scatter - sum(c[1] for c in contributions.values()),
        nu_fission - sum(c[2] for c in contributions.values()),
    )
    return {
        "k_state": state["k_inf"],
        "k_matrix": float(nu_fission @ phi),
        "k_history": matched[0]["k_inf"],
        "dilute_check": float(dilute_check),
        "transported": list(names),
        "dilute": dilute_names,
        "masks": range_masks(library["group_bounds"], library["pointwise_range"]),
        "total": total,
        "scatter": scatter,
        "chi": chi,
        "nu_fission": nu_fission,
        "contributions": contributions,
        "rest": rest,
    }


def kind(name, a, b):
    ta, tb = name in a["transported"], name in b["transported"]
    if ta and tb:
        return "transported"
    if ta:
        return "to_dilute"
    if tb:
        return "to_transported"
    return "dilute"


def step(a, b):
    """A → B の差と、その分け方"""
    groups = len(a["total"])
    removal_a = np.diag(a["total"]) - a["scatter"].T
    removal_b = np.diag(b["total"]) - b["scatter"].T
    adjoint = np.linalg.solve(removal_a.T, a["nu_fission"])
    forward = np.linalg.solve(removal_b, b["chi"])
    lambda_a = 1.0 / float(a["nu_fission"] @ np.linalg.solve(removal_a, a["chi"]))
    lambda_b = 1.0 / float(b["nu_fission"] @ forward)
    adjoint_chi = float(adjoint @ a["chi"])
    denominator = adjoint_chi * float(a["nu_fission"] @ forward)

    def term(d_total, d_scatter, d_production, d_chi=None):
        """Δρ(B − A) の項 (pcm)"""
        loss = float(adjoint @ (d_total * forward - d_scatter.T @ forward))
        gain = adjoint_chi * float(d_production @ forward)
        if d_chi is not None:
            gain += float(adjoint @ d_chi) * float(b["nu_fission"] @ forward)
        return -(loss - lambda_b * gain) / denominator * 1e5

    zero = np.zeros(groups)
    zero_scatter = np.zeros((groups, groups))
    masks = a["masks"]

    def effects(ca, cb):
        dt, ds, df = (y - x for x, y in zip(ca, cb))
        out = {}
        for effect in EFFECTS:
            out[effect] = {}
            for r in RANGES:
                m = masks[r]
                if effect == "absorption":
                    out[effect][r] = term((dt - ds.sum(axis=1)) * m, zero_scatter, zero)
                elif effect == "scattering":
                    out[effect][r] = term(ds.sum(axis=1) * m, ds * m[:, None], zero)
                else:
                    out[effect][r] = term(zero, zero_scatter, df * m)
        return out

    empty = (zero, zero_scatter, zero)
    names = sorted(set(a["contributions"]) | set(b["contributions"]))
    nuclides = {}
    for n in names:
        nuclides[n] = {
            "kind": kind(n, a, b),
            **effects(
                a["contributions"].get(n, empty), b["contributions"].get(n, empty)
            ),
        }
    rest = effects(a["rest"], b["rest"])
    chi = term(zero, zero_scatter, zero, b["chi"] - a["chi"])
    total = (lambda_a - lambda_b) * 1e5

    def summed(effect, ranges, select):
        return sum(
            v[effect][r]
            for n, v in nuclides.items()
            if select(n, v["kind"])
            for r in ranges
        )

    moved = ("to_dilute", "to_transported")
    parts = {
        "u238_resonance": summed(
            "absorption", ["ufg"], lambda n, k: k == "transported" and n == "U-238"
        ),
        "resonance_other": summed(
            "absorption", ["ufg"], lambda n, k: k == "transported" and n != "U-238"
        ),
        "absorption_other_ranges": summed(
            "absorption", ["fast", "thermal"], lambda n, k: k == "transported"
        ),
        "transported_scattering_production": summed(
            "scattering", RANGES, lambda n, k: k == "transported"
        )
        + summed("production", RANGES, lambda n, k: k == "transported"),
        "moved_absorption": summed("absorption", RANGES, lambda n, k: k in moved),
        "moved_scattering_fast": summed(
            "scattering", ["fast"], lambda n, k: k in moved
        ),
        "moved_scattering_below": summed(
            "scattering", ["ufg", "thermal"], lambda n, k: k in moved
        ),
        "moved_production": summed("production", RANGES, lambda n, k: k in moved),
        "rest": sum(summed(e, RANGES, lambda n, k: k == "dilute") for e in EFFECTS)
        + sum(rest[e][r] for e in EFFECTS for r in RANGES)
        + chi,
    }
    both = (
        "u238_resonance",
        "resonance_other",
        "absorption_other_ranges",
        "transported_scattering_production",
    )
    return {
        "total_pcm": total,
        "closure_pcm": sum(parts.values()) - total,
        "parts": parts,
        "transported_pcm": sum(parts[k] for k in both),
        "moved_absorption_ranges": {
            r: summed("absorption", [r], lambda n, k: k in moved) for r in RANGES
        },
        "moved": sorted(n for n, v in nuclides.items() if v["kind"] in moved),
        "nuclides": nuclides,
        "rest_effects": rest,
        "chi_pcm": chi,
    }


def pruned(nuclides):
    """すべての項の絶対値の最大が PRUNE_PCM 以上の核種と、外した核種の数と和"""
    kept, count, dropped = {}, 0, 0.0
    for n, v in nuclides.items():
        values = [v[e][r] for e in EFFECTS for r in RANGES]
        if max(abs(x) for x in values) >= PRUNE_PCM:
            kept[n] = v
        else:
            count += 1
            dropped += sum(values)
    return kept, {"count": count, "pcm": dropped}


def summarize(cases):
    """`cases`: (名前, case) の列"""
    steps = []
    for (name_a, a), (name_b, b) in itertools.pairwise(cases):
        s = step(a, b)
        kept, dropped = pruned(s.pop("nuclides"))
        steps.append(
            {"from": name_a, "to": name_b, **s, "nuclides": kept, "pruned": dropped}
        )
    return {
        "format": FORMAT,
        "cases": [
            {
                "name": name,
                "k_state": c["k_state"],
                "k_matrix": c["k_matrix"],
                "k_history": c["k_history"],
                "dilute_check": c["dilute_check"],
                "transported": len(c["transported"]),
                "dilute": len(c["dilute"]),
            }
            for name, c in cases
        ],
        "steps": steps,
    }


def report(summary):
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    for c in summary["cases"]:
        lines.append(
            f"{c['name']}: k（状態 {c['k_state']:.6f}、群定数 {c['k_matrix']:.6f}、履歴 {c['k_history']:.6f}）、"
            f"輸送 {c['transported']}、無限希釈 {c['dilute']}、無限希釈の吸収の冷却の断面積との相対差 {c['dilute_check']:+.1e}"
        )
    for s in summary["steps"]:
        parts = "、".join(f"{k} {v:+.1f}" for k, v in s["parts"].items())
        lines.append(
            f"{s['from']} → {s['to']}: 差 {s['total_pcm']:+.1f} pcm（和との差 {s['closure_pcm']:+.1e}）。{parts}"
        )
        ranges = s["moved_absorption_ranges"]
        lines.append(
            f"  両方で輸送する核種の和 {s['transported_pcm']:+.1f}、扱いを変えた核種の吸収の範囲ごと: "
            + "、".join(f"{r} {ranges[r]:+.1f}" for r in RANGES)
        )
        largest = sorted(
            s["nuclides"].items(),
            key=lambda kv: -abs(sum(kv[1][e][r] for e in EFFECTS for r in RANGES)),
        )[:8]
        lines.append(
            "  核種の和の大きい順: "
            + "、".join(
                f"{n}（{v['kind']}）{sum(v[e][r] for e in EFFECTS for r in RANGES):+.1f}"
                for n, v in largest
            )
        )
    total = sum(s["total_pcm"] for s in summary["steps"])
    first, last = summary["cases"][0], summary["cases"][-1]
    lines.append(
        f"段の和 {total:+.1f} pcm、最初と最後の履歴の k の差 "
        f"{(1.0 / first['k_history'] - 1.0 / last['k_history']) * 1e5:+.1f} pcm"
    )
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json")
    parser.add_argument("--report")
    parser.add_argument("--burnup", type=float)
    parser.add_argument("--case", nargs=4, action="append", default=[])
    args = parser.parse_args()
    if args.report:
        print("\n".join(report(json.loads(Path(args.report).read_text()))))
        return
    if not args.json or args.burnup is None or len(args.case) < 2:
        raise SystemExit(__doc__)
    cases = [
        (name, case(read(b), read(h), read(lib), args.burnup))
        for name, b, h, lib in args.case
    ]
    files = [f for _, *paths in args.case for f in paths]
    summary = {
        **summarize(cases),
        "burnup": args.burnup,
        "inputs": [
            {"path": f, "md5": hashlib.md5(Path(f).read_bytes()).hexdigest()}
            for f in dict.fromkeys(files)
        ],
    }
    Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
