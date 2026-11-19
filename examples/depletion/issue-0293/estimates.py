"""#0293 の項目 6、7、3 の大きさを、xlc の結果（urr/out/<燃料>.msgpack の 69 群のテーブルと micro）から 1 次の摂動で見積もる。

中性子束は、テーブルの 69 群の定数の無限媒質の固有値問題（(diag Σt − Sᵀ) φ = χ (νΣf·φ) / k）を解いて求める（鏡面反射の格子の
体積で平均した定数なので、k∞ と群の中性子束の比を保つ。解いた k と表の k∞ の差を書いて確かめる）。核種ごとの反応率は均質化した
数密度 N̄ と微視的断面積と φ の積。

- 項目 6（非分離共鳴の自己遮蔽）: UFG の範囲（4 eV〜15.03 keV。UO2 と MOX の燃料の上端は U-238 の 20 keV 以下で最も高い群の境界。
  ADR-0069）の群のうち、核種の非分離共鳴の範囲（ENDF の MF2 の LRU = 2 の EL〜EH）に入る群で、UFG は無限希釈の平均の断面積を使う。
  その群の自己遮蔽因子 f = σ(σ0) / σ(∞) をライブラリの σ0 の表（ラダーのサンプリングで作った非分離共鳴の自己遮蔽。ADR-0177）から
  取り、吸収と生成をそれぞれ (1 − f) 倍だけ減らしたときの k の変化を出す。σ0 は燃料の領域の組成（steps/out/<燃料>-40.msgpack の
  燃焼点の数密度）の均質な Bondarenko の背景 Σ_{j≠i} N_j σt_j(∞) / N_i で、燃料の外への逃げの項を入れない（σ0 が小さくなり、遮蔽を
  大きく見積もる）。被覆管にもある核種（Sn）は、燃料と被覆管の組成の背景の小さい方を使う。境界の群は、群の下端が EL 以上の群だけを数える
- 項目 7（収率の 2 点の重み）: 核分裂率のうち収率の境界（112.47 eV。ADR-0183）より上の割合（`thermal_fission` から）と、その
  核分裂の収率を 500 keV の点から熱の点に置き換えたときの、主な核分裂生成物の累積の収率の相対変化 δ と、吸収の割合 a から
  Δk/k ≈ −a δ（数密度が収率に比例するとみなす上限。熱外の核分裂の本当の収率は 2 点の間にあると見込む）
- 項目 3（出力の規格化）: 核分裂の Q（チェーンの値）の和に対する、捕獲（MT102）の Q の和の割合。micro の核種は核種ごとの
  MT102 の Q（ENDF の MF3 MT102 の QI）で数え、micro にない吸収（減速材の H-1、B-10、被覆管、micro に挙げない核種）は吸収の残りに
  Q の範囲（REST_Q）を掛けて幅で書く

使い方（python-tests の devShell。endf と numpy を使う）: python3 estimates.py <燃料> [--json <出力>]
"""

import argparse
import json
import math
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib

HERE = Path(__file__).resolve().parent
DEPLETION = HERE.parent

# UFG の範囲（ADR-0069。WIMS の 69 群の境界）
UFG_LOW = 4.0
UFG_HIGH = 15030.0
# 収率の境界（ADR-0183。√(0.0253 × 5e5) eV）
YIELD_CUTOFF = math.sqrt(0.0253 * 5e5)
FUEL_TEMPERATURE = 900.0
CLAD_TEMPERATURE = 565.0
# micro にない吸収の捕獲の Q の範囲 (eV): H-1(n,γ) 2.22 MeV から、被覆管と核分裂生成物の (n,γ) の典型の 8 MeV まで
REST_Q = (2.22e6, 8.0e6)


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def flux(state: dict) -> tuple[np.ndarray, float]:
    """少数群定数 `state` の無限媒質の中性子束（νΣf·φ = k になる規格化）と k"""
    total = np.array(state["total"])
    scatter = np.array(state["scatter"])
    matrix = np.diag(total) - scatter.T
    phi = np.linalg.solve(matrix, np.array(state["chi"]))
    return phi, float(np.dot(state["nu_fission"], phi))


def shielding_factor(table: list[float], sigma0: list[float], value: float) -> float:
    """σ0 の表 `table`（`sigma0` の順。最初が無限希釈）を ln σ0 に線形に内挿した値の、無限希釈の値に対する比"""
    if table[0] == 0.0:
        return 1.0
    if value >= sigma0[0]:
        return 1.0
    for i in range(len(sigma0) - 1):
        high, low = sigma0[i], sigma0[i + 1]
        if low <= value <= high:
            t = (math.log(high) - math.log(value)) / (math.log(high) - math.log(low))
            return (table[i] + t * (table[i + 1] - table[i])) / table[0]
    return table[-1] / table[0]


def urr_groups(bounds: list[float], urr: tuple[float, float]) -> list[int]:
    """UFG の範囲の中で、下端が非分離共鳴の範囲の下端 EL 以上の群（降順の境界 `bounds` の群の番号）"""
    low, high = urr
    groups = []
    for g in range(len(bounds) - 1):
        top, bottom = bounds[g], bounds[g + 1]
        if bottom >= low and top <= UFG_HIGH and bottom >= UFG_LOW and bottom < high:
            groups.append(g)
    return groups


def cumulative_yields(chain: dict, independent: dict[int, float]) -> dict[int, float]:
    """独立収率（核種の番号ごと）から、崩壊の分岐を辿った累積の収率"""
    nuclides = chain["nuclides"]
    parents: dict[int, list[tuple[int, float]]] = {}
    for p, n in enumerate(nuclides):
        for mode in (n.get("decay") or {}).get("modes", []):
            for product in mode["products"]:
                parents.setdefault(product["nuclide"], []).append(
                    (p, mode["ratio"] * product["count"])
                )
    memo: dict[int, float] = {}

    def cumulative(i: int, path: frozenset) -> float:
        if i in memo:
            return memo[i]
        value = independent.get(i, 0.0)
        for p, ratio in parents.get(i, []):
            if p not in path:
                value += ratio * cumulative(p, path | {i})
        memo[i] = value
        return value

    return {i: cumulative(i, frozenset()) for i in range(len(nuclides))}


def yield_points(fission: dict) -> tuple[int, int] | None:
    """熱（0.0253 eV）と高速（500 keV）の点の番号。2 点が無い核種は None（収率の重みが点を 1 つしか使わない）"""
    energies = fission["energies"]
    try:
        return energies.index(0.0253), energies.index(5e5)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fuel", choices=["u", "mox"])
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    import endf

    table = read(HERE / "urr" / "out" / f"{args.fuel}.msgpack")
    history = read(HERE / "steps" / "out" / f"{args.fuel}-40.msgpack")
    library = read(DEPLETION / "depletion.msgpack")
    chain = read(DEPLETION / "chain.msgpack")
    xs = tomllib.loads((DEPLETION / "xs.toml").read_text())
    endf_directory = DEPLETION / xs["endf-directory"]

    bounds = table["bounds"]
    sigma0 = library["sigma0"]
    by_name = {n["name"]: n for n in library["nuclides"]}
    chain_index = {n["name"]: i for i, n in enumerate(chain["nuclides"])}

    def at_temperature(name: str, temperature: float = FUEL_TEMPERATURE) -> dict:
        for t in by_name[name]["temperatures"]:
            if t["temperature"] == temperature:
                return t
        sys.exit(f"{name}: no {temperature} K in the library")

    # 核種ごとの ENDF の非分離共鳴の範囲と MT102 の Q
    micro = table["micro"]
    names = micro["nuclides"]
    urr: dict[str, tuple[float, float]] = {}
    capture_q: dict[str, float] = {}
    for name in names:
        material = endf.Material(endf_directory / xs["nuclide"][name]["endf"])
        for iso in material.section_data.get((2, 151), {"isotopes": []})["isotopes"]:
            for r in iso["ranges"]:
                if r["LRU"] == 2:
                    urr[name] = (r["EL"], r["EH"])
        if (3, 102) not in material.section_data:
            sys.exit(f"{name}: no MF3 MT102")
        capture_q[name] = material.section_data[3, 102]["QI"]
    missing = [
        n
        for n in by_name
        if n not in names
        and n in history["transportable"]
        and _has_urr_in_ufg(endf, endf_directory, xs, n)
    ]
    if missing:
        sys.exit(
            f"nuclides with unresolved resonances in the UFG range are not in micro: {missing}"
        )

    # 被覆管の組成（Sn の非分離共鳴の背景。被覆管の Sn は micro の均質化した反応率に含まれる）
    branch = tomllib.loads((HERE / "urr" / f"{args.fuel}.toml").read_text())
    clad = branch["material"]["zircaloy"]["nuclides"]
    burnup_axis = [a["name"] for a in table["axes"]].index("burnup")
    burnups = table["axes"][burnup_axis]["points"]
    fuel_region = history["regions"][0]
    if fuel_region["material"] != "fuel" or len(history["regions"]) != 1:
        sys.exit("expected one depletion region `fuel`")
    results = []
    for point, state, nuclides in zip(
        table["points"], table["states"], micro["states"]
    ):
        b = burnups[point[burnup_axis]]
        densities = micro["densities"][point[burnup_axis]]
        phi, k = flux(state)
        absorption = float(np.dot(state["absorption"], phi))
        production = float(np.dot(state["nu_fission"], phi))

        # 燃料の領域の組成（σ0 の背景）
        hp = next(p for p in history["points"] if p["burnup"] == b)
        fuel = {
            n: d
            for n, d in zip(history["nuclides"], hp["densities"][0])
            if d > 0 and n in by_name and n in history["transportable"]
        }
        materials = [fuel, clad]
        # 背景の全断面積は物質の温度で（燃料 900 K、被覆管 565 K）
        totals = {}
        for m, temperature in ((clad, CLAD_TEMPERATURE), (fuel, FUEL_TEMPERATURE)):
            for n in m:
                totals[n] = np.array(
                    at_temperature(n, temperature)["shielded"][0]["total"]
                )

        # 項目 6
        d_absorption = d_production = 0.0
        urr_rows = {}
        for i, name in enumerate(names):
            if name not in urr or not any(m.get(name, 0.0) > 0 for m in materials):
                continue
            groups = urr_groups(bounds, urr[name])
            if not groups:
                continue
            shielded = at_temperature(name)["shielded"]
            s = nuclides[i]
            da = dp = 0.0
            for g in groups:
                # 燃料と被覆管（Sn）のうち、σ0 の小さい（遮蔽の強い）方
                background = min(
                    sum(m[j] * totals[j][g] for j in m if j != name) / m[name]
                    for m in materials
                    if m.get(name, 0.0) > 0
                )
                absorbed = [t["capture"][g] + t["fission"][g] for t in shielded]
                f_a = shielding_factor(absorbed, sigma0, background)
                f_nf = shielding_factor(
                    [t["nu_fission"][g] for t in shielded], sigma0, background
                )
                sigma_a = s["total"][g] - sum(s["scatter"][g])
                da += densities[i] * phi[g] * sigma_a * (1 - f_a)
                dp += densities[i] * phi[g] * s["nu_fission"][g] * (1 - f_nf)
            d_absorption += da
            d_production += dp
            urr_rows[name] = {
                "groups": [bounds[groups[-1] + 1], bounds[groups[0]]],
                "dk_pcm": (
                    (production - dp) / (absorption - da) / (production / absorption)
                    - 1
                )
                * 1e5,
            }
        urr_pcm = (
            (production - d_production)
            / (absorption - d_absorption)
            / (production / absorption)
            - 1
        ) * 1e5

        # 項目 7
        fissions = {}
        thermal = {}
        for i, name in enumerate(names):
            s = nuclides[i]
            if s["thermal_fission"] is None:
                continue
            fissions[name] = densities[i] * float(np.dot(s["fission"], phi))
            thermal[name] = s["thermal_fission"]
        total_fission = sum(fissions.values())
        epithermal = (
            sum(f * (1 - thermal[n]) for n, f in fissions.items()) / total_fission
        )
        cumulative = {}
        for name in fissions:
            fission = chain["nuclides"][chain_index[name]]["fission"]
            points = yield_points(fission) if fission else None
            if points is None:
                continue
            cumulative[name] = [
                cumulative_yields(
                    chain, {y["nuclide"]: y["value"] for y in fission["yields"][p]}
                )
                for p in points
            ]
        yield_rows = {}
        for i, name in enumerate(names):
            if name in fissions or name not in chain_index:
                continue
            j = chain_index[name]
            current = alternative = 0.0
            for f_name, f in fissions.items():
                if f_name not in cumulative:
                    continue
                th, fast = cumulative[f_name][0][j], cumulative[f_name][1][j]
                current += f * (thermal[f_name] * th + (1 - thermal[f_name]) * fast)
                alternative += f * th
            if current == 0.0:
                continue
            s = nuclides[i]
            share = (
                densities[i]
                * sum(
                    (s["total"][g] - sum(s["scatter"][g])) * phi[g]
                    for g in range(len(phi))
                )
                / absorption
            )
            delta = alternative / current - 1
            yield_rows[name] = {
                "share": share,
                "delta": delta,
                "dk_pcm": -share * delta * 1e5,
            }

        # 項目 3
        q_fission = sum(
            f * chain["nuclides"][chain_index[n]]["fission"]["q"]
            for n, f in fissions.items()
        )
        q_capture = 0.0
        micro_absorption = 0.0
        for i, name in enumerate(names):
            s = nuclides[i]
            micro_absorption += densities[i] * sum(
                (s["total"][g] - sum(s["scatter"][g])) * phi[g] for g in range(len(phi))
            )
            capture = sum(
                densities[i] * float(np.dot(channel, phi))
                for reaction in s["reactions"]
                if reaction["mt"] == 102
                for channel in reaction["channels"]
            )
            q_capture += capture * capture_q[name]
        rest = absorption - micro_absorption
        energy = [(q_capture + rest * q) / q_fission for q in REST_Q]

        results.append(
            {
                "burnup": b,
                "k_table": state["k_inf"],
                "k_flux": k,
                "urr": {"dk_pcm": urr_pcm, "nuclides": urr_rows},
                "yields": {
                    "epithermal_fission": epithermal,
                    "dk_pcm": sum(r["dk_pcm"] for r in yield_rows.values()),
                    "abs_dk_pcm": sum(abs(r["dk_pcm"]) for r in yield_rows.values()),
                    "nuclides": yield_rows,
                },
                "energy": {
                    "capture_micro": q_capture / q_fission,
                    "rest_absorption_share": rest / absorption,
                    "capture_over_fission": energy,
                },
            }
        )
    add_burnup_shift(results, history)
    for r in results:
        print(
            f"{r['burnup']:g} MWd/kg: k {r['k_table']:.5f} (flux {r['k_flux']:.5f}); "
            f"URR {r['urr']['dk_pcm']:+.1f} pcm; epithermal fission {r['yields']['epithermal_fission'] * 100:.1f}%, "
            f"yields {r['yields']['dk_pcm']:+.1f} pcm (|sum| {r['yields']['abs_dk_pcm']:.1f}); "
            f"capture energy / fission energy {r['energy']['capture_over_fission'][0] * 100:.2f}–"
            f"{r['energy']['capture_over_fission'][1] * 100:.2f}%"
        )
        if r["energy"]["burnup_shift"]:
            lo, hi = r["energy"]["burnup_shift"]
            print(
                f"  normalization: mean delta {lo['mean_delta'] * 100:.2f}–{hi['mean_delta'] * 100:.2f}%, "
                f"burnup {lo['burnup']:.2f}–{hi['burnup']:.2f} MWd/kg, k-inf {lo['dk_pcm']:+.0f} to {hi['dk_pcm']:+.0f} pcm"
            )
        top = sorted(r["urr"]["nuclides"].items(), key=lambda x: x[1]["dk_pcm"])[:4]
        print("  URR: " + ", ".join(f"{n} {v['dk_pcm']:+.1f}" for n, v in top))
        top = sorted(
            r["yields"]["nuclides"].items(), key=lambda x: -abs(x[1]["dk_pcm"])
        )[:4]
        print(
            "  yields: "
            + ", ".join(
                f"{n} {v['delta'] * 100:+.2f}% ({v['dk_pcm']:+.1f})" for n, v in top
            )
        )
    if args.json:
        args.json.write_text(
            json.dumps({"fuel": args.fuel, "results": results}, indent=1) + "\n"
        )


def add_burnup_shift(results: list[dict], history: dict) -> None:
    """項目 3 の k∞ への効きの見込み: 出力を捕獲の分も含むエネルギーで規格化すると、同じ燃焼度までの核分裂が 1/(1 + δ̄) に減る
    （δ̄ は燃焼度 0 からの δ の燃焼度の平均。燃焼点の間は台形）。今の規格化の燃焼度 B の組成は、その規格化の B (1 + δ̄) の組成に
    当たるので、同じ燃焼度の k∞ の差（捕獲を含む規格化 − 今の規格化）は −(dk/dB) B δ̄ / k（dk/dB は履歴の 1 つ前の燃焼点との差分）。δ は REST_Q の 2 つの端ごと"""
    ks = {p["burnup"]: p["k_inf"] for p in history["points"]}
    points = sorted(ks)
    burnups = [r["burnup"] for r in results]
    for i, r in enumerate(results):
        b = r["burnup"]
        if b == 0:
            r["energy"]["burnup_shift"] = None
            continue
        j = points.index(b)
        slope = (ks[b] - ks[points[j - 1]]) / (b - points[j - 1])
        shift = []
        for end in range(2):
            deltas = [
                results[m]["energy"]["capture_over_fission"][end] for m in range(i + 1)
            ]
            area = sum(
                (deltas[m] + deltas[m + 1]) / 2 * (burnups[m + 1] - burnups[m])
                for m in range(i)
            )
            mean = area / b
            shift.append(
                {
                    "mean_delta": mean,
                    "burnup": b * mean,
                    "dk_pcm": -slope * b * mean / ks[b] * 1e5,
                }
            )
        r["energy"]["burnup_shift"] = shift


def _has_urr_in_ufg(endf, directory: Path, xs: dict, name: str) -> bool:
    material = endf.Material(directory / xs["nuclide"][name]["endf"])
    if (2, 151) not in material.section_data:
        return False
    for iso in material.section_data[2, 151]["isotopes"]:
        for r in iso["ranges"]:
            if r["LRU"] == 2 and r["EL"] < UFG_HIGH:
                return True
    return False


if __name__ == "__main__":
    main()
