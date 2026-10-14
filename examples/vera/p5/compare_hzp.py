"""VERA Problem 5（零出力試験）の xlc の結果を、仕様書の実測値と KENO-VI の値と比べる（ADR-0145、ADR-0146）。

使い方: python3 compare_hzp.py <spec.txt> <ケースの TOML> <core の出力>

ケースの状態は、仕様書の表 P5-2 の場合の番号で `case-<n>`（ホウ素濃度を与える）か `case-<n>-boron`（k = 1 になる
ホウ素濃度を探す。ADR-0145）と名付ける。道具は、状態のホウ素濃度、温度、B-10 の割合、バンクの位置（ステップ。
ADR-0127）が表の場合と同じことを確かめてから、core の出力で次を求める。反応度の差は仕様書と同じく ρ = (1/k₁ − 1/k₂) × 10⁵ pcm。

- 臨界ホウ素濃度（`case-<n>-boron`）: 表の場合のホウ素濃度（実測の臨界）との差
- k（`case-<n>`）: 臨界の場合 1〜10 は実測の臨界（k = 1）と、すべての場合は表 P5-6 の KENO-VI の k と、反応度で比べる
- バンクの価値: 場合 11（1170 ppm の ARO）と場合 12〜19 の差。合計も
- 微分ホウ素価値（DBW）: 場合 2 と 11 の差をホウ素濃度の差で割る（pcm/ppm）
- 等温温度係数（ITC）: 場合 20 と 21 の差を温度の差で割る（pcm/°F。表 P5-7 の単位）
- D バンクの積分価値: 場合 22（全挿入）からの各場合の反応度
- 出力分布（`case-1` の状態）: 付録 F の KENO-VI の軸方向の出力（49 層）と集合体の径方向の出力（1/8）。仕様書は付録 F の
  状態を明記しないが、本文の図 P5-8、P5-9 は D バンクを部分挿入した初臨界（場合 1）の形で、計装用シンブルを入れずに解いた
  もの（k は 1.000072）。軸方向は core の出力の z の面ごとの核分裂率の割合を線密度にし（燃料の高さの平均で 1）、参照の層
  ごとに厚さで重みを付けて平均する（compare_p3.py と同じ）。集合体は core の集合体ごとの核分裂率（193 体の平均で 1）を、
  付録 F の 1/8 の表を鏡映で広げた全炉心と比べる（compare_p5.py と同じ広げ方）

KENO-VI の k は、仕様書のとおり H-1 の S(α,β) の温度の補正（約 −43 pcm）を含み、熱膨張（約 −57 pcm）を含まない。
求めるのに要る状態がなければ、何が足りないかを出す。定義の取り違えを防ぐため、最初に表 P5-6 の KENO-VI の k から
表 P5-7 の KENO-VI の価値と DBW を計算し直し、一致しなければ止める（KENO-VI の ITC は付録 G の当てはめで求めた値
なので、2 点の差とは比べるだけにする）。
"""

import re
import sys
from pathlib import Path

import numpy as np
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_p3 import axial_report, coarsen
from compare_p5 import fuel_assemblies, summary, unfold

BANKS = ["A", "B", "C", "D", "SA", "SB", "SC", "SD"]
# 表 P5-2 の `-`（全引き抜き）のステップ
WITHDRAWN = 230
B10 = 0.199
NUMBER = r"\d+(?:\.\d+)?"
P56 = "Table P5-6: Problem 5 Reference Solution Eigenvalue Results"
P57 = "Table P5-7: Problem 5 Measured and Reference Solution ZPPT Results"
ROW = re.compile(
    rf"^\s*(\d+)\s+(\d+|↓)\s+(\d+|↓)\s+((?:(?:\d+|-)\s+){{8}})({NUMBER}) ± ({NUMBER})"
)


# 付録 F の見出しと、その中の表の見出し
AXIAL_TABLE = "Problem 5 Average Axial Power Results"
RADIAL_TABLE = "Problem 5 Assembly Radial Power Results"
RADIAL_ROWS = "Assembly Radial Powers"
# 付録 F の 1/8 の表の行（8〜15 行）ごとの集合体の数（H 列から。compare_p5.py の OCTANT_ROWS と同じ）
OCTANT_ROWS = [1, 2, 3, 4, 5, 6, 6, 4]
# 燃料の下端と上端（仕様書の高さ、cm。下部炉心板の上面から）
FUEL_BOTTOM = 11.951
FUEL_TOP = 377.711
# core の z の原点（下の水の外端）の、仕様書の高さ（p5_inputs.py の軸方向の層）
CORE_ORIGIN = -55.0
POWER_STATE = "case-1"


def table_lines(spec: list[str], caption: str) -> list[str]:
    """本文の表の見出し（目次の行でなく最後の出現）から後の行。"""
    found = [i for i, line in enumerate(spec) if caption in line]
    if not found:
        raise SystemExit(f"caption not found: {caption}")
    return spec[found[-1] + 1 :]


def keno_cases(spec: list[str]) -> dict[int, dict]:
    """表 P5-6: 場合ごとのホウ素濃度、温度、バンクの位置（ステップ）、KENO-VI の k と標準偏差。"""
    cases: dict[int, dict] = {}
    boron = temperature = None
    for line in table_lines(spec, P56):
        if P57 in line:
            break
        m = ROW.match(line)
        if not m:
            continue
        n = int(m.group(1))
        boron = boron if m.group(2) == "↓" else float(m.group(2))
        temperature = temperature if m.group(3) == "↓" else float(m.group(3))
        steps = [WITHDRAWN if s == "-" else int(s) for s in m.group(4).split()]
        cases[n] = {
            "boron": boron,
            "temperature": temperature,
            "banks": dict(zip(BANKS, steps)),
            "k": float(m.group(5)),
            "sigma": float(m.group(6)),
        }
    if sorted(cases) != list(range(1, 33)):
        raise SystemExit(f"{P56}: expected cases 1 to 32, found {sorted(cases)}")
    return cases


def zppt_results(spec: list[str]) -> dict[str, tuple[float, float]]:
    """表 P5-7: 量ごとの (実測値, KENO-VI の値)。"""
    results = {}
    patterns = {
        **{bank: rf"Bank {bank} Worth" for bank in BANKS},
        "total": r"Total Bank Worths",
        "DBW": r"DBW \(pcm/ppm\)",
        "ITC": r"ITC \(pcm/F\)",
    }
    lines = table_lines(spec, P57)[:40]
    for name, label in patterns.items():
        for line in lines:
            m = re.search(
                rf"{label}(?: \(pcm\))?\s+(-?{NUMBER})\s+(-?{NUMBER}) ±", line
            )
            if m:
                results[name] = (float(m.group(1)), float(m.group(2)))
                break
        else:
            raise SystemExit(f"{P57}: no row for {label}")
    return results


def rho(k_from: float, k_to: float) -> float:
    """k_from から k_to への反応度の変化 (pcm)。"""
    return (1.0 / k_from - 1.0 / k_to) * 1.0e5


def worths(k: dict[int, float], cases: dict[int, dict]) -> dict[str, float]:
    """場合ごとの k から、バンクの価値 (pcm)、DBW (pcm/ppm)、ITC (pcm/°F)。割る量は表の場合の値から取る。要る場合が
    なければその量を含めない。"""
    out = {}
    for i, bank in enumerate(BANKS):
        if 11 in k and 12 + i in k:
            out[bank] = -rho(k[11], k[12 + i])
    if all(bank in out for bank in BANKS):
        out["total"] = sum(out[bank] for bank in BANKS)
    if 2 in k and 11 in k:
        out["DBW"] = rho(k[11], k[2]) / (cases[2]["boron"] - cases[11]["boron"])
    if 20 in k and 21 in k:
        # 1 K = 1.8 °F
        out["ITC"] = (
            rho(k[20], k[21])
            / (cases[21]["temperature"] - cases[20]["temperature"])
            / 1.8
        )
    return out


def check_keno_definitions(
    cases: dict[int, dict], table: dict[str, tuple[float, float]]
) -> None:
    """表 P5-6 の KENO-VI の k から、表 P5-7 の KENO-VI の価値と DBW を再現することを確かめる。"""
    computed = worths({n: c["k"] for n, c in cases.items()}, cases)
    # k は 6 桁で丸めてあり（反応度で 0.1 pcm 以下）、価値は整数、DBW は小数 2 桁で丸めてある
    margins = {name: 1.0 for name in [*BANKS, "total"]} | {"DBW": 0.01}
    for name, margin in margins.items():
        keno = table[name][1]
        if abs(computed[name] - keno) > margin:
            raise SystemExit(
                f"{name}: KENO-VI k of {P56} give {computed[name]:.2f}, but {P57} has {keno}"
            )
    print(
        f"KENO-VI definitions: worths and DBW of {P57} reproduced from {P56}; "
        f"ITC from cases 20 and 21 is {computed['ITC']:.2f} pcm/F "
        f"(the table's {table['ITC'][1]} is fitted over 550 to 580 K, Appendix G)"
    )


def keno_axial(spec: list[str]) -> tuple[list[float], list[float]]:
    """付録 F の 49 層: 燃料の下端からの境界（下から）と出力"""
    start = next((i for i, line in enumerate(spec) if AXIAL_TABLE in line), None)
    if start is None:
        raise SystemExit(f"{AXIAL_TABLE} not found")
    rows = []
    for line in spec[start : start + 60]:
        m = re.match(
            rf"\s*(\d+)\s+({NUMBER})\s+({NUMBER})\s+({NUMBER}) \+/- ({NUMBER})%", line
        )
        if m:
            rows.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
    if len(rows) != 49:
        raise SystemExit(f"{AXIAL_TABLE}: expected 49 levels, found {len(rows)}")
    rows.sort()
    bounds = [rows[0][0] - rows[0][1] - FUEL_BOTTOM]
    for top, thickness, _ in rows:
        if abs(top - thickness - FUEL_BOTTOM - bounds[-1]) > 1e-3:
            raise SystemExit(
                f"{AXIAL_TABLE}: level ending at {top} does not follow the one below"
            )
        bounds.append(top - FUEL_BOTTOM)
    if abs(bounds[0]) > 1e-3 or abs(bounds[-1] - (FUEL_TOP - FUEL_BOTTOM)) > 1e-3:
        raise SystemExit(
            f"{AXIAL_TABLE}: levels span {bounds[0]}-{bounds[-1]} cm, not the fuel height"
        )
    return bounds, [r[2] for r in rows]


def keno_radial(spec: list[str]) -> np.ndarray:
    """付録 F の 1/8 の表（行 8〜15、H 列から）を全炉心の 15 × 15 に広げる"""
    caption = next((i for i, line in enumerate(spec) if RADIAL_TABLE in line), None)
    if caption is None:
        raise SystemExit(f"{RADIAL_TABLE} not found")
    start = next(
        (i for i in range(caption, len(spec)) if spec[i].strip() == RADIAL_ROWS), None
    )
    if start is None:
        raise SystemExit(f"{RADIAL_ROWS} not found after {RADIAL_TABLE}")
    quarter = np.full((8, 8), np.nan)
    for r, line in enumerate(spec[start + 2 : start + 10]):
        values = [float(x) for x in re.findall(r"(?<![\d.])(\d\.\d{5})(?![\d%])", line)]
        if len(values) != OCTANT_ROWS[r]:
            raise SystemExit(
                f"{RADIAL_ROWS}: row {8 + r} has {len(values)} values, expected {OCTANT_ROWS[r]}"
            )
        quarter[r, : len(values)] = values
    # 1/8 の表の欠けた位置（下三角の外）を、対角の鏡映で埋める
    quarter = np.where(np.isnan(quarter), quarter.T, quarter)
    if np.count_nonzero(~np.isnan(quarter)) != 56:
        raise SystemExit(f"{RADIAL_ROWS}: expected 56 assemblies in the quarter")
    return unfold(quarter)


def state_power(path: Path, name: str) -> tuple[list[float], list[float], np.ndarray]:
    """core の出力の状態 `name` の、z の面の境界（燃料の下端から）と燃料の高さの平均で 1 の線密度の出力、集合体ごとの
    核分裂率（行は上から。反射体と `-` の位置は NaN）"""
    lines = path.read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if line == f"state {name}"), None)
    if start is None:
        raise SystemExit(f"{path}: no state {name}")
    end = next(
        (i for i in range(start + 1, len(lines)) if lines[i].startswith("state ")),
        len(lines),
    )
    block = lines[start:end]
    heads = {}
    for head in ("axial fission fractions", "assembly fission rates"):
        heads[head] = next(
            (i for i, line in enumerate(block) if line.startswith(head)), None
        )
        if heads[head] is None:
            raise SystemExit(f"{path}: state {name} has no `{head}` table")
    a, b = heads.values()
    nodes = [tuple(float(x) for x in line.split()) for line in block[a + 1 : b]]
    bottom, top = FUEL_BOTTOM - CORE_ORIGIN, FUEL_TOP - CORE_ORIGIN
    tolerance = 1e-4
    fuel = [
        (z0, z1, f)
        for z0, z1, f in nodes
        if z0 >= bottom - tolerance and z1 <= top + tolerance
    ]
    outside = sum(f for z0, z1, f in nodes) - sum(f for *_, f in fuel)
    if (
        not fuel
        or abs(fuel[0][0] - bottom) > tolerance
        or abs(fuel[-1][1] - top) > tolerance
        or abs(outside) > 1e-12
    ):
        raise SystemExit(
            f"{path}: state {name}: fission outside the fuel height or no node bound at its ends"
        )
    total = sum(f for *_, f in fuel)
    height = top - bottom
    bounds = [fuel[0][0] - bottom] + [z1 - bottom for _, z1, _ in fuel]
    axial = [f / (z1 - z0) / (total / height) for z0, z1, f in fuel]
    rows = [
        [np.nan if x in (".", "-") else float(x) for x in line.split()]
        for line in block[b + 2 :]
        if line.strip()
    ]
    return bounds, axial, np.array(rows)


def power_report(spec: list[str], core_out: Path) -> None:
    print(
        f"\npower distribution ({POWER_STATE} against Appendix F, KENO-VI without the instrument thimbles)"
    )
    bounds, ours, rates = state_power(core_out, POWER_STATE)
    keno_bounds, keno = keno_axial(spec)
    axial_report("KENO-VI", keno_bounds, coarsen(bounds, ours, keno_bounds), keno)
    rates = fuel_assemblies(rates)
    reference = keno_radial(spec)
    if rates.shape != reference.shape or (np.isnan(rates) != np.isnan(reference)).any():
        raise SystemExit(
            "assemblies: the xlc and KENO-VI cores have different fuel positions"
        )
    fuel = ~np.isnan(reference)
    rates = rates / rates[fuel].mean()
    summary(
        "assembly relative difference vs KENO-VI",
        (rates[fuel] / reference[fuel] - 1),
    )
    print(
        "octant (rows 8-15 from column H): xlc (mean of the 8 symmetric positions), KENO-VI, relative difference"
    )
    images = [
        rates,
        rates.T,
        rates[::-1],
        rates[:, ::-1],
        rates[::-1, ::-1],
        rates.T[::-1],
        rates.T[:, ::-1],
        rates.T[::-1, ::-1],
    ]
    folded = np.mean(images, axis=0)
    print(
        f"largest spread of the 8 symmetric positions in xlc: {100 * np.nanmax(np.max(images, axis=0) / np.min(images, axis=0) - 1):.3f}%"
    )
    for r, n in enumerate(OCTANT_ROWS):
        cells = [
            f"{folded[7 + r, 7 + c]:.4f} {reference[7 + r, 7 + c]:.4f} {100 * (folded[7 + r, 7 + c] / reference[7 + r, 7 + c] - 1):+6.2f}%"
            for c in range(n)
        ]
        print(f"  {8 + r:2d}  " + "   ".join(cells))


def core_results(path: Path) -> dict[str, dict]:
    """core の出力の、状態ごとの k と探索したホウ素濃度。"""
    results: dict[str, dict] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith("state "):
            name = line.split()[1]
            if name in results:
                raise SystemExit(f"{path}: state {name} appears twice")
            results[name] = {}
        elif line.startswith("boron "):
            if name is None or "boron" in results[name]:
                raise SystemExit(f"{path}: boron line outside a state block or twice")
            results[name]["boron"] = float(line.split()[1])
        elif line.startswith("k-eff "):
            if name is None:
                raise SystemExit(f"{path}: k-eff before any `state <name>` line")
            if "k" in results[name]:
                raise SystemExit(f"{path}: state {name} has two k-eff lines")
            results[name]["k"] = float(line.split()[1])
    if not results or any("k" not in r for r in results.values()):
        raise SystemExit(f"{path}: expected `state <name>` blocks with k-eff")
    return results


def case_states(path: Path, cases: dict[int, dict]) -> dict[str, tuple[int, bool]]:
    """ケースの状態の名前ごとに (場合の番号, 探索か)。状態が表の場合と違えば止める。"""
    states = tomllib.loads(path.read_text()).get("state", {})
    out = {}
    for name, state in states.items():
        if not isinstance(state, dict):
            raise SystemExit(f"{path}: state {name}: name each state as [state.<name>]")
        m = re.fullmatch(r"case-([1-9]\d*)(-boron)?", name)
        if not m or int(m.group(1)) not in cases:
            raise SystemExit(
                f"{path}: state {name}: name it case-<n> or case-<n>-boron (Table P5-2)"
            )
        n, search = int(m.group(1)), m.group(2) is not None
        case = cases[n]
        if state.get("banks") != case["banks"]:
            raise SystemExit(
                f"{path}: state {name}: banks {state.get('banks')} differ from case {n} ({case['banks']})"
            )
        if search and n > 10:
            raise SystemExit(
                f"{path}: state {name}: only the critical cases 1 to 10 search the boron"
            )
        target = state.get("search")
        if search != (target is not None) or (
            search and (not isinstance(target, dict) or target.get("k") != 1.0)
        ):
            raise SystemExit(
                f"{path}: state {name}: -boron states search the boron for k = 1, others give it"
            )
        wanted = {
            "b10": B10,
            "moderator": case["temperature"],
            "fuel": case["temperature"],
        }
        if not search:
            wanted["boron"] = case["boron"]
        for key, value in wanted.items():
            if state.get(key) != value:
                raise SystemExit(
                    f"{path}: state {name}: {key} {state.get(key)} differs from case {n} ({value})"
                )
        out[name] = (n, search)
    return out


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    spec = Path(sys.argv[1]).read_text().splitlines()
    cases = keno_cases(spec)
    table = zppt_results(spec)
    check_keno_definitions(cases, table)
    states = case_states(Path(sys.argv[2]), cases)
    results = core_results(Path(sys.argv[3]))
    if set(results) != set(states):
        raise SystemExit(
            f"states of the case {sorted(states)} differ from the output {sorted(results)}"
        )

    searches = {
        name: state["search"]
        for name, state in tomllib.loads(Path(sys.argv[2]).read_text())["state"].items()
        if "search" in state
    }
    for name in states:
        if (name in searches) != ("boron" in results[name]):
            raise SystemExit(
                f"state {name}: the output has a boron line only for -boron states"
            )
        # 探索の状態の k は、ケースに書いた許容の中で 1（ADR-0145）。出力の k は小数 8 桁に丸めてある
        rounding = 0.5e-8
        if (
            name in searches
            and abs(results[name]["k"] - 1.0) > searches[name]["tolerance"] + rounding
        ):
            raise SystemExit(
                f"state {name}: searched k {results[name]['k']} is not within the tolerance of 1"
            )

    print("\ncritical boron (ppm)        xlc  measured   diff")
    if not any(search for _, search in states.values()):
        print("not computed (needs case-<n>-boron states)")
    for name, (n, search) in sorted(states.items(), key=lambda s: s[1]):
        if search:
            ppm = results[name]["boron"]
            measured = cases[n]["boron"]
            print(
                f"case {n:2d} {'':14s} {ppm:9.2f} {measured:9.2f} {ppm - measured:+6.2f}"
            )

    print(
        "\nk by case                      xlc    KENO-VI   xlc-KENO (pcm)  xlc-critical (pcm)"
    )
    k = {}
    for name, (n, search) in sorted(states.items(), key=lambda s: s[1]):
        if search:
            continue
        k[n] = results[name]["k"]
        keno = cases[n]["k"]
        critical = f"{rho(1.0, k[n]):+8.1f}" if n <= 10 else ""
        # 反応度 1/k の標準偏差は σ_k / k²
        sigma = cases[n]["sigma"] / keno**2 * 1e5
        print(
            f"case {n:2d} {'':12s} {k[n]:.8f}  {keno:.6f}  {rho(keno, k[n]):+8.1f} ± {sigma:.1f}  {critical}"
        )

    ours = worths(k, cases)
    print(
        "\nquantity                    xlc  measured   KENO-VI  vs measured  vs KENO-VI"
    )
    for name in [*BANKS, "total", "DBW", "ITC"]:
        label = {
            "DBW": "DBW (pcm/ppm)",
            "ITC": "ITC (pcm/F)",
            "total": "total bank worth (pcm)",
        }.get(name, f"bank {name} worth (pcm)")
        if name not in ours:
            print(f"{label:24s} not computed (needs its cases as states)")
            continue
        measured, keno_table = table[name]
        # 表 P5-7 の KENO-VI の値（ITC は付録 G の当てはめ）。表と同じく、価値は比、DBW と ITC は差で比べる
        if name in ("DBW", "ITC"):
            vs = [f"{ours[name] - x:+8.2f}" for x in (measured, keno_table)]
        else:
            vs = [f"{ours[name] / x - 1:+8.1%}" for x in (measured, keno_table)]
        print(
            f"{label:24s} {ours[name]:8.2f} {measured:9.2f} {keno_table:9.2f}"
            f"  {vs[0]}        {vs[1]}"
        )

    missing = [n for n in range(22, 33) if n not in k]
    print("\nbank D integral worth (pcm from case 22)")
    if missing:
        print(f"not computed (needs cases {missing} as states)")
    else:
        for n in range(23, 33):
            steps = cases[n]["banks"]["D"]
            keno = rho(cases[22]["k"], cases[n]["k"])
            print(
                f"D at {steps:3d} steps  {rho(k[22], k[n]):8.1f}  KENO-VI {keno:8.1f}"
            )

    if POWER_STATE in results:
        power_report(spec, Path(sys.argv[3]))
    else:
        print(f"\npower distribution: not computed (needs the state {POWER_STATE})")


if __name__ == "__main__":
    main()
