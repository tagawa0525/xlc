"""VERA Problem 5（零出力試験）の xlc の結果を、仕様書の実測値と KENO-VI の値と比べる（ADR-0145、ADR-0146）。

使い方: python3 compare_hzp.py <spec.txt> <ケースの TOML> <core の出力>

ケースの状態は、仕様書の表 P5-2 の場合の番号で `case-<n>`（ホウ素濃度を与える）か `case-<n>-boron`（k = 1 になる
ホウ素濃度を探す。ADR-0145）と名付ける。道具は、状態のホウ素濃度、温度、B-10 の割合が表の場合と同じことを確かめて
から、core の出力で次を求める。反応度の差は仕様書と同じく ρ = (1/k₁ − 1/k₂) × 10⁵ pcm。

- 臨界ホウ素濃度（`case-<n>-boron`）: 表の場合のホウ素濃度（実測の臨界）との差
- k（`case-<n>`）: 臨界の場合 1〜10 は実測の臨界（k = 1）と、すべての場合は表 P5-6 の KENO-VI の k と、反応度で比べる
- バンクの価値: 場合 11（1170 ppm の ARO）と場合 12〜19 の差。合計も
- 微分ホウ素価値（DBW）: 場合 2 と 11 の差をホウ素濃度の差で割る（pcm/ppm）
- 等温温度係数（ITC）: 場合 20 と 21 の差を温度の差で割る（pcm/°F。表 P5-7 の単位）
- D バンクの積分価値: 場合 22（全挿入）からの各場合の反応度

KENO-VI の k は、仕様書のとおり H-1 の S(α,β) の温度の補正（約 −43 pcm）を含み、熱膨張（約 −57 pcm）を含まない。
求めるのに要る状態がなければ、何が足りないかを出す。定義の取り違えを防ぐため、最初に表 P5-6 の KENO-VI の k から
表 P5-7 の KENO-VI の価値と DBW を計算し直し、一致しなければ止める（KENO-VI の ITC は付録 G の当てはめで求めた値
なので、2 点の差とは比べるだけにする）。
"""

import re
import sys
from pathlib import Path

import tomllib

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
        if any(steps != WITHDRAWN for steps in case["banks"].values()):
            # バンクの位置を状態に書く形（ADR-0127、#0063）が入るまで、挿入した場合は状態で表せない
            raise SystemExit(
                f"{path}: state {name}: case {n} inserts banks, which states cannot hold yet (#0063)"
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

    print("\ncritical boron (ppm)        xlc  measured   diff")
    for name, (n, search) in sorted(states.items(), key=lambda s: s[1]):
        if search:
            if "boron" not in results[name]:
                raise SystemExit(
                    f"state {name} searches the boron, but the output has no boron line"
                )
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


if __name__ == "__main__":
    main()
