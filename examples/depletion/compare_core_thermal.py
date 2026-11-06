"""熱水力のある core の燃焼（ADR-0229 の燃焼のステップ）を、同じサイクルを入口の温度の一様な状態（熱水力なし）で解いた燃焼と比べ、
点と中点の k とホウ素の差、熱水力の効きと、テーブルの k∞ の係数とノードの温度から見込んだ効き、固有値計算の数、外側の反復の回数、
時間、熱水力の繰り返しの回数を出す。

使い方: python3 compare_core_thermal.py --json <要約の出力先> --table <2A のテーブル（msgpack）> <名前>=<熱水力なしの履歴>,<そのログ>,<熱水力の履歴>,<そのログ> ...
        python3 compare_core_thermal.py --report <要約の JSON>

- 履歴とログは compare_core_guess.py と同じ（差は熱水力 − 熱水力なし）。2 つの燃焼は、同じ燃焼点（EFPD と日数）、状態の名前、定格の出力
  密度、軸方向のメッシュ、集合体を持ち、熱水力なしの状態の減速材と燃料の温度は、熱水力の状態の入口の温度と同じこと
- 熱水力の効き `effect` (ppm): 点のホウ素の差から、k の差の反応度 (1 − 1/k) を k∞ のホウ素の係数でホウ素に直した分を引いた値（探索した点は
  k が同じでホウ素の差、ホウ素を与えた点は k の差）
- 見込み `estimate` (ppm): テーブルの k∞ の反応度の 1 軸の分岐の係数（制御棒の値 out、ほかの軸は公称の値。燃焼度は点の炉心の平均の燃焼度で、
  燃焼度の点の間は線形）に、点のノードの温度の入口の温度からの上がりの、出力の割合の 2 乗の重みの平均（燃料は √T。テーブルの燃料の軸の
  内挿の変数）を掛け、ホウ素の係数（熱水力なしの点のホウ素を挟む格子点の割線）で割った値の負。2 乗の重みは、随伴の中性子束を中性子束に
  比例とみた 1 次の摂動の重み（近似。確かめていない）
- `--json`: 名前ごとに `points`（compare_core_guess.py の項目に、熱水力なしと熱水力の固有値計算の回数 `solves`、ノードの減速材と燃料の温度の
  最大、重みの平均の上がり `moderator_rise`、`fuel_rise`（K）、`effect`、`estimate` の燃料と減速材と和）、`midpoints`、`runs` の `isothermal` と
  `thermal`（compare_core_guess.py の `run_of`）、`thermal`（ログの熱水力の行の数、繰り返しの回数の [最小, 最大]、出口の温度と燃料の温度の
  最大）。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import compare_core_guess
import msgpack

THERMAL = re.compile(
    r"thermal-hydraulics iterations=(\d+) change=\S+ outlet=([0-9.eE+-]+) fuel_average=\S+ fuel_max=([0-9.eE+-]+)"
)
# 2 つの燃焼で同じであるべき履歴の項目（状態は温度の書き方が違うので別に比べる）
SHARED = ("specific_power", "axial_mesh", "assemblies", "heavy_metal")


def rho(k: float) -> float:
    return 1.0 - 1.0 / k


def rest(state: dict) -> dict:
    """状態の、温度のほかの項目"""
    return {k: v for k, v in state.items() if k not in ("moderator", "fuel", "inlet")}


def check_states(isothermal: dict, thermal: dict) -> None:
    if isothermal.keys() != thermal.keys():
        raise SystemExit("the histories name different states")
    for name, t in thermal.items():
        i = isothermal[name]
        inlet = t.get("inlet")
        if inlet is None or (i.get("moderator"), i.get("fuel")) != (inlet, inlet):
            raise SystemExit(
                f"state {name}: the isothermal temperatures {i.get('moderator')}, {i.get('fuel')} are not the inlet {inlet}"
            )
        if rest(i) != rest(t):
            raise SystemExit(f"state {name} differs beyond the temperatures")


def table_rho(table: dict, burnup: float, at: dict[str, int]) -> float:
    """テーブルの k∞ の反応度を、燃焼度 `burnup`（燃焼度の点の間は線形）、制御棒の値 out、軸 `at`（名前 → 格子点の番号。ほかの軸は公称）で"""
    axes = table["axes"]
    (points,) = [a["points"] for a in axes if a["name"] == "burnup"]
    if not points[0] <= burnup <= points[-1]:
        raise SystemExit(f"burnup {burnup} outside the table points {points}")
    j = min(max(i for i, u in enumerate(points) if u <= burnup), len(points) - 2)
    w = (burnup - points[j]) / (points[j + 1] - points[j])

    def at_point(i: int) -> float:
        index = []
        for a in axes:
            if a["name"] == "burnup":
                index.append(i)
            elif a["name"] == "rod":
                index.append(a["values"].index("out"))
            else:
                index.append(at.get(a["name"], a["points"].index(a["nominal"])))
        try:
            k = table["states"][table["points"].index(index)]["k_inf"]
        except ValueError:
            raise SystemExit(f"the table has no point {index}") from None
        return rho(k)

    return (1.0 - w) * at_point(j) + w * at_point(j + 1)


def coefficients(table: dict, burnup: float, boron: float) -> dict[str, float]:
    """燃焼度 `burnup` の k∞ の反応度の係数: 燃料（/√K）、減速材（/K）は公称から 1 つ上の格子点の割線、ホウ素（/ppm）は `boron` を挟む格子点の割線"""
    axes = {a["name"]: a for a in table["axes"]}

    def rho_at(name: str, i: int) -> float:
        return table_rho(table, burnup, {name: i})

    def step(name: str, transform) -> float:
        a = axes[name]
        i = a["points"].index(a["nominal"])
        if i + 1 >= len(a["points"]):
            raise SystemExit(f"the {name} axis has no point above the nominal")
        x0, x1 = a["points"][i], a["points"][i + 1]
        return (rho_at(name, i + 1) - rho_at(name, i)) / (transform(x1) - transform(x0))

    b = axes["boron"]["points"]
    segments = [i for i in range(len(b) - 1) if b[i] <= boron <= b[i + 1]]
    if not segments:
        raise SystemExit(f"boron {boron} outside the table points {b}")
    i = segments[0]
    return {
        "fuel": step("fuel", math.sqrt),
        "moderator": step("moderator", lambda x: x),
        "boron": (rho_at("boron", i + 1) - rho_at("boron", i)) / (b[i + 1] - b[i]),
    }


def rises(point: dict, inlet: float) -> tuple[float, float, float]:
    """点のノードの温度の、出力の割合の 2 乗の重みの平均の上がり: 減速材 (K)、燃料 (K)、燃料 (√K)"""
    nodes = point.get("temperatures")
    if nodes is None:
        raise SystemExit(
            f"the thermal point at {point['efpd']} EFPD has no node temperatures"
        )
    total = moderator = fuel = root = 0.0
    for temperatures, powers in zip(nodes, point["powers"]):
        for t, p in zip(temperatures, powers):
            if p == 0.0:
                continue
            w = p * p
            total += w
            moderator += w * (t[0] - inlet)
            fuel += w * (t[1] - inlet)
            root += w * (math.sqrt(t[1]) - math.sqrt(inlet))
    return moderator / total, fuel / total, root / total


def compare(
    iso_core: dict, iso_log: str, th_core: dict, th_log: str, table: dict
) -> dict:
    for key in SHARED:
        if iso_core.get(key) != th_core.get(key):
            raise SystemExit(f"the histories differ in {key}")
    check_states(iso_core["states"], th_core["states"])
    a, b = th_core["points"], iso_core["points"]
    if [(p["efpd"], p["days"], p["state"]) for p in a] != [
        (p["efpd"], p["days"], p["state"]) for p in b
    ]:
        raise SystemExit("the histories have different points")
    points = []
    for p, q in zip(a, b):
        inlet = th_core["states"][p["state"]]["inlet"]
        moderator_rise, fuel_rise, root_rise = rises(p, inlet)
        c = coefficients(table, q["burnup"], q["boron"])
        estimate = {
            "fuel": -c["fuel"] * root_rise / c["boron"],
            "moderator": -c["moderator"] * moderator_rise / c["boron"],
        }
        estimate["total"] = estimate["fuel"] + estimate["moderator"]
        nodes = [t for row in p["temperatures"] for t in row]
        points.append(
            {
                "efpd": q["efpd"],
                "days": q["days"],
                "state": q["state"],
                "k": p["k"] - q["k"],
                "boron": p["boron"] - q["boron"],
                "burnup": compare_core_guess.largest_relative(
                    compare_core_guess.flat(p["burnups"]),
                    compare_core_guess.flat(q["burnups"]),
                ),
                "power": compare_core_guess.largest_relative(
                    compare_core_guess.flat(p["powers"]),
                    compare_core_guess.flat(q["powers"]),
                ),
                "isothermal_boron": q["boron"],
                "solves": [q["solves"], p["solves"]],
                "moderator_max": max(t[0] for t in nodes),
                "fuel_max": max(t[1] for t in nodes),
                "moderator_rise": moderator_rise,
                "fuel_rise": fuel_rise,
                "effect": p["boron"]
                - q["boron"]
                - (rho(p["k"]) - rho(q["k"])) / c["boron"],
                "estimate": estimate,
            }
        )
    mids = compare_core_guess.MIDPOINT.findall(th_log)
    iso_mids = compare_core_guess.MIDPOINT.findall(iso_log)
    if len(mids) != len(iso_mids):
        raise SystemExit(f"{len(mids)} and {len(iso_mids)} midpoints")
    midpoints = [
        {"k": float(m[0]) - float(f[0]), "boron": float(m[1]) - float(f[1])}
        for m, f in zip(mids, iso_mids)
    ]
    lines = THERMAL.findall(th_log)
    if len(lines) != len(a) + len(mids):
        raise SystemExit(
            f"{len(lines)} thermal-hydraulic lines for {len(a)} points and {len(mids)} midpoints"
        )
    iterations = [int(n) for n, _, _ in lines]
    return {
        "points": points,
        "midpoints": midpoints,
        "runs": {
            "isothermal": compare_core_guess.run_of(iso_core, iso_log),
            "thermal": compare_core_guess.run_of(th_core, th_log),
        },
        "thermal": {
            "solves": len(lines),
            "iterations": [min(iterations), max(iterations)],
            "outlet": max(float(o) for _, o, _ in lines),
            "fuel": max(float(f) for _, _, f in lines),
        },
    }


def report(summary: dict) -> list[str]:
    lines = []
    for name, c in summary["cases"].items():
        for p in c["points"]:
            e = p["estimate"]
            lines.append(
                f"{name}: {p['efpd']} EFPD（{p['state']}）: ホウ素 {p['isothermal_boron']:.2f} → {p['isothermal_boron'] + p['boron']:.2f} ppm、"
                f"熱水力の効き {p['effect']:+.2f} ppm、見込み {e['total']:+.2f} ppm（燃料 {e['fuel']:+.2f}、減速材 {e['moderator']:+.2f}）、"
                f"重みの平均の上がり 減速材 {p['moderator_rise']:.2f} K、燃料 {p['fuel_rise']:.2f} K、最大 減速材 {p['moderator_max']:.2f} K、"
                f"燃料 {p['fuel_max']:.2f} K、固有値計算 {p['solves'][0]} → {p['solves'][1]} 回"
            )
        i, t = c["runs"]["isothermal"], c["runs"]["thermal"]
        th = c["thermal"]
        lines.append(
            f"{name}: 固有値計算 {i['solves']} → {t['solves']} 回（{t['solves'] / i['solves']:.3f} 倍）、外側の反復 "
            f"{i['iterations']:,} → {t['iterations']:,} 回、固有値計算の時間 {i['eigenvalue_seconds']:.1f} → "
            f"{t['eigenvalue_seconds']:.1f} s、全体 {i['seconds']:.1f} → {t['seconds']:.1f} s（{t['seconds'] / i['seconds']:.3f} 倍）。"
            f"熱水力の繰り返し {th['solves']} 回の解で {th['iterations'][0]}〜{th['iterations'][1]} 回、出口の温度の最大 {th['outlet']:.2f} K、"
            f"燃料の温度の最大 {th['fuel']:.2f} K"
        )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    paths = [Path(p) for p in files.split(",")]
    if not name or len(paths) != 4 or not all(str(p) for p in paths):
        raise SystemExit(
            f"give <name>=<isothermal history>,<isothermal log>,<thermal history>,<thermal log>, got {arg}"
        )
    return name, paths


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 5 and args[0] == "--json" and args[2] == "--table":
        table_path = Path(args[3])
        table = msgpack.unpackb(table_path.read_bytes())
        cases = {}
        files = [table_path]
        for name, paths in map(parse, args[4:]):
            if name in cases:
                raise SystemExit(f"a second case named {name}")
            iso_history, iso_log, history, log = paths
            cases[name] = compare(
                msgpack.unpackb(iso_history.read_bytes()),
                iso_log.read_text(),
                msgpack.unpackb(history.read_bytes()),
                log.read_text(),
                table,
            )
            files += paths
        summary = {
            "format": "xlc-compare-core-thermal",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            "cases": cases,
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
