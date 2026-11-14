"""#0285 の段 1（ADR-0250）: 本番の形の 1 燃焼点のテーブル（inputs.py の PRODUCTION）を作った時間の内訳と、テーブルの差。

時間は、ログの状態の区切り（`lattice state=` の行）の間の区間を状態ごとに足し、状態を自己遮蔽の種類で分ける:

- ufg: UFG の自己遮蔽（`shielding` の区間）だけを解いた状態（今の解き方の参照と、比の基準）
- ratio-first: 比の自己遮蔽（`ratio-shielding` の区間）を、その燃料の温度の点で最初に解いた状態（基準の静止した標的の減速計算を
  解く。基準の UFG を同じ状態で先に解いたときは、その区間も入る）
- ratio: ほかの、比の自己遮蔽を解いた状態
- remix: どちらも解かない状態（ホウ素の軸の、参照の自己遮蔽を作り直す状態）

種類ごとに、状態の数、状態の `seconds` の和、UFG と比の自己遮蔽の区間の数と和、状態の最初の MOC の区間（状態の固有値計算）と
ほかの MOC の区間（微分の固有値計算）の数と和、微分（`derivative` の行の `seconds`。作り直しと MOC を含む）の数と和。全体はログの
`branches` の行の `total_seconds`。状態ごとの値も `states` に置く。

テーブルの差は、PAIRS の組（解き方, 基準）ごとの、状態ごとの k∞ の差 (k − k_基準) × 1e5 と反応度の差 (1/k_基準 − 1/k) × 1e5
（ADR-0200 の判定の値は反応度）、少数群定数（`ITEMS`）と微視的断面積（`MICRO_ITEMS`。核種ごと）の状態ごとの相対差の最大。
微視的断面積の反応（`reactions`）と、数密度の微分は比べない。

使い方: python3 compare_production.py --json <要約の出力先>（このディレクトリで）
"""

import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# 値の読み方と相対差（非有限の値で止める）は compare.py と同じものを使う
COMPARE = load("compare", HERE / "compare.py")
fields, md5, relative = COMPARE.fields, COMPARE.md5, COMPARE.relative

VARIANTS = ["nominal", "nominal-cp6", "f1-static-cp6"]
# テーブルの差: (解き方, 比べる基準)。求積 6 点の効きと、同じ 6 点の上での比の方式の効き
PAIRS = [
    ("nominal-cp6", "nominal"),
    ("f1-static-cp6", "nominal-cp6"),
    ("f1-static-cp6", "nominal"),
]
ITEMS = ["absorption", "nu_fission", "kappa_fission", "diffusion", "scatter", "adf"]
MICRO_ITEMS = ["total", "nu_fission", "fission", "scatter", "thermal_fission"]
SPAN = re.compile(
    r" INFO (shielding|ratio-shielding|moc)(?:\{[^}]*\})?: xlc_lattice::(?:shielding|moc): close "
    r"time\.busy=([\d.]+)(ms|s|µs|ns)"
)


def log(text: str) -> dict:
    """状態ごとの区間の秒（`states`）と、`branches` の行（`branches`）"""
    states, branches = [], None

    def fresh() -> dict:
        return {"shielding": [], "ratio-shielding": [], "moc": [], "derivative": []}

    current = fresh()
    for line in text.splitlines():
        span = SPAN.search(line)
        if span:
            current[span.group(1)].append(
                float(span.group(2)) * COMPARE.UNITS[span.group(3)]
            )
        elif " derivative " in line and "nuclide=" in line:
            current["derivative"].append(float(fields(line)["seconds"]))
        elif "lattice state=" in line:
            f = fields(line)
            states.append(
                {**current, "state": int(f["state"]), "seconds": float(f["seconds"])}
            )
            current = fresh()
        elif " branches shieldings=" in line:
            branches = {k: float(v) for k, v in fields(line).items()}
    if branches is None:
        sys.exit("the log has no branches line (the run did not finish)")
    return {"states": states, "branches": branches}


def kinds(states: list[dict], fuels: int) -> list[str]:
    """状態の種類（モジュールの説明）。状態の番号は軸の直積の順で、燃料の温度が最も速く変わる（`fuels` はその点の数）"""
    seen, out = set(), []
    for s in states:
        if s["ratio-shielding"]:
            fuel = s["state"] % fuels
            out.append("ratio" if fuel in seen else "ratio-first")
            seen.add(fuel)
        elif s["shielding"]:
            out.append("ufg")
        else:
            out.append("remix")
    return out


def breakdown(parsed: dict, fuels: int) -> dict:
    out: dict = {}
    rows = []
    for s, name in zip(parsed["states"], kinds(parsed["states"], fuels)):
        k = out.setdefault(
            name,
            {
                "states": 0,
                "seconds": 0.0,
                "ufg": [0, 0.0],
                "ratio": [0, 0.0],
                "moc_state": [0, 0.0],
                "moc_other": [0, 0.0],
                "derivative": [0, 0.0],
            },
        )
        parts = {
            "ufg": s["shielding"],
            "ratio": s["ratio-shielding"],
            "moc_state": s["moc"][:1],
            "moc_other": s["moc"][1:],
            "derivative": s["derivative"],
        }
        k["states"] += 1
        k["seconds"] += s["seconds"]
        for key, values in parts.items():
            k[key][0] += len(values)
            k[key][1] += sum(values)
        rows.append(
            {"state": s["state"], "kind": name, "seconds": s["seconds"]}
            | {key: sum(values) for key, values in parts.items()}
        )
    return {
        "kinds": out,
        "state_seconds": sum(s["seconds"] for s in parsed["states"]),
        "total_seconds": parsed["branches"]["total_seconds"],
        "shieldings": parsed["branches"]["shieldings"],
        "ratios": parsed["branches"].get("ratios"),
        "states": rows,
    }


def table(path: Path) -> dict:
    """`points`（[状態] の軸の格子点の添字: 燃焼度、ホウ素、減速材、燃料）、`states`（[状態] の群定数）、`fuels`（燃料の温度の
    点の数）、`micro`（微視的断面積: `nuclides` と [状態][核種]）"""
    data = msgpack.unpackb(path.read_bytes(), strict_map_key=False)
    # 比の最初の状態（`kinds`）とログの状態の番号の読み替えは、燃焼度 1 点で燃料の温度が最も速く変わる軸の順を前提にする
    names = [a["name"] for a in data["axes"]]
    counts = [len(a["points"]) for a in data["axes"]]
    if names != ["burnup", "boron", "moderator", "fuel"] or counts[0] != 1:
        sys.exit(f"{path}: axes {names} with {counts} points, not one burnup point")
    return {
        "points": data["points"],
        "states": data["states"],
        "counts": counts,
        "fuels": counts[3],
        "micro": data["micro"],
    }


def point(state: int, counts: list[int]) -> list[int]:
    """ログの状態の番号（軸の直積の順。最後の軸が最も速く変わる）の格子点の添字"""
    out = []
    for n in reversed(counts):
        out.append(state % n)
        state //= n
    return out[::-1]


def worst(per: list[float], points: list) -> dict:
    at = max(range(len(per)), key=lambda i: abs(per[i]))
    return {"value": per[at], "state": at, "point": points[at]}


def differences(variant: dict, reference: dict) -> dict:
    points, states, base = variant["points"], variant["states"], reference["states"]
    if points != reference["points"]:
        sys.exit("the tables have different states")
    k = [(s["k_inf"], b["k_inf"]) for s, b in zip(states, base)]
    if not all(math.isfinite(x) and math.isfinite(y) for x, y in k):
        sys.exit("k-inf is not finite")
    pcm = [(x - y) * 1e5 for x, y in k]
    rho = [(1.0 / y - 1.0 / x) * 1e5 for x, y in k]
    items = {
        item: worst([relative(s[item], b[item]) for s, b in zip(states, base)], points)
        for item in ITEMS
    }
    micro, base_micro = variant["micro"], reference["micro"]
    if micro["nuclides"] != base_micro["nuclides"]:
        sys.exit("the tables have different micro nuclides")
    micro_items = {}
    for item in MICRO_ITEMS:
        per = []
        for s, b in zip(micro["states"], base_micro["states"]):
            # 核分裂しない核種の thermal_fission は値が無い（両方に無いことだけを許す）
            if any((x[item] is None) != (y[item] is None) for x, y in zip(s, b)):
                sys.exit(f"micro {item} is missing from one table")
            values = [
                0.0 if x[item] is None else relative(x[item], y[item])
                for x, y in zip(s, b)
            ]
            n = max(range(len(values)), key=lambda i: values[i])
            per.append((values[n], micro["nuclides"][n]))
        found = worst([v for v, _ in per], points)
        micro_items[item] = found | {"nuclide": per[found["state"]][1]}
    return {
        "k_pcm": pcm,
        "rho_pcm": rho,
        "worst_k": worst(pcm, points),
        "worst_rho": worst(rho, points),
        "items": items,
        "micro": micro_items,
    }


def main() -> None:
    args = sys.argv[1:]
    if len(args) != 2 or args[0] != "--json":
        raise SystemExit(__doc__)
    stem = "production-2A-60-"
    summary: dict = {"format": "xlc-compare-production-shielding", "variants": {}}
    inputs, tables = [], {}
    for v in VARIANTS:
        tpath, lpath = HERE / f"table-{stem}{v}.msgpack", HERE / f"{stem}{v}.log"
        inputs += [
            {"path": tpath.name, "md5": md5(tpath)},
            {"path": lpath.name, "md5": md5(lpath)},
        ]
        tables[v] = table(tpath)
        timing = breakdown(log(lpath.read_text()), tables[v]["fuels"])
        # 時間の行の `state` はログの番号（軸の直積の順）、差の `state` は表の状態の並び（`points`）の番号
        for row in timing["states"]:
            row["point"] = point(row["state"], tables[v]["counts"])
        summary["variants"][v] = {"timing": timing}
    summary["differences"] = {
        f"{v} - {b}": differences(tables[v], tables[b]) for v, b in PAIRS
    }
    summary["points"] = tables[VARIANTS[0]]["points"]
    summary["inputs"] = inputs
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
    for v, r in summary["variants"].items():
        t = r["timing"]
        print(
            f"{v}: total {t['total_seconds']:.0f} s, states {t['state_seconds']:.0f} s"
        )
        for name, k in t["kinds"].items():
            print(
                f"  {name}: {k['states']} states {k['seconds']:.0f} s; "
                + ", ".join(
                    f"{key} {k[key][0]} × {k[key][1]:.0f} s"
                    for key in ["ufg", "ratio", "moc_state", "moc_other", "derivative"]
                )
            )
    for pair, r in summary["differences"].items():
        print(
            f"{pair}: k {r['worst_k']['value']:+.3f} pcm {r['worst_k']['point']}, "
            f"rho {r['worst_rho']['value']:+.3f} pcm {r['worst_rho']['point']}; "
            + " ".join(
                f"{i} {x['value']:.1e}{x['point']}" for i, x in r["items"].items()
            )
            + "; micro "
            + " ".join(
                f"{i} {x['value']:.1e}{x['point']}{x['nuclide']}"
                for i, x in r["micro"].items()
            )
        )


if __name__ == "__main__":
    main()
