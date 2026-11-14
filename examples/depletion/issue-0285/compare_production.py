"""#0285 の段 1（ADR-0250）: 本番の形の 1 燃焼点のテーブル（inputs.py の PRODUCTION）を作った時間の内訳と、テーブルの差。

時間は、ログの状態の区切り（`lattice state=` の行）の間の区間を状態ごとに足し、状態を自己遮蔽の種類で分ける:

- ufg: UFG の自己遮蔽（`shielding` の区間）を解いた状態（`every-state` と比の基準）
- ratio: 比の自己遮蔽（`ratio-shielding` の区間）を解いた状態
- remix: どちらも解かない状態（ホウ素の軸の、参照の自己遮蔽を作り直す状態）

種類ごとに、状態の数、状態の `seconds` の和、自己遮蔽の区間の和、MOC の区間（`moc` の close。微分の MOC を含む）の数と和、
微分（`derivative` の行の `seconds`。自己遮蔽の作り直しと MOC を含む）の数と和。全体はログの `branches` の行の `total_seconds`。

テーブルの差は、PAIRS の組（解き方, 基準）ごとの、状態ごとの k∞ の差 (k − k_基準) × 1e5（pcm）と、少数群定数（`ITEMS`）の
状態ごとの相対差の最大。

使い方: python3 compare_production.py --json <要約の出力先>（このディレクトリで）
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
VARIANTS = ["nominal", "nominal-cp6", "f1-static-cp6"]
# テーブルの差: (解き方, 比べる基準)。求積 6 点の効きと、同じ 6 点の上での比の方式の効き
PAIRS = [
    ("nominal-cp6", "nominal"),
    ("f1-static-cp6", "nominal-cp6"),
    ("f1-static-cp6", "nominal"),
]
ITEMS = ["absorption", "nu_fission", "kappa_fission", "diffusion", "scatter", "adf"]
UNITS = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "ns": 1e-9}
FIELD = re.compile(r'(\w+)=("[^"]*"|\S+)')
SPAN = re.compile(
    r" INFO (shielding|ratio-shielding|moc)(?:\{[^}]*\})?: xlc_lattice::(?:shielding|moc): close "
    r"time\.busy=([\d.]+)(ms|s|µs|ns)"
)


def fields(line: str) -> dict:
    return {k: v.strip('"') for k, v in FIELD.findall(line)}


def log(text: str) -> dict:
    """状態ごとの区間の秒（`states`）と、`branches` の行（`branches`）"""
    states, branches = [], None

    def fresh() -> dict:
        return {"shielding": [], "ratio-shielding": [], "moc": [], "derivative": []}

    current = fresh()
    for line in text.splitlines():
        span = SPAN.search(line)
        if span:
            current[span.group(1)].append(float(span.group(2)) * UNITS[span.group(3)])
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


def kind(state: dict) -> str:
    """状態の種類（モジュールの説明。ufg+ratio は比の基準を先に解いた比の状態）"""
    if state["shielding"] and state["ratio-shielding"]:
        return "ufg+ratio"
    if state["shielding"]:
        return "ufg"
    if state["ratio-shielding"]:
        return "ratio"
    return "remix"


def breakdown(parsed: dict) -> dict:
    out: dict = {}
    for s in parsed["states"]:
        k = out.setdefault(
            kind(s),
            {
                "states": 0,
                "seconds": 0.0,
                "ufg": [0, 0.0],
                "ratio": [0, 0.0],
                "moc": [0, 0.0],
                "derivative": [0, 0.0],
            },
        )
        k["states"] += 1
        k["seconds"] += s["seconds"]
        for name, key in [
            ("shielding", "ufg"),
            ("ratio-shielding", "ratio"),
            ("moc", "moc"),
            ("derivative", "derivative"),
        ]:
            k[key][0] += len(s[name])
            k[key][1] += sum(s[name])
    return {
        "kinds": out,
        "state_seconds": sum(s["seconds"] for s in parsed["states"]),
        "total_seconds": parsed["branches"]["total_seconds"],
        "shieldings": parsed["branches"]["shieldings"],
        "ratios": parsed["branches"].get("ratios"),
    }


def relative(a, b) -> float:
    x, y = flat(a), flat(b)
    if len(x) != len(y):
        sys.exit(f"{len(x)} vs {len(y)} values")
    worst = 0.0
    for u, v in zip(x, y):
        scale = max(abs(u), abs(v))
        if scale > 0.0:
            worst = max(worst, abs(u - v) / scale)
    return worst


def flat(value) -> list[float]:
    if value is None:
        return []
    if isinstance(value, list):
        return [x for v in value for x in flat(v)]
    return [float(value)]


def table(path: Path) -> tuple[list, list[dict]]:
    """（[状態] の軸の格子点の添字（燃焼度、ホウ素、減速材、燃料）、[状態] の群定数）"""
    data = msgpack.unpackb(path.read_bytes(), strict_map_key=False)
    return data["points"], data["states"]


def differences(variant: tuple, reference: tuple) -> dict:
    (points, states), (base_points, base) = variant, reference
    if points != base_points:
        sys.exit("the tables have different states")
    k = [(s["k_inf"] - b["k_inf"]) * 1e5 for s, b in zip(states, base)]
    worst = max(range(len(k)), key=lambda i: abs(k[i]))
    items = {}
    for item in ITEMS:
        per = [relative(s[item], b[item]) for s, b in zip(states, base)]
        at = max(range(len(per)), key=lambda i: per[i])
        items[item] = {"relative": per[at], "state": at}
    return {
        "k_pcm": k,
        "worst_k": {"pcm": k[worst], "state": worst, "point": points[worst]},
        "items": {item: x | {"point": points[x["state"]]} for item, x in items.items()},
    }


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def main() -> None:
    args = sys.argv[1:]
    if len(args) != 2 or args[0] != "--json":
        raise SystemExit(__doc__)
    stem = "production-2A-60-"
    summary: dict = {"format": "xlc-compare-production-shielding", "variants": {}}
    inputs = []
    for v in VARIANTS:
        tpath, lpath = HERE / f"table-{stem}{v}.msgpack", HERE / f"{stem}{v}.log"
        inputs += [
            {"path": tpath.name, "md5": md5(tpath)},
            {"path": lpath.name, "md5": md5(lpath)},
        ]
        summary["variants"][v] = {"timing": breakdown(log(lpath.read_text()))}
    summary["differences"] = {
        f"{v} - {b}": differences(
            table(HERE / f"table-{stem}{v}.msgpack"),
            table(HERE / f"table-{stem}{b}.msgpack"),
        )
        for v, b in PAIRS
    }
    summary["points"] = table(HERE / f"table-{stem}{VARIANTS[0]}.msgpack")[0]
    summary["inputs"] = inputs
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
    for v, r in summary["variants"].items():
        t = r["timing"]
        print(
            f"{v}: total {t['total_seconds']:.0f} s, states {t['state_seconds']:.0f} s"
        )
        for name, k in t["kinds"].items():
            print(
                f"  {name}: {k['states']} states {k['seconds']:.0f} s; UFG {k['ufg'][0]} × {k['ufg'][1]:.0f} s, "
                f"ratio {k['ratio'][0]} × {k['ratio'][1]:.0f} s, MOC {k['moc'][0]} × {k['moc'][1]:.0f} s, "
                f"derivative {k['derivative'][0]} × {k['derivative'][1]:.0f} s"
            )
    for pair, r in summary["differences"].items():
        print(
            f"{pair}: k {r['worst_k']['pcm']:+.2f} pcm (point {r['worst_k']['point']}), "
            + " ".join(
                f"{i} {x['relative']:.1e}({x['point']})" for i, x in r["items"].items()
            )
        )


if __name__ == "__main__":
    main()
