"""#0285 の段 1（ADR-0250）: 比で求めた分岐のテーブル（inputs.py の f1、f2 と候補）を、状態ごとに UFG で解いたテーブル（every）と
比べる。

燃料の種類と燃焼度（inputs.py の <種類>-<燃焼度>。今の既定で比べ直す入力は six-<種類>-<燃焼度>）ごとに、every と同じ 9 状態（減速材 3 × 燃料の温度 3）で:

- k∞ の差 (k − k_every) × 1e5（pcm）
- 少数群定数（`ITEMS`）の相対差の最大
- ピン出力: 2 群の形状関数を、無限媒質の 2 群の中性子束の比 φ2/φ1 = Σ(1→2) / (Σa2 + Σ(2→1)) と κΣf で重みを付けて足した
  セルの出力（ノードの平均で 1）の相対差の最大
- 係数の比: 燃料の温度の反応度 1/k(公称の燃料) − 1/k(燃料) を減速材の点ごとに、減速材の反応度 1/k(公称の減速材) − 1/k(減速材)
  を燃料の温度の点ごとに求め、every の値に対する比 − 1 の最大（絶対値）
- リングの反応率: ログ（target ring_rates の DEBUG）の、燃焼の領域の元の物質（リングの位置）ごとの捕獲（U-238、Pu-239、Gd-155、
  Gd-157。Gd は核分裂しないので捕獲が吸収）の相対差の最大
- 時間: ログの自己遮蔽の区間（UFG の `shielding` と比の `ratio-shielding`）の回数と平均の秒、状態の `seconds` の和

使い方: python3 compare.py --json <要約の出力先>（このディレクトリで）
"""

import hashlib
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
ITEMS = ["absorption", "nu_fission", "diffusion", "scatter", "adf"]
RING_NUCLIDES = ["U-238", "Pu-239", "Gd-155", "Gd-157"]
UNITS = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "ns": 1e-9}
FIELD = re.compile(r'(\w+)=("[^"]*"|\S+)')
SPAN = re.compile(
    r" INFO (shielding|ratio-shielding): xlc_lattice::shielding: close time\.busy=([\d.]+)(ms|s|µs|ns)"
)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


INPUTS = load("inputs", HERE / "inputs.py")


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def flat(value) -> list[float]:
    if isinstance(value, list):
        return [x for v in value for x in flat(v)]
    return [float(value)]


def relative(a, b) -> float:
    """b に対する a の、0 でない要素の相対差の最大"""
    x, y = flat(a), flat(b)
    if len(x) != len(y):
        sys.exit(f"{len(x)} vs {len(y)} values")
    worst = 0.0
    for u, v in zip(x, y):
        if not (math.isfinite(u) and math.isfinite(v)):
            sys.exit(f"values {u} and {v} are not finite")
        scale = max(abs(u), abs(v))
        if scale > 0.0:
            worst = max(worst, abs(u - v) / scale)
    return worst


def power(state: dict) -> list[float]:
    """セルの 2 群の出力（ノードの平均で 1。モジュールの説明）"""
    a, s, kappa = state["absorption"], state["scatter"], state["kappa_fission"]
    flux = [1.0, s[0][1] / (a[1] + s[1][0])]
    rates = [kappa[g] * flux[g] for g in range(2)]
    form = state["pins"]["form"]
    cells = [
        sum(form[g][r][c] * rates[g] for g in range(2)) / sum(rates)
        for r in range(len(form[0]))
        for c in range(len(form[0][0]))
    ]
    return cells


def table(path: Path) -> list[dict]:
    data = msgpack.unpackb(path.read_bytes(), strict_map_key=False)
    counts = [len(a["points"]) for a in data["axes"]]
    if counts != [1, 1, 3, 3]:
        sys.exit(
            f"{path}: axes {counts}, not burnup 1 × boron 1 × moderator 3 × fuel 3"
        )
    return data["states"]


def fields(line: str) -> dict:
    return {k: v.strip('"') for k, v in FIELD.findall(line)}


def log(path: Path) -> dict:
    """ログの自己遮蔽の区間の秒、状態の秒、状態ごとの（物質, 核種）の捕獲"""
    spans = {"shielding": [], "ratio-shielding": []}
    seconds, rings, current = [], [], {}
    for line in path.read_text().splitlines():
        busy = SPAN.search(line)
        if busy:
            spans[busy.group(1)].append(float(busy.group(2)) * UNITS[busy.group(3)])
        elif "ring reaction rates per neutron produced" in line:
            f = fields(line)
            current[(f["material"], f["nuclide"])] = float(f["capture"])
        elif "lattice state=" in line:
            seconds.append(float(fields(line)["seconds"]))
            rings.append(current)
            current = {}
    return {"spans": spans, "seconds": seconds, "rings": rings}


def coefficients(states: list[dict]) -> dict:
    """[減速材][燃料] の状態の、燃料の温度と減速材の反応度（pcm。モジュールの説明）"""
    k = [[states[m * 3 + f]["k_inf"] for f in range(3)] for m in range(3)]
    fuel = {
        f"m{m}f{f}": (1.0 / k[m][1] - 1.0 / k[m][f]) * 1e5
        for m in range(3)
        for f in (0, 2)
    }
    moderator = {
        f"m{m}f{f}": (1.0 / k[1][f] - 1.0 / k[m][f]) * 1e5
        for m in (0, 2)
        for f in range(3)
    }
    return {"fuel": fuel, "moderator": moderator}


def worst_ratio(variant: dict, every: dict) -> dict:
    """係数ごとの、every に対する比 − 1 の絶対値の最大と、その点"""
    out = {}
    for kind in ("fuel", "moderator"):
        at, worst = None, 0.0
        for key, value in every[kind].items():
            r = variant[kind][key] / value - 1.0
            if abs(r) >= abs(worst):
                at, worst = key, r
        out[kind] = {"ratio_minus_one": worst, "at": at}
    return out


def ring_worst(variant: list[dict], every: list[dict]) -> dict:
    """核種ごとの、リングの捕獲の相対差の最大と、その状態と物質"""
    out = {}
    for nuclide in RING_NUCLIDES:
        worst = None
        for s, (v, e) in enumerate(zip(variant, every)):
            for (material, n), value in e.items():
                if n != nuclide or value == 0.0:
                    continue
                if (material, n) not in v:
                    sys.exit(f"state {s}: {material} {n} is missing from the ratio log")
                d = v[(material, n)] / value - 1.0
                if worst is None or abs(d) > abs(worst["relative"]):
                    worst = {"relative": d, "state": s, "material": material}
        if worst is not None:
            out[nuclide] = worst
    return out


def timing(parsed: dict) -> dict:
    spans = parsed["spans"]

    def mean(v: list[float]) -> float | None:
        return sum(v) / len(v) if v else None

    return {
        "ufg": {
            "count": len(spans["shielding"]),
            "mean_seconds": mean(spans["shielding"]),
        },
        "ratio": {
            "count": len(spans["ratio-shielding"]),
            "mean_seconds": mean(spans["ratio-shielding"]),
        },
        "state_seconds": sum(parsed["seconds"]),
    }


def compare(case: str, variants: list[str], inputs: list) -> dict:
    paths = {v: HERE / f"table-{case}-{v}.msgpack" for v in ["every", *variants]}
    logs = {v: HERE / f"{case}-{v}.log" for v in ["every", *variants]}
    for v in ["every", *variants]:
        inputs.append({"path": paths[v].name, "md5": md5(paths[v])})
    every = table(paths["every"])
    every_log = log(logs["every"])
    if len(every_log["rings"]) != 9:
        sys.exit(f"{logs['every']}: {len(every_log['rings'])} states")
    every_coefficients = coefficients(every)
    result = {
        "every": {"k_inf": [s["k_inf"] for s in every], "timing": timing(every_log)}
    }
    for v in variants:
        states = table(paths[v])
        parsed = log(logs[v])
        if len(parsed["rings"]) != 9:
            sys.exit(f"{logs[v]}: {len(parsed['rings'])} states")
        k = [(s["k_inf"] - e["k_inf"]) * 1e5 for s, e in zip(states, every)]
        worst_k = max(range(9), key=lambda i: abs(k[i]))
        result[v] = {
            "k_pcm": k,
            "worst_k": {"pcm": k[worst_k], "state": worst_k},
            "items": {
                item: max(relative(s[item], e[item]) for s, e in zip(states, every))
                for item in ITEMS
            },
            "power": max(relative(power(s), power(e)) for s, e in zip(states, every)),
            "coefficients": worst_ratio(coefficients(states), every_coefficients),
            "rings": ring_worst(parsed["rings"], every_log["rings"]),
            "timing": timing(parsed),
        }
    return result


def main() -> None:
    args = sys.argv[1:]
    if len(args) != 2 or args[0] != "--json":
        raise SystemExit(__doc__)
    cases: dict[str, list[str]] = {}
    for stem in INPUTS.names():
        name, burnup, shielding = stem.split("-", 2)
        if shielding != "every":
            cases.setdefault(f"{name}-{burnup}", []).append(shielding)
    # 今の既定で比べ直す入力（inputs.py の SIX）
    for name, burnup in INPUTS.SIX:
        cases[f"six-{name}-{burnup:g}"] = [
            v for v in INPUTS.SIX_SHIELDINGS if v != "every"
        ]
    present = {
        case: [v for v in variants if (HERE / f"table-{case}-{v}.msgpack").exists()]
        for case, variants in cases.items()
        if (HERE / f"table-{case}-every.msgpack").exists()
    }
    missing = sorted(
        f"{case}-{v}"
        for case, variants in cases.items()
        for v in ["every", *variants]
        if not (HERE / f"table-{case}-{v}.msgpack").exists()
    )
    inputs: list = []
    summary = {
        "format": "xlc-compare-ratio-shielding",
        "cases": {case: compare(case, v, inputs) for case, v in present.items()},
        "missing": missing,
        "inputs": inputs,
    }
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
    for case, result in summary["cases"].items():
        for v, r in result.items():
            if v == "every":
                t = r["timing"]
                print(
                    f"{case} every: UFG {t['ufg']['count']} × {t['ufg']['mean_seconds']:.1f} s"
                )
                continue
            c, t = r["coefficients"], r["timing"]
            rings = " ".join(
                f"{n} {x['relative']:+.2e}({x['state']},{x['material']})"
                for n, x in r["rings"].items()
            )
            ratio = t["ratio"]["mean_seconds"]
            print(
                f"{case} {v}: k {r['worst_k']['pcm']:+.2f} pcm (state {r['worst_k']['state']}), "
                + " ".join(f"{i} {x:.1e}" for i, x in r["items"].items())
                + f" power {r['power']:.1e}, fuel coef {c['fuel']['ratio_minus_one']:+.2%} ({c['fuel']['at']}), "
                f"moderator coef {c['moderator']['ratio_minus_one']:+.2%} ({c['moderator']['at']}), rings {rings}; "
                f"UFG {t['ufg']['count']} × {t['ufg']['mean_seconds'] or 0:.1f} s, ratio {t['ratio']['count']} × "
                f"{ratio or 0:.2f} s, states {t['state_seconds']:.0f} s"
            )
    if missing:
        print("missing: " + " ".join(missing))


if __name__ == "__main__":
    main()
