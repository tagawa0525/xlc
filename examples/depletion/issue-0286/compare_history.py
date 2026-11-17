"""同じ入力を既定の自己遮蔽（ADR-0251 の自由ガスの核種の絞り込み、ADR-0252 の求積 6 点）と厳密な経路（`ufg-free-gas =
{ share = false }`、`ufg-cp-points = 16`）で解いた 2 つの燃焼の履歴（lattice の `[depletion]` の出力、ADR-0173）を比べる。

使い方: python3 compare_history.py --json <要約の出力先> <ケース> <既定の履歴> <厳密の履歴> [<ケース> <既定> <厳密>...]
        python3 compare_history.py --report <要約の JSON>

- `--json`: ケースごとに、核種の並び、燃焼する物質（種類、セル、体積）、燃焼度の点が同じことを確かめてから、燃焼点ごとの k∞ の
  反応度の差 1/k(既定) − 1/k(厳密)（pcm）、物質の種類ごとの核種の数密度の相対の差 N(既定) / N(厳密) − 1（種類の領域の体積の
  重みの平均と、領域ごとの差の絶対値の最大）、ピン出力の相対の差の絶対値の最大を書く。比べる核種は検証の記録の集合体の表と同じ
  （compare_depletion_assembly.py の `nuclides_of`。MOX のゾーン `mox-*` は燃料の核種）。`inputs` は読んだファイルのパス（実行した
  cwd から）と md5。履歴が少数群定数を持てば、比べずに止める
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ASSEMBLY = load(
    "compare_depletion_assembly", HERE.parent / "compare_depletion_assembly.py"
)


def nuclides(kind: str) -> list[str]:
    """物質の種類の比べる核種。"""
    return (
        ASSEMBLY.FUEL_NUCLIDES
        if kind.startswith("mox-")
        else ASSEMBLY.nuclides_of(kind)
    )


def relative(ours: float, theirs: float, what: str) -> float | None:
    """厳密 `theirs` に対する既定 `ours` の相対の差。両方 0 なら比べない（None）。片方だけ 0 は止める。"""
    if theirs == 0.0:
        if ours == 0.0:
            return None
        raise SystemExit(f"{what}: {ours} on the defaults against 0 on the exact path")
    if math.isnan(ours) or math.isnan(theirs):
        raise SystemExit(f"{what}: NaN")
    return ours / theirs - 1.0


def check_same_run(default: dict, exact: dict) -> None:
    for h in (default, exact):
        if h["format"] != "xlc-depletion-history":
            raise SystemExit(f"unexpected format {h['format']}")
    if default["nuclides"] != exact["nuclides"]:
        raise SystemExit("the histories list other nuclides")
    keys = ("material", "cell", "volume")
    ours = [tuple(r[k] for k in keys) for r in default["regions"]]
    theirs = [tuple(r[k] for k in keys) for r in exact["regions"]]
    if ours != theirs:
        raise SystemExit("the histories burn other regions")
    ours = [p["burnup"] for p in default["points"]]
    theirs = [p["burnup"] for p in exact["points"]]
    if ours != theirs:
        raise SystemExit(f"burnup points differ: {ours} vs {theirs}")
    for h in (default, exact):
        if any(p["few_group"] is not None for p in h["points"]):
            raise SystemExit("few-group constants in the history are not compared")


def largest(values: list[float | None], burnup: list[float]) -> dict | None:
    """絶対値が最大の値とその燃焼度。比べた値が無ければ None。"""
    found = [(v, b) for v, b in zip(values, burnup) if v is not None]
    if not found:
        return None
    value, at = max(found, key=lambda x: abs(x[0]))
    return {"value": value, "burnup": at}


def pin_difference(
    ours: list[list[float]], theirs: list[list[float]], burnup: float
) -> float:
    """ピン出力の相対の差の絶対値の最大（核分裂の無いセルは両方 0 であることを確かめて除く）。"""
    if len(ours) != len(theirs) or any(len(a) != len(b) for a, b in zip(ours, theirs)):
        raise SystemExit(f"pin maps of other shapes at {burnup} MWd/kg")
    worst = 0.0
    for i, (a, b) in enumerate(zip(ours, theirs)):
        for j, (x, y) in enumerate(zip(a, b)):
            d = relative(x, y, f"pin ({i}, {j}) at {burnup} MWd/kg")
            if d is not None:
                worst = max(worst, abs(d))
    return worst


def summarize(default: dict, exact: dict) -> dict:
    """2 つの履歴を比べた要約（モジュールの説明の `--json`）。`max_density` は種類と核種をまたいだ体積の平均の差の絶対値の最大。"""
    check_same_run(default, exact)
    burnup = [p["burnup"] for p in exact["points"]]
    index = {name: i for i, name in enumerate(exact["nuclides"])}
    regions = exact["regions"]
    k = []
    for ours, theirs in zip(default["points"], exact["points"]):
        k.append(
            {
                "burnup": theirs["burnup"],
                "default": ours["k_inf"],
                "exact": theirs["k_inf"],
                "drho_pcm": (1.0 / ours["k_inf"] - 1.0 / theirs["k_inf"]) * 1e5,
            }
        )
    kinds: dict[str, dict] = {}
    for kind in dict.fromkeys(r["material"] for r in regions):
        members = [r for r, region in enumerate(regions) if region["material"] == kind]
        total = sum(regions[r]["volume"] for r in members)
        kinds[kind] = {}
        for name in nuclides(kind):
            n = index[name]
            series, per_region = [], []
            for ours, theirs in zip(default["points"], exact["points"]):
                what = f"{kind} {name} at {theirs['burnup']} MWd/kg"
                mean = [
                    sum(p["densities"][r][n] * regions[r]["volume"] for r in members)
                    / total
                    for p in (ours, theirs)
                ]
                series.append(relative(*mean, what))
                diffs = [
                    relative(
                        ours["densities"][r][n],
                        theirs["densities"][r][n],
                        f"{what}, region {r}",
                    )
                    for r in members
                ]
                found = [d for d in diffs if d is not None]
                per_region.append(max(found, key=abs) if found else None)
            kinds[kind][name] = {
                "series": series,
                "max": largest(series, burnup),
                "region_max": largest(per_region, burnup),
            }
    pins = [
        pin_difference(ours["pin_powers"], theirs["pin_powers"], theirs["burnup"])
        for ours, theirs in zip(default["points"], exact["points"])
    ]
    found = [
        {"kind": kind, "nuclide": name, **v["max"]}
        for kind, by_name in kinds.items()
        for name, v in by_name.items()
        if v["max"] is not None
    ]
    return {
        "burnup": burnup,
        "k": k,
        "max_drho_pcm": largest([x["drho_pcm"] for x in k], burnup),
        "kinds": kinds,
        "max_density": max(found, key=lambda x: abs(x["value"])),
        "pins": pins,
        "max_pins": largest(pins, burnup),
    }


def read(path: Path) -> tuple[dict, dict]:
    data = path.read_bytes()
    return msgpack.unpackb(data, strict_map_key=False), {
        "path": str(path),
        "md5": hashlib.md5(data).hexdigest(),
    }


def write_json(out: Path, cases: list[tuple[str, Path, Path]]) -> None:
    summary = {}
    for name, default_path, exact_path in cases:
        default, default_input = read(default_path)
        exact, exact_input = read(exact_path)
        summary[name] = {
            "inputs": {"default": default_input, "exact": exact_input},
            **summarize(default, exact),
        }
    out.write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")


def at(found: dict | None, form: str, unit: str = "") -> str:
    return (
        "–"
        if found is None
        else f"{format(found['value'], form)}{unit}（{found['burnup']:g} MWd/kg）"
    )


def report(summary: dict) -> list[str]:
    lines = []
    for name, case in summary.items():
        first, last = case["k"][0], case["k"][-1]
        d = case["max_density"]
        lines.append(
            f"{name}: 0〜{case['burnup'][-1]:g} MWd/kg の {len(case['burnup'])} 点。"
            f"1/k(既定) − 1/k(厳密) の絶対値の最大 {at(case['max_drho_pcm'], '+.4f', ' pcm')}、"
            f"最後の点 {last['drho_pcm']:+.4f} pcm（k∞ 既定 {last['default']:.6f}、厳密 {last['exact']:.6f}）。"
            f"0 の点 {first['drho_pcm']:+.4f} pcm。"
            f"ピン出力の相対の差の最大 {at(case['max_pins'], '.2e')}。"
            f"数密度（体積の平均）の相対の差の最大 {d['kind']} の {d['nuclide']} {at(d, '+.2e')}"
        )
        for kind, by_name in case["kinds"].items():
            parts = [
                f"{n} {at(v['max'], '+.2e')}、領域の最大 {at(v['region_max'], '+.2e')}"
                for n, v in by_name.items()
            ]
            lines.append(
                f"  {kind} の数密度の相対の差（体積の平均）: " + "; ".join(parts)
            )
    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--json", type=Path)
    group.add_argument("--report", type=Path)
    parser.add_argument("cases", nargs="*")
    args = parser.parse_args()
    if args.report:
        if args.cases:
            parser.error("--report takes no cases")
        print("\n".join(report(json.loads(args.report.read_text()))))
        return
    if not args.cases or len(args.cases) % 3:
        parser.error(
            "give the cases as <name> <default history> <exact history> triples"
        )
    triples = [
        (args.cases[i], Path(args.cases[i + 1]), Path(args.cases[i + 2]))
        for i in range(0, len(args.cases), 3)
    ]
    write_json(args.json, triples)


if __name__ == "__main__":
    main()
