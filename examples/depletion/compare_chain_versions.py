"""2 つの版の ENDF から bin/chain で作った燃焼チェーン（ADR-0164、ADR-0262）を、値で比べる（#0303）。

使い方: python3 compare_chain_versions.py --json <要約の出力先> <チェーン a（msgpack）> <チェーン b（msgpack）>
        python3 compare_chain_versions.py --report <要約の JSON>

xlc のチェーンと OpenMC のチェーンを比べるのは compare_chain.py で、こちらは同じ形式の 2 つのチェーンを比べる。

- `--json`: 核種を名前で対応させ、量ごとに比べた項目の数、違う項目の数、相対差がしきい値を超える項目の数、相対差の大きい順の
  先頭 `LARGEST` 件を書く。量は次のとおりで、項目は両方のチェーンにある核種のもの:
  - `excitation`: 核種の励起エネルギー
  - `half_life`: 半減期。片方だけが安定なら、安定の側を null にして相対差を null（いちばん先に並べる）
  - `decay_energy`: 崩壊の平均のエネルギーの 3 成分（両方で不安定な核種）
  - `decay_modes`: 崩壊の（型、生まれる核種）ごとの分岐比（両方で不安定な核種）。同じ型と生成物を 2 度書いたものは足す
  - `fission_q`: 核分裂の Q 値（両方で核分裂する核種）
  - `yield_sources`: 収率を取った核種（両方で核分裂する核種）
  - `yield_energies`: 収率の入射エネルギーの並び（両方で核分裂する核種）
  - `fission_yields`: （核分裂する核種、両方にある入射エネルギー、生成物）ごとの独立収率。片方に無い生成物の収率は 0
  - `reactions`: （核種、MT）ごとの経路（生成物の状態と生まれる核種）の集合。片方だけの反応も違いに数える
  数値は値が等しいかで比べ（ENDF の書き方と、チェーンの中の並びには依らない）、相対差は |a − b| / max(|a|, |b|)（片方が 0 なら 1）。
  あわせて、（核分裂する核種、入射エネルギー）ごとの収率の差の絶対値の和と最大（`yield_tables`）、`WATCHED_PARENTS` の収率と
  `WATCHED_NUCLIDES` の半減期、分岐比、収率の両方の値（等しくても。`watched`）を書く。`inputs` は読んだファイルのパス（ステージの
  cwd から）と md5
- `--report`: 要約の JSON を文章で描く
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

# 量ごとに一覧に書く、相対差の大きい項目の数（JSON は git で追跡するので、全部は書かない）
LARGEST = 30
# 違う項目を、相対差がこれを超えるかで数える（書き方の桁の違いと値の違いを分けて見る目安）
THRESHOLDS = [1e-6, 1e-3, 1e-2]
# 燃焼の反応度に効きうる量を、等しくても両方の値で書く（#0303。どれが効くかはこのツールでは計算しない）
WATCHED_PARENTS = ["U-235", "U-238", "Pu-239", "Pu-241"]
WATCHED_NUCLIDES = [
    "Xe-135",
    "I-135",
    "Sm-149",
    "Pm-149",
    "Nd-143",
    "Rh-103",
    "Cs-133",
]
ENERGY_PARTS = ["light", "electromagnetic", "heavy"]
QUANTITIES = [
    "excitation",
    "half_life",
    "decay_energy",
    "decay_modes",
    "fission_q",
    "yield_sources",
    "yield_energies",
    "fission_yields",
    "reactions",
]


def relative(a, b) -> float | None:
    """|a − b| / max(|a|, |b|)。片方が None なら None"""
    if a is None or b is None:
        return None
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else 0.0


def named(chain: dict) -> dict:
    """名前から核種。生成物の添字を名前にした写し"""
    if chain.get("format") != "xlc-chain":
        raise SystemExit("not an xlc chain")
    names = [n["name"] for n in chain["nuclides"]]
    if len(set(names)) != len(names):
        raise SystemExit("a nuclide appears twice in the chain")

    def products(ps):
        return tuple(sorted((names[p["nuclide"]], p["count"]) for p in ps))

    result = {}
    for n in chain["nuclides"]:
        decay = None
        if n["decay"] is not None:
            modes = {}
            for m in n["decay"]["modes"]:
                key = (m["mode"], products(m["products"]))
                modes[key] = modes.get(key, 0.0) + m["ratio"]
            decay = {
                "half_life": n["decay"]["half_life"],
                "energy": n["decay"]["energy"],
                "modes": modes,
            }
        reactions = {}
        for r in n["reactions"]:
            channels = []
            for c in r["channels"]:
                fs = c["final_state"]
                state = None if fs is None else (fs["mf"], fs["lfs"])
                channels.append((state, products(c["products"])))
            reactions[r["mt"]] = tuple(sorted(channels, key=repr))
        fission = None
        if n["fission"] is not None:
            f = n["fission"]
            tables = {}
            for e, table in zip(f["energies"], f["yields"]):
                ys = {}
                for y in table:
                    product = names[y["nuclide"]]
                    ys[product] = ys.get(product, 0.0) + y["value"]
                tables[e] = ys
            fission = {
                "q": f["q"],
                "yield_source": f["yield_source"],
                "energies": list(f["energies"]),
                "tables": tables,
            }
        result[n["name"]] = {
            "excitation": n["excitation"],
            "decay": decay,
            "reactions": reactions,
            "fission": fission,
        }
    return result


def products_text(products) -> str:
    return " + ".join(f"{c}×{p}" if c != 1 else p for p, c in products) or "なし"


def channel_text(channels) -> str:
    parts = []
    for state, products in channels:
        where = "" if state is None else f"MF{state[0]} LFS {state[1]}: "
        parts.append(where + products_text(products))
    return "; ".join(parts)


def collect(a: dict, b: dict) -> dict:
    """量ごとの比べた項目。数値の項目は `a`、`b`、`relative` を持ち、ほかは `a`、`b` の文字"""
    items = {q: [] for q in QUANTITIES}
    common = sorted(set(a) & set(b))
    for name in common:
        x, y = a[name], b[name]
        items["excitation"].append(
            {"nuclide": name, "a": x["excitation"], "b": y["excitation"]}
        )
        xd, yd = x["decay"], y["decay"]
        if xd is not None or yd is not None:
            items["half_life"].append(
                {
                    "nuclide": name,
                    "a": None if xd is None else xd["half_life"],
                    "b": None if yd is None else yd["half_life"],
                }
            )
        if xd is not None and yd is not None:
            for part in ENERGY_PARTS:
                items["decay_energy"].append(
                    {
                        "nuclide": name,
                        "part": part,
                        "a": xd["energy"][part],
                        "b": yd["energy"][part],
                    }
                )
            for key in sorted(set(xd["modes"]) | set(yd["modes"])):
                items["decay_modes"].append(
                    {
                        "nuclide": name,
                        "mode": key[0],
                        "products": products_text(key[1]),
                        "a": xd["modes"].get(key, 0.0),
                        "b": yd["modes"].get(key, 0.0),
                    }
                )
        for mt in sorted(set(x["reactions"]) | set(y["reactions"])):
            xr, yr = x["reactions"].get(mt), y["reactions"].get(mt)
            items["reactions"].append(
                {
                    "nuclide": name,
                    "mt": mt,
                    "a": None if xr is None else channel_text(xr),
                    "b": None if yr is None else channel_text(yr),
                }
            )
        xf, yf = x["fission"], y["fission"]
        if xf is None or yf is None:
            continue
        items["fission_q"].append({"nuclide": name, "a": xf["q"], "b": yf["q"]})
        items["yield_sources"].append(
            {"nuclide": name, "a": xf["yield_source"], "b": yf["yield_source"]}
        )
        items["yield_energies"].append(
            {
                "nuclide": name,
                "a": " ".join(f"{e:g}" for e in xf["energies"]),
                "b": " ".join(f"{e:g}" for e in yf["energies"]),
            }
        )
        for e in xf["energies"]:
            if e not in yf["tables"]:
                continue
            xt, yt = xf["tables"][e], yf["tables"][e]
            for product in sorted(set(xt) | set(yt)):
                items["fission_yields"].append(
                    {
                        "nuclide": name,
                        "energy": e,
                        "product": product,
                        "a": xt.get(product, 0.0),
                        "b": yt.get(product, 0.0),
                    }
                )
    for q in QUANTITIES:
        for item in items[q]:
            if isinstance(item["a"], str) or isinstance(item["b"], str):
                item["relative"] = None if item["a"] != item["b"] else 0.0
            else:
                item["relative"] = relative(item["a"], item["b"])
    return items


def differs(item) -> bool:
    return item["a"] != item["b"]


def order(item):
    """相対差の大きい順。相対差の無い（片方だけ、文字の）違いを先に"""
    r = item["relative"]
    return (0, 0.0) if r is None else (1, -r)


def summarize(a: dict, b: dict) -> dict:
    """2 つのチェーン（msgpack を読んだもの）を比べた要約（モジュールの説明の `--json`）"""
    na, nb = named(a), named(b)
    items = collect(na, nb)
    quantities = {}
    for q in QUANTITIES:
        diff = sorted((i for i in items[q] if differs(i)), key=order)
        quantities[q] = {
            "compared": len(items[q]),
            "differ": len(diff),
            "above": {
                f"{t:g}": sum(
                    1 for i in diff if i["relative"] is None or i["relative"] > t
                )
                for t in THRESHOLDS
            },
            "largest": diff[:LARGEST],
        }
    by_table = {}
    for i in items["fission_yields"]:
        by_table.setdefault((i["nuclide"], i["energy"]), []).append(i)
    tables = []
    for (name, e), ys in by_table.items():
        worst = max(ys, key=lambda i: abs(i["a"] - i["b"]))
        tables.append(
            {
                "nuclide": name,
                "energy": e,
                "products": len(ys),
                "differ": sum(1 for i in ys if differs(i)),
                "sum_a": sum(i["a"] for i in ys),
                "sum_b": sum(i["b"] for i in ys),
                "sum_abs": sum(abs(i["a"] - i["b"]) for i in ys),
                "max_abs": abs(worst["a"] - worst["b"]),
                "max_product": worst["product"],
            }
        )
    watched_nuclides = []
    for name in WATCHED_NUCLIDES:
        if name not in na or name not in nb:
            watched_nuclides.append({"nuclide": name, "missing": True})
            continue
        watched_nuclides.append(
            {
                "nuclide": name,
                # 両方で安定な核種は half_life の項目を持たない
                "half_life": next(
                    (
                        {"a": i["a"], "b": i["b"]}
                        for i in items["half_life"]
                        if i["nuclide"] == name
                    ),
                    {"a": None, "b": None},
                ),
                "modes": [i for i in items["decay_modes"] if i["nuclide"] == name],
            }
        )
    watched_yields = [
        i
        for i in items["fission_yields"]
        if i["nuclide"] in WATCHED_PARENTS and i["product"] in WATCHED_NUCLIDES
    ]
    return {
        "format": "xlc-compare-chain-versions",
        "nuclides": {
            "a": len(na),
            "b": len(nb),
            "common": len(set(na) & set(nb)),
            "only_a": sorted(set(na) - set(nb)),
            "only_b": sorted(set(nb) - set(na)),
        },
        "quantities": quantities,
        "yield_tables": tables,
        "watched": {"nuclides": watched_nuclides, "yields": watched_yields},
    }


def value(v) -> str:
    if v is None:
        return "なし"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def describe(q: str, item: dict) -> str:
    where = item["nuclide"]
    if q == "decay_energy":
        where += f" {item['part']}"
    elif q == "decay_modes":
        where += f" {item['mode']} → {item['products']}"
    elif q == "fission_yields":
        where += f" E={item['energy']:g} eV → {item['product']}"
    elif q == "reactions":
        where += f" MT{item['mt']}"
    r = item["relative"]
    rel = "" if r is None else f"（相対差 {r:.3g}）"
    return f"{where}: a {value(item['a'])}、b {value(item['b'])}{rel}"


def report(summary: dict) -> list[str]:
    """要約の JSON を文章にする"""
    lines = []
    if summary["inputs"]:
        lines.append(
            "inputs: "
            + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
        )
    n = summary["nuclides"]
    lines.append(
        f"核種: a {n['a']}、b {n['b']}、共通 {n['common']}、a だけ {len(n['only_a'])}、b だけ {len(n['only_b'])}"
    )
    for side in ("only_a", "only_b"):
        if n[side]:
            lines.append(f"  {side}: {' '.join(n[side])}")
    lines.append("")
    lines.append(
        "| 量 | 比べた項目 | 違う項目 | "
        + " | ".join(f"相対差 > {t:g}" for t in THRESHOLDS)
        + " |"
    )
    lines.append("| --- | --- | --- | " + " | ".join("---" for _ in THRESHOLDS) + " |")
    for q, s in summary["quantities"].items():
        above = " | ".join(str(s["above"][f"{t:g}"]) for t in THRESHOLDS)
        lines.append(f"| {q} | {s['compared']} | {s['differ']} | {above} |")
    for q, s in summary["quantities"].items():
        if not s["largest"]:
            continue
        lines.append("")
        lines.append(
            f"## {q}（違う {s['differ']} 項目のうち、相対差の大きい {len(s['largest'])}）"
        )
        lines.extend(f"- {describe(q, i)}" for i in s["largest"])
    lines.append("")
    lines.append("## 収率の表（核分裂する核種、入射エネルギーごと）")
    lines.append(
        "| 核種 | E (eV) | 生成物 | 違う | 和 a | 和 b | 差の絶対値の和 | 差の最大 | その生成物 |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for t in summary["yield_tables"]:
        lines.append(
            f"| {t['nuclide']} | {t['energy']:g} | {t['products']} | {t['differ']} | {t['sum_a']:.6f} | "
            f"{t['sum_b']:.6f} | {t['sum_abs']:.3e} | {t['max_abs']:.3e} | {t['max_product']} |"
        )
    lines.append("")
    lines.append("## 注目する核種（等しくても両方の値）")
    for w in summary["watched"]["nuclides"]:
        if w.get("missing"):
            lines.append(f"- {w['nuclide']}: 片方に無い")
            continue
        h = w["half_life"]
        lines.append(
            f"- {w['nuclide']} 半減期: a {value(h['a'])} s、b {value(h['b'])} s（相対差 {value(relative(h['a'], h['b']))}）"
        )
        lines.extend(f"  - 分岐 {describe('decay_modes', m)}" for m in w["modes"])
    for y in summary["watched"]["yields"]:
        lines.append(f"- 収率 {describe('fission_yields', y)}")
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 4 and args[0] == "--json":
        files = [Path(args[2]), Path(args[3])]
        a, b = (msgpack.unpackb(f.read_bytes(), raw=False) for f in files)
        summary = {
            "format": "xlc-compare-chain-versions",
            # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(a, b),
        }
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
