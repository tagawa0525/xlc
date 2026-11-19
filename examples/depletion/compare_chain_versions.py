"""2 つの版の ENDF から bin/chain で作った燃焼チェーン（ADR-0164、ADR-0262）を、値で比べる（#0303）。

使い方: python3 compare_chain_versions.py --json <要約の出力先> <チェーン a（msgpack）> <チェーン b（msgpack）>
        python3 compare_chain_versions.py --report <要約の JSON>

xlc のチェーンと OpenMC のチェーンを比べるのは compare_chain.py で、こちらは同じ形式の 2 つのチェーンを比べる。

- `--json`: 核種を名前で対応させ、量ごとに比べた項目の数、違う項目の数、相対差がしきい値を超える項目の数、相対差の大きい順の
  先頭 `LARGEST` 件を書き、違う項目のうち片方だけが値を持つもの（`one_sided`）を別に数える。量は次のとおりで、項目は両方のチェーンにある核種のもの:
  - `excitation`: 核種の励起エネルギー
  - `half_life`: 半減期。片方だけが安定なら、安定の側を null にして相対差を null（いちばん先に並べる）
  - `decay_energy`: 崩壊の平均のエネルギーの 3 成分（両方で不安定な核種）
  - `decay_modes`: 崩壊の（型、生まれる核種）ごとの分岐比（両方で不安定な核種）。同じ型と生成物を 2 度書いたものは足す。
    片方に無い分岐は値なし
  - `fission_q`: 核分裂の Q 値（どちらかで核分裂する核種。片方だけなら、もう片方は値なし）
  - `yield_sources`: 収率を取った核種（どちらかで核分裂する核種。片方だけなら、もう片方は値なし）
  - `yield_energies`: 収率の入射エネルギーの並び（どちらかで核分裂する核種。片方だけなら、もう片方は値なし）
  - `fission_yields`: （核分裂する核種、入射エネルギー、生成物）ごとの独立収率。両方にある表で片方に無い生成物の
    収率は 0、片方の版だけにある入射エネルギーの表の項目は、もう片方を値なしにする
  - `reactions`: （核種、MT）ごとの経路（生成物の状態と生まれる核種）の集合。片方だけの反応も違いに数える
  数値は値が等しいかで比べ（ENDF の書き方と、チェーンの中の並びには依らない）、相対差は |a − b| / max(|a|, |b|)（片方が 0 なら 1）。
  あわせて、（核分裂する核種、入射エネルギー）ごとの収率の差の絶対値の和と最大（`yield_tables`）、片方の版だけにある表の側と生成物の数（`one_sided_tables`）と、`WATCHED_PARENTS` と
  `WATCHED_NUCLIDES` の半減期、分岐比、核分裂の Q 値、`WATCHED_PARENTS` から `WATCHED_NUCLIDES` への収率の両方の値（等しくても。
  `watched`）を書く。`inputs` は読んだファイルのパス（ステージの cwd から）と md5
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
    # I-135 の β− の分岐の片方の娘（版で分岐比が違う。#0303）
    "Xe-135m",
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
        if state is None:
            where = ""
        elif state[1] is None:
            # MF8 が基底状態を挙げない反応の、挙げた状態を除いた残り
            where = f"MF{state[0]} の残り（基底状態）: "
        else:
            where = f"MF{state[0]} LFS {state[1]}: "
        parts.append(where + products_text(products))
    return "; ".join(parts)


def no_fission() -> dict:
    """片方の版だけで核分裂する核種の、もう片方の核分裂（どれも値なし。核種ごとに作る）"""
    return {"q": None, "yield_source": None, "energies": None, "tables": {}}


def energies_text(energies) -> str | None:
    # repr は float を往復できる桁で書くので、値が違えば文字も違う
    return None if energies is None else " ".join(repr(e) for e in energies)


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
                        # チェーンの分岐比は正なので、片方に無い分岐は 0 でなく値なし
                        "a": xd["modes"].get(key),
                        "b": yd["modes"].get(key),
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
        if xf is None and yf is None:
            continue
        # 片方の版だけで核分裂する核種は、もう片方の Q、収率を取った核種、入射エネルギー、表を値なしにする
        xf = xf or no_fission()
        yf = yf or no_fission()
        items["fission_q"].append({"nuclide": name, "a": xf["q"], "b": yf["q"]})
        items["yield_sources"].append(
            {"nuclide": name, "a": xf["yield_source"], "b": yf["yield_source"]}
        )
        items["yield_energies"].append(
            {
                "nuclide": name,
                "a": energies_text(xf["energies"]),
                "b": energies_text(yf["energies"]),
            }
        )
        for e in sorted(set(xf["tables"]) | set(yf["tables"])):
            xt, yt = xf["tables"].get(e), yf["tables"].get(e)
            for product in sorted(set(xt or {}) | set(yt or {})):
                items["fission_yields"].append(
                    {
                        "nuclide": name,
                        "energy": e,
                        "product": product,
                        # 片方の版だけにある入射エネルギーの表は、もう片方を値なしにする
                        "a": None if xt is None else xt.get(product, 0.0),
                        "b": None if yt is None else yt.get(product, 0.0),
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
            # 片方だけが値を持つ項目（片方だけの反応、分岐、収率の表、片方だけ安定な核種）
            "one_sided": sum(1 for i in diff if i["a"] is None or i["b"] is None),
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
    one_sided_tables = []
    for (name, e), ys in by_table.items():
        missing_a = [i["a"] is None for i in ys]
        missing_b = [i["b"] is None for i in ys]
        if any(missing_a) or any(missing_b):
            # 値なしは表の単位で入る（collect）。生成物ごとに混ざれば、片方を決められない
            if not (all(missing_a) or all(missing_b)):
                raise ValueError(f"{name} E={e}: a table partly in one chain")
            only = "b" if all(missing_a) else "a"
            one_sided_tables.append(
                {"nuclide": name, "energy": e, "only": only, "products": len(ys)}
            )
            continue
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
    for name in WATCHED_PARENTS + WATCHED_NUCLIDES:
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
                # どちらかの版で核分裂する核種だけが fission_q の項目を持つ（片方だけなら、もう片方は None）
                "fission_q": next(
                    (
                        {"a": i["a"], "b": i["b"]}
                        for i in items["fission_q"]
                        if i["nuclide"] == name
                    ),
                    None,
                ),
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
        "one_sided_tables": one_sided_tables,
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
        "| 量 | 比べた項目 | 違う項目 | 片方だけ | "
        + " | ".join(f"相対差 > {t:g}" for t in THRESHOLDS)
        + " |"
    )
    lines.append(
        "| --- | --- | --- | --- | " + " | ".join("---" for _ in THRESHOLDS) + " |"
    )
    for q, s in summary["quantities"].items():
        above = " | ".join(str(s["above"][f"{t:g}"]) for t in THRESHOLDS)
        lines.append(
            f"| {q} | {s['compared']} | {s['differ']} | {s['one_sided']} | {above} |"
        )
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
    only = summary["one_sided_tables"]
    if only:
        lines.append("")
        lines.append("片方の版だけにある収率の表:")
        lines.extend(
            f"- {t['nuclide']} E={t['energy']:g} eV: {t['only']} だけ（生成物 {t['products']}）"
            for t in only
        )
    lines.append("")
    lines.append("## 注目する核種（等しくても両方の値）")
    for w in summary["watched"]["nuclides"]:
        if w.get("missing"):
            lines.append(f"- {w['nuclide']}: 片方に無い")
            continue
        h = w["half_life"]
        if h["a"] is None and h["b"] is None:
            lines.append(f"- {w['nuclide']}: 両方で安定")
            continue
        lines.append(
            f"- {w['nuclide']} 半減期: a {value(h['a'])} s、b {value(h['b'])} s（相対差 {value(relative(h['a'], h['b']))}）"
        )
        lines.extend(f"  - 分岐 {describe('decay_modes', m)}" for m in w["modes"])
        q = w["fission_q"]
        if q is not None:
            lines.append(
                f"  - 核分裂の Q: a {value(q['a'])} eV、b {value(q['b'])} eV（相対差 {value(relative(q['a'], q['b']))}）"
            )
    for y in summary["watched"]["yields"]:
        lines.append(f"- 収率 {describe('fission_yields', y)}")
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 4 and args[0] == "--json":
        files = [Path(args[2]), Path(args[3])]
        data = [f.read_bytes() for f in files]
        summary = summarize(*(msgpack.unpackb(d, raw=False) for d in data))
        # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
        summary["inputs"] = [
            {"path": str(f), "md5": hashlib.md5(d).hexdigest()}
            for f, d in zip(files, data)
        ]
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
