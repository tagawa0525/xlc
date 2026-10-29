"""前のサイクルの燃焼の履歴から restart した core の燃焼（ADR-0220）を、切らずに燃焼させた燃焼と比べ、装荷替えの 1 サイクルの点を並べる。

使い方: python3 compare_core_restart.py --json <要約の出力先> <切らない燃焼の履歴> [<名前>=<前の半分の履歴>,<後の半分の履歴> ...]
            [--reload <装荷替えの履歴>]
        python3 compare_core_restart.py --report <要約の JSON>

- `--json`: core の履歴（msgpack）を読んで次を書く。`inputs` は読んだファイルのパスと md5
  - 切る点ごと（`cuts`）: 切らない燃焼を前の半分の最後の点で切り、後の半分（前の半分の履歴から同じ装荷で restart した燃焼）の点を、
    EFPD と日数に切った点の値を足して切らない燃焼の点と組にする。前の半分、後の半分、切らない燃焼の軸方向のメッシュ、集合体、重金属の
    質量、追う核種と補助の核種、Xe の扱いが同じで、後の半分の燃焼点 0 のノードの燃焼度と数密度が前の半分の最後の点と同じで、組にした点の
    EFPD と日数が相対 1e-12 で同じであること（違えば止める）。点ごとに、k と切らない燃焼の k、反応度の差
    Δρ = (1/k_切らない − 1/k) × 1e5 pcm、ホウ素の差（後の半分 − 切らない）、ノードの燃焼度の差の絶対値の最大、追う核種の数密度の差の
    絶対値の最大（核種ごとの、切らない燃焼のその点の炉心の最大の数密度に対する比。下限は炉心のすべての核種の最大の数密度の 1e-12 倍）。
    最後の点の補助の核種の数密度の差も同じ比で（下限は追う核種と補助の核種の最大の 1e-12 倍。差の最大の核種とその核種の炉心の最大の数密度も）。前の半分の点と切らない燃焼の初めの点の
    k の差の絶対値の最大も書く（同じ計算なので 0 の見込み）
  - 装荷替え（`reload`）: 切らない燃焼を前のサイクルとして、装荷替えの燃焼の燃焼点 0 の、前のサイクルの履歴にある集合体（引き継いだ集合体）の
    ノードの燃焼度と数密度が前のサイクルの最後の点のその集合体と同じで、前のサイクルにない集合体（新しい集合体）のノードの燃焼度が 0 である
    こと（違えば止める）。点ごとの EFPD、日数、状態、k、ホウ素、固有値計算の回数、炉心の平均の燃焼度、集合体ごとの燃焼度（`assemblies` の順）
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

# 組にする点の EFPD と日数の相対の許容（後の半分は切った点から数え直すので、丸めの分だけ違う）
TIME_TOLERANCE = 1e-12


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def check_format(name: str, core: dict) -> None:
    if core.get("format") != "xlc-core-history":
        raise SystemExit(f"{name}: unexpected core format {core.get('format')}")


def check_same_core(name: str, core: dict, full: dict) -> None:
    """`core` が `full` と同じ炉心（メッシュ、集合体、重金属、核種、Xe の扱い）の燃焼であること"""
    check_format(name, core)
    for key in [
        "axial_mesh",
        "assemblies",
        "heavy_metal",
        "nuclides",
        "auxiliary",
        "xenon",
    ]:
        if core.get(key) != full.get(key):
            raise SystemExit(f"{name}: the {key} differs from the uncut depletion")


def close(a: float, b: float) -> bool:
    return abs(a - b) <= TIME_TOLERANCE * max(abs(a), abs(b), 1.0)


# CRAM の近似の誤差の水準の、ノードの最大の数密度に対する大きさ（crates/core/src/depletion.rs の `NEGATIVE_TOLERANCE`）。これより小さい差は
# 0 に近い値の揺らぎ（30 日の停止で 1e-27 ほどに消えた補助の核種など）で、核種ごとの最大に対する比にすると 1 の大きさになる
CRAM_ERROR = 1e-12


def largest_density(*densities) -> float:
    """数密度 [集合体][軸方向のノード][核種] の絶対値の最大"""
    return max(
        (abs(v) for d in densities for a in d for n in a for v in n), default=0.0
    )


def relative_density_differences(what: str, got, want, scale: float) -> float:
    return largest_density_difference(what, got, want, scale)[0]


def largest_density_difference(
    what: str, got, want, scale: float
) -> tuple[float, int | None, float]:
    """ノードの数密度 [集合体][軸方向のノード][核種] の差の絶対値の最大。核種ごとに、`want` の炉心の最大の数密度と、`scale`（炉心の最大の
    数密度）の `CRAM_ERROR` 倍の大きい方に対する比。どちらかが数密度を持たないか、形が違えば止める（比べずに差 0 と書かない）。比、その
    核種の番号（比べるノードがなければ `None`）、その核種の `want` の炉心の最大の数密度を返す"""
    if got is None or want is None:
        raise SystemExit(f"{what}: the densities are missing")
    shape = lambda d: [[len(n) for n in a] for a in d]
    if shape(got) != shape(want):
        raise SystemExit(f"{what}: the densities differ in shape")
    nodes = [(g, w) for ga, wa in zip(got, want) for g, w in zip(ga, wa) if w]
    if not nodes:
        return 0.0, None, 0.0
    count = len(nodes[0][1])
    floor = CRAM_ERROR * scale
    peak = [max(abs(w[i]) for _, w in nodes) for i in range(count)]
    largest = [max(floor, p) for p in peak]
    if not count:
        return 0.0, None, 0.0
    # 同じ比なら先の核種。差がすべて 0 なら、差の最大の核種はない
    pairs = [
        (abs(g[i] - w[i]) / (largest[i] or 1.0), i)
        for g, w in nodes
        for i in range(count)
    ]
    ratio = max(r for r, _ in pairs)
    if ratio == 0.0:
        return 0.0, None, 0.0
    i = min(i for r, i in pairs if r == ratio)
    return ratio, i, peak[i]


def compare_cut(name: str, full: dict, first: dict, second: dict) -> dict:
    """前の半分 `first` と後の半分 `second` を、切らない燃焼 `full` と比べた要約"""
    check_same_core(f"{name} (first half)", first, full)
    check_same_core(f"{name} (second half)", second, full)
    cut = first["points"][-1]
    start = second["points"][0]
    if start["burnups"] != cut["burnups"] or start.get("densities") != cut.get(
        "densities"
    ):
        raise SystemExit(
            f"{name}: the second half does not start from the last point of the first half"
        )
    n = len(first["points"]) - 1
    paired = full["points"][n:]
    if len(paired) != len(second["points"]):
        raise SystemExit(
            f"{name}: the second half has {len(second['points'])} points for the {len(paired)} uncut points after the cut"
        )
    points = []
    for p, q in zip(second["points"], paired):
        efpd, days = p["efpd"] + cut["efpd"], p["days"] + cut["days"]
        if not (close(efpd, q["efpd"]) and close(days, q["days"])):
            raise SystemExit(
                f"{name}: the point at {p['efpd']} EFPD, {p['days']} days after the cut is not the uncut point at "
                f"{q['efpd']} EFPD, {q['days']} days"
            )
        points.append(
            {
                "efpd": q["efpd"],
                "days": q["days"],
                "state": p["state"],
                "k": p["k"],
                "k_full": q["k"],
                "rho_pcm": rho_pcm(q["k"], p["k"]),
                "boron": p["boron"],
                "boron_full": q["boron"],
                "boron_diff": p["boron"] - q["boron"],
                "burnup_max_diff": max(
                    abs(a - b)
                    for ra, rb in zip(p["burnups"], q["burnups"])
                    for a, b in zip(ra, rb)
                ),
                "density_max_relative": relative_density_differences(
                    f"{name} at {q['efpd']} EFPD",
                    p.get("densities"),
                    q.get("densities"),
                    largest_density(q.get("densities") or []),
                ),
            }
        )
    auxiliary = largest_density_difference(
        f"{name} auxiliary",
        second.get("auxiliary_densities"),
        full.get("auxiliary_densities"),
        largest_density(
            full["points"][-1].get("densities") or [],
            full.get("auxiliary_densities") or [],
        ),
    )
    return {
        "name": name,
        "cut": {"efpd": cut["efpd"], "days": cut["days"]},
        "first_k_max_diff": max(
            abs(a["k"] - b["k"]) for a, b in zip(first["points"], full["points"])
        ),
        "points": points,
        "auxiliary_max_relative": auxiliary[0],
        "auxiliary_max_nuclide": None
        if auxiliary[1] is None
        else {"name": full["auxiliary"][auxiliary[1]], "largest": auxiliary[2]},
    }


def compare_reload(previous: dict, reload: dict) -> dict:
    """前のサイクル `previous` の履歴から装荷を替えた燃焼 `reload` の要約"""
    check_format("previous", previous)
    check_format("reload", reload)
    for key in ["axial_mesh", "nuclides", "auxiliary"]:
        if reload.get(key) != previous.get(key):
            raise SystemExit(f"reload: the {key} differs from the previous cycle")
    last = previous["points"][-1]
    start = reload["points"][0]
    inherited, fresh = [], []
    for a, id in enumerate(reload["assemblies"]):
        if id in previous["assemblies"]:
            b = previous["assemblies"].index(id)
            if reload["heavy_metal"][a] != previous["heavy_metal"][b]:
                raise SystemExit(
                    f"reload: assembly {id} has other heavy metal than in the previous cycle"
                )
            # ミクロ燃焼の履歴は、引き継いだ追う核種の数密度も前のサイクルの最後の点と同じ（ミクロ燃焼でなければどちらも持たない）
            densities = (start.get("densities"), last.get("densities"))
            if (densities[0] is None) != (densities[1] is None):
                raise SystemExit("reload: the densities are missing on one side")
            same = start["burnups"][a] == last["burnups"][b] and (
                densities[0] is None or densities[0][a] == densities[1][b]
            )
            if not same:
                raise SystemExit(
                    f"reload: assembly {id} does not start from its state at the last point of the previous cycle"
                )
            inherited.append(id)
        else:
            if any(x != 0.0 for x in start["burnups"][a]):
                raise SystemExit(
                    f"reload: the fresh assembly {id} does not start at burnup 0"
                )
            fresh.append(id)
    return {
        "assemblies": reload["assemblies"],
        "inherited": inherited,
        "fresh": fresh,
        "discharged": [
            a for a in previous["assemblies"] if a not in reload["assemblies"]
        ],
        "previous": {"efpd": last["efpd"], "days": last["days"]},
        "points": [
            {
                "efpd": p["efpd"],
                "days": p["days"],
                "state": p["state"],
                "k": p["k"],
                "boron": p["boron"],
                "solves": p["solves"],
                "burnup": p["burnup"],
                "assembly_burnups": p["assembly_burnups"],
            }
            for p in reload["points"]
        ],
    }


def summarize(
    full: dict, cuts: dict[str, tuple[dict, dict]], reload: dict | None
) -> dict:
    """要約（モジュールの説明の `--json`）。`cuts` は名前 → (前の半分, 後の半分)"""
    check_format("uncut", full)
    return {
        "format": "xlc-compare-core-restart",
        "cuts": [compare_cut(name, full, a, b) for name, (a, b) in cuts.items()],
        "reload": None if reload is None else compare_reload(full, reload),
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    lines = [
        "inputs: " + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
    ]
    for cut in summary["cuts"]:
        points = cut["points"]
        rho = [abs(p["rho_pcm"]) for p in points]
        boron = [abs(p["boron_diff"]) for p in points]
        n = cut["auxiliary_max_nuclide"]
        aux = "" if n is None else f"（{n['name']}、炉心の最大 {n['largest']:.3e}）"
        lines.append(
            f"切る点 {cut['name']}（{cut['cut']['efpd']:.4f} EFPD、{cut['cut']['days']:.4f} 日）: 後の半分の点 {len(points)}、"
            f"|Δρ| の最大 {max(rho):.3e} pcm、|ホウ素の差| の最大 {max(boron):.3e} ppm、"
            f"ノードの燃焼度の差の最大 {max(p['burnup_max_diff'] for p in points):.3e} MWd/kg、"
            f"追う核種の数密度の差の最大 {max(p['density_max_relative'] for p in points):.3e}、"
            f"最後の点の補助の核種 {cut['auxiliary_max_relative']:.3e}{aux}、前の半分の k の差の最大 {cut['first_k_max_diff']:.3e}"
        )
        for p in points:
            lines.append(
                f"  {p['efpd']:.4f} EFPD {p['days']:.4f} 日 {p['state']}: k {p['k']:.8f}（切らない {p['k_full']:.8f}）、"
                f"Δρ {p['rho_pcm']:+.3e} pcm、ホウ素 {p['boron']:.4f}（差 {p['boron_diff']:+.3e}） ppm"
            )
    reload = summary.get("reload")
    if reload is not None:
        lines.append(
            f"装荷替え: 引き継いだ {', '.join(reload['inherited'])}、新しい {', '.join(reload['fresh'])}、"
            f"取り出した {', '.join(reload['discharged'])}（前のサイクルの最後の点 {reload['previous']['efpd']:.4f} EFPD、"
            f"{reload['previous']['days']:.4f} 日）"
        )
        lines.append(
            "EFPD | 日数 | 状態 | k | ホウ素 (ppm) | 炉心の燃焼度 | "
            + " | ".join(reload["assemblies"])
        )
        for p in reload["points"]:
            lines.append(
                f"{p['efpd']:.4f} | {p['days']:.4f} | {p['state']} | {p['k']:.8f} | {p['boron']:.4f} | {p['burnup']:.4f} | "
                + " | ".join(f"{b:.4f}" for b in p["assembly_burnups"])
            )
    return lines


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes())


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
        return
    if len(args) < 3 or args[0] != "--json":
        raise SystemExit(__doc__)
    out, full_path, rest = Path(args[1]), Path(args[2]), args[3:]
    reload_path = None
    if "--reload" in rest:
        i = rest.index("--reload")
        if i + 2 != len(rest):
            raise SystemExit("give --reload <reload history> last")
        reload_path = Path(rest[i + 1])
        rest = rest[:i]
    pairs = {}
    for a in rest:
        name, sep, paths = a.partition("=")
        first, comma, second = paths.partition(",")
        if not (sep and comma and name and first and second):
            raise SystemExit(
                f"give a cut as <name>=<first half>,<second half>, got {a}"
            )
        if name in pairs:
            raise SystemExit(f"the cut {name} appears twice")
        pairs[name] = (Path(first), Path(second))
    files = [full_path, *(p for pair in pairs.values() for p in pair)]
    if reload_path is not None:
        files.append(reload_path)
    summary = {
        # 読んだファイルのパス（実行した cwd から）と md5
        "inputs": [
            {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
            for f in files
        ],
        **summarize(
            read(full_path),
            {n: (read(a), read(b)) for n, (a, b) in pairs.items()},
            None if reload_path is None else read(reload_path),
        ),
    }
    out.write_text(json.dumps(summary, indent=1) + "\n")
    print("\n".join(report(summary)))


if __name__ == "__main__":
    main()
