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
    絶対値の最大（核種ごとの、切らない燃焼のその点の炉心の最大の数密度に対する比）。最後の点の補助の核種の数密度の差も同じ比で。前の半分の
    点と切らない燃焼の初めの点の k の差の絶対値の最大も書く（同じ計算なので 0 の見込み）
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


def relative_density_differences(got: list, want: list) -> float:
    """ノードの数密度 [集合体][軸方向のノード][核種] の差の絶対値の最大（核種ごとの `want` の炉心の最大の数密度に対する比。最大が 0 の核種は
    差の絶対値）"""
    nodes = [(g, w) for ga, wa in zip(got, want) for g, w in zip(ga, wa) if w]
    if not nodes:
        return 0.0
    count = len(nodes[0][1])
    largest = [max(abs(w[i]) for _, w in nodes) for i in range(count)]
    return max(
        (
            abs(g[i] - w[i]) / (largest[i] or 1.0)
            for g, w in nodes
            for i in range(count)
        ),
        default=0.0,
    )


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
                    p.get("densities") or [], q.get("densities") or []
                ),
            }
        )
    return {
        "name": name,
        "cut": {"efpd": cut["efpd"], "days": cut["days"]},
        "first_k_max_diff": max(
            abs(a["k"] - b["k"]) for a, b in zip(first["points"], full["points"])
        ),
        "points": points,
        "auxiliary_max_relative": relative_density_differences(
            second.get("auxiliary_densities") or [],
            full.get("auxiliary_densities") or [],
        ),
    }


def compare_reload(previous: dict, reload: dict) -> dict:
    """前のサイクル `previous` の履歴から装荷を替えた燃焼 `reload` の要約"""
    check_format("previous", previous)
    check_format("reload", reload)
    if reload["axial_mesh"] != previous["axial_mesh"]:
        raise SystemExit("reload: the axial mesh differs from the previous cycle")
    last = previous["points"][-1]
    start = reload["points"][0]
    inherited, fresh = [], []
    for a, id in enumerate(reload["assemblies"]):
        if id in previous["assemblies"]:
            b = previous["assemblies"].index(id)
            same = start["burnups"][a] == last["burnups"][b] and (
                start.get("densities") is None
                or start["densities"][a] == last["densities"][b]
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
        lines.append(
            f"切る点 {cut['name']}（{cut['cut']['efpd']:.4f} EFPD、{cut['cut']['days']:.4f} 日）: 後の半分の点 {len(points)}、"
            f"|Δρ| の最大 {max(rho):.3e} pcm、|ホウ素の差| の最大 {max(boron):.3e} ppm、"
            f"ノードの燃焼度の差の最大 {max(p['burnup_max_diff'] for p in points):.3e} MWd/kg、"
            f"追う核種の数密度の差の最大 {max(p['density_max_relative'] for p in points):.3e}、"
            f"最後の点の補助の核種 {cut['auxiliary_max_relative']:.3e}、前の半分の k の差の最大 {cut['first_k_max_diff']:.3e}"
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
