"""xlc の格子計算の燃焼の履歴（lattice の `[depletion]` の出力、ADR-0173）と、集合体の OpenMC の燃焼
（openmc_depletion_assembly.py extract。ADR-0178）を比べる。

使い方: python3 compare_depletion_assembly.py --json <要約の出力先> <xlc の履歴（msgpack）> <OpenMC の結果（msgpack）>
        python3 compare_depletion_assembly.py --report <要約の JSON>

- `--json`: 輸送する核種の集合、燃焼度の点、燃焼する物質（種類と体積）の並びが同じことを確かめてから、燃焼点ごとの k∞ の差
  (k_xlc − k_OpenMC) × 1e5 pcm（OpenMC の標準偏差つき）とその 0 の点からの変化、物質の種類ごとの主要な核種の数密度の相対差
  (N_xlc − N_OpenMC) / N_OpenMC（種類の領域の体積の重みの平均）、最後の点の種類ごとの領域の間のばらつき、OpenMC のピン出力の
  統計誤差の rms（対称操作の組で平均したピン出力の、標準偏差 / 平均の二乗平均平方根。xlc は対称な格子を 1/8 で解くので組で平均して比べる。
  ADR-0032 の目標との比較に使う）を書く。xlc の燃焼の領域と OpenMC の物質は、
  並びの順に対応させる（種類と体積で確かめる）。xlc の燃焼の履歴はピン出力を持たない（ADR-0173）ので、ピン出力の比較はまだない。
  `inputs` は読んだファイルのパス（段階の cwd から）と md5。時間の内訳（xlc の履歴の `seconds`、OpenMC の輸送ごとの時間と CRAM の
  時間）もコピーする。AWS の機種はまちまちなので、OpenMC の時間は xlc との時間の比較に使わない（ADR-0174）
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import msgpack

FUEL_NUCLIDES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"]
# 燃焼する物質の種類ごとに比べる核種。Gd 燃料のリング（gd0, gd1, ...）は Gd-155 と Gd-157、
# 被覆（IFBA）、Pyrex、WABA の吸収体は B-10
GD_NUCLIDES = ["Gd-155", "Gd-157"]
ABSORBER_NUCLIDES = ["B-10"]
# 燃焼度と体積の一致の許容（燃焼度は入力の同じ値、体積は xlc の格子の円の面積と OpenMC の面積）
BURNUP_TOLERANCE = 1e-9
VOLUME_TOLERANCE = 1e-6


def nuclides_of(kind: str) -> list[str]:
    """物質の種類（xlc の入力の物質の名前）の比べる核種。"""
    if kind == "fuel":
        return FUEL_NUCLIDES
    if re.fullmatch(r"gd\d+", kind):
        return GD_NUCLIDES
    if kind in ("ifba", "pyrex", "waba"):
        return ABSORBER_NUCLIDES
    raise SystemExit(f"unknown kind of burnable material {kind}")


def openmc_name(name: str) -> str:
    """xlc の核種名（`Am-242m`）を OpenMC の名前（`Am242_m1`）にする。"""
    m = re.fullmatch(r"([A-Z][a-z]?)-(\d+)(m(\d*))?", name)
    if m is None:
        raise SystemExit(f"unexpected nuclide name {name}")
    symbol, mass, isomer, level = m.groups()
    return f"{symbol}{mass}" + (f"_m{level or 1}" if isomer else "")


def relative(ours: float, theirs: float) -> float | None:
    """相対差。OpenMC が 0 なら比べない（None）"""
    return None if theirs == 0.0 else ours / theirs - 1.0


def check_same_run(xlc: dict, omc: dict) -> None:
    if (
        xlc["format"] != "xlc-depletion-history"
        or omc["format"] != "openmc-depletion-assembly"
    ):
        raise SystemExit("unexpected formats")
    ours = {openmc_name(n) for n in xlc["transportable"]}
    theirs = set(omc["transportable"])
    if ours != theirs:
        raise SystemExit(
            f"transported nuclides differ: xlc only {sorted(ours - theirs)}, "
            f"OpenMC only {sorted(theirs - ours)}"
        )
    burnup = [p["burnup"] for p in xlc["points"]]
    if len(burnup) != len(omc["burnup"]) or any(
        abs(a - b) > BURNUP_TOLERANCE * max(b, 1.0)
        for a, b in zip(burnup, omc["burnup"])
    ):
        raise SystemExit(f"burnup points differ: {burnup} vs {omc['burnup']}")
    regions, materials = xlc["regions"], omc["regions"]
    if len(regions) != len(materials):
        raise SystemExit(
            f"{len(regions)} xlc regions for {len(materials)} OpenMC materials"
        )
    for i, (r, m) in enumerate(zip(regions, materials)):
        if r["material"] != m["material"]:
            raise SystemExit(f"region {i}: {r['material']} vs {m['material']}")
        if abs(r["volume"] - m["volume"]) > VOLUME_TOLERANCE * m["volume"]:
            raise SystemExit(
                f"region {i} ({r['material']}): volume {r['volume']} vs {m['volume']}"
            )


def orbit_std_rms(mean: list[list[float]], std: list[list[float]]) -> float:
    """ピン出力（[行][列]）を対称操作の組（中央からの距離の対）で平均したときの統計誤差の rms。組の平均は Σ 平均 / 個数、その誤差は
    sqrt(Σ σ²) / 個数（組のピンの誤差を独立とする）で、平均が正の組の（誤差 / 平均）の二乗平均平方根。"""
    n = len(mean)
    centre = (n - 1) / 2
    orbits: dict[tuple[float, float], list[tuple[float, float]]] = {}
    for i in range(n):
        for j in range(n):
            key = tuple(sorted((abs(i - centre), abs(j - centre)), reverse=True))
            orbits.setdefault(key, []).append((mean[i][j], std[i][j]))
    relative_errors = []
    for members in orbits.values():
        count = len(members)
        average = sum(m for m, _ in members) / count
        if average > 0:
            error = math.sqrt(sum(s * s for _, s in members)) / count
            relative_errors.append(error / average)
    return math.sqrt(sum(e * e for e in relative_errors) / len(relative_errors))


def summarize(xlc: dict, omc: dict) -> dict:
    """2 つの履歴を比べた要約（モジュールの説明の `--json`）。"""
    check_same_run(xlc, omc)
    burnup = [p["burnup"] for p in xlc["points"]]
    index = {name: i for i, name in enumerate(xlc["nuclides"])}
    regions, materials = xlc["regions"], omc["regions"]
    kinds = list(dict.fromkeys(r["material"] for r in regions))

    def difference(r: int, point: int, name: str) -> float | None:
        return relative(
            xlc["points"][point]["densities"][r][index[name]],
            materials[r]["densities"][openmc_name(name)][point],
        )

    first = (xlc["points"][0]["k_inf"] - omc["k"][0]) * 1e5
    k = []
    for i, point in enumerate(xlc["points"]):
        diff = (point["k_inf"] - omc["k"][i]) * 1e5
        k.append(
            {
                "xlc": point["k_inf"],
                "openmc": omc["k"][i],
                "openmc_std_pcm": omc["k-std"][i] * 1e5,
                "diff_pcm": diff,
                "change_pcm": diff - first,
            }
        )
    by_kind: dict[str, dict[str, list]] = {}
    spread: dict[str, dict[str, dict]] = {}
    last = len(burnup) - 1
    for kind in kinds:
        members = [r for r, region in enumerate(regions) if region["material"] == kind]
        total = sum(regions[r]["volume"] for r in members)
        by_kind[kind], spread[kind] = {}, {}
        for name in nuclides_of(kind):
            series = []
            for point in range(len(burnup)):
                ours = sum(
                    xlc["points"][point]["densities"][r][index[name]]
                    * regions[r]["volume"]
                    for r in members
                )
                theirs = sum(
                    materials[r]["densities"][openmc_name(name)][point]
                    * regions[r]["volume"]
                    for r in members
                )
                series.append(relative(ours / total, theirs / total))
            by_kind[kind][name] = series
            diffs = [d for r in members if (d := difference(r, last, name)) is not None]
            spread[kind][name] = {
                "regions": len(diffs),
                "max_abs": max((abs(d) for d in diffs), default=None),
                "rms": math.sqrt(sum(d * d for d in diffs) / len(diffs))
                if diffs
                else None,
            }
    pin_std_rms = [
        orbit_std_rms(mean, std) for mean, std in zip(omc["pins"], omc["pins-std"])
    ]
    transport = omc["transport-seconds"]
    seconds = {
        "xlc": {**xlc["seconds"], "solves": xlc["solves"]},
        "openmc": {
            "transport": sum(transport),
            "fastest": min(transport),
            "slowest": max(transport),
            "solves": len(transport),
            "cram": sum(omc["cram-seconds"]),
            "total": omc["total-seconds"],
        },
    }
    return {
        "problem": omc["problem"],
        "burnup": burnup,
        "transported": len(omc["transportable"]),
        "regions": len(regions),
        "seconds": seconds,
        "k": k,
        "kinds": by_kind,
        "spread": spread,
        "pin_std_rms": pin_std_rms,
    }


def percent(x: float | None) -> str:
    return "–" if x is None else f"{x * 100:+.2f}"


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the stage")
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        (
            f"問題 {summary['problem']}、輸送する核種: 両方で同じ {summary['transported']} 核種、"
            f"燃焼する物質 {summary['regions']}"
        ),
        "燃焼度 (MWd/kg) | k xlc | k OpenMC | 差 (pcm) | 0 からの変化 (pcm) | OpenMC のピン出力の統計誤差 rms (%)",
    ]
    for b, k, s in zip(summary["burnup"], summary["k"], summary["pin_std_rms"]):
        lines.append(
            f"{b:g} | {k['xlc']:.6f} | {k['openmc']:.6f} ± {k['openmc_std_pcm'] / 1e5:.6f} | "
            f"{k['diff_pcm']:+.0f} ± {k['openmc_std_pcm']:.0f} | {k['change_pcm']:+.0f} | {s * 100:.3f}"
        )
    for kind, table in summary["kinds"].items():
        names = list(table)
        lines.append(f"{kind} の数密度の相対差 (%): 燃焼度 | " + " | ".join(names))
        for i, b in enumerate(summary["burnup"]):
            lines.append(f"{b:g} | " + " | ".join(percent(table[n][i]) for n in names))
        lines.append(
            f"{kind} の最後の点の領域間のばらつき (%): 核種 | 領域数 | 最大の絶対値 | rms"
        )
        for n in names:
            sp = summary["spread"][kind][n]
            lines.append(
                f"{n} | {sp['regions']} | "
                f"{percent(sp['max_abs']).lstrip('+')} | {percent(sp['rms']).lstrip('+')}"
            )
    x, o = summary["seconds"]["xlc"], summary["seconds"]["openmc"]
    lines.append(
        f"時間 xlc（輸送 {x['solves']} 回）: 自己遮蔽 {x['shielding']:.1f} 秒、MOC {x['moc']:.1f} 秒、"
        f"反応率 {x['rates']:.2f} 秒、燃焼行列 {x['matrix']:.2f} 秒、CRAM {x['cram']:.2f} 秒"
    )
    lines.append(
        f"時間 OpenMC（輸送 {o['solves']} 回、AWS の機種による。xlc との比較に使わない）: 全体 {o['total']:.0f} 秒、"
        f"輸送 {o['transport']:.0f} 秒（1 回 {o['fastest']:.0f}〜{o['slowest']:.0f} 秒）、CRAM {o['cram']:.1f} 秒"
    )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 4 and args[0] == "--json":
        files = [Path(args[2]), Path(args[3])]
        xlc, omc = (msgpack.unpackb(f.read_bytes()) for f in files)
        summary = {
            "format": "xlc-compare-depletion-assembly",
            # 読んだファイルのパス（段階の cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(xlc, omc),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
