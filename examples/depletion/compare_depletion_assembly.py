"""xlc の格子計算の燃焼の履歴（lattice の `[depletion]` の出力、ADR-0173）と、集合体の OpenMC の燃焼
（openmc_depletion_assembly.py extract。ADR-0178）を比べる。

使い方: python3 compare_depletion_assembly.py --json <要約の出力先> <xlc の履歴（msgpack）> <OpenMC の結果（msgpack）>
        python3 compare_depletion_assembly.py --report <要約の JSON>

- `--json`: 輸送する核種の集合、燃焼度の点、燃焼する物質（種類と体積）の並びが同じことを確かめてから、燃焼点ごとの k∞ の差
  (k_xlc − k_OpenMC) × 1e5 pcm（OpenMC の標準偏差つき）とその 0 の点からの変化、物質の種類ごとの主要な核種の数密度の相対差
  (N_xlc − N_OpenMC) / N_OpenMC（種類の領域の体積の重みの平均）、最後の点の種類ごとの領域の間のばらつき、OpenMC のピン出力の
  統計誤差の rms（対称操作の組で平均したピン出力の、標準偏差 / 平均の二乗平均平方根。xlc は対称な格子を 1/8 で解くので組で平均して比べる。
  ADR-0032 の目標との比較に使う）、燃焼点ごとのピン出力の差（`pins`。xlc と OpenMC のピン出力をそれぞれ組で平均した相対差
  (xlc / OpenMC − 1) の、核分裂する物質（燃料と Gd 燃料のリング）のある組の rms と、絶対値が最大の組の差、位置（組のうち上の行から
  最初のセルの [行, 列]、0 から）、その組の OpenMC の統計誤差、燃焼する吸収体だけのセルの OpenMC の核分裂の絶対値の最大
  `absorber_max`。ADR-0181）を書く。OpenMC のピン出力は、核分裂する物質のあるセルの平均で 1 に割り直してから比べる。xlc の燃焼の領域と OpenMC の物質は、並びの順に対応させる（種類と体積で確かめる）。
  `inputs` は読んだファイルのパス（ステージの cwd から）と md5。時間の内訳（xlc の履歴の `seconds`、OpenMC の輸送ごとの時間と CRAM の
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
# 燃焼する吸収体のセルの OpenMC の核分裂（燃料のピンの平均 1 に対する値）の上限。CRAM の丸めの誤差でできるアクチニドの核分裂は
# 2F と 2N で 5e-16 以下で、燃料の取り違えなら 1 の桁になる。その間の、どちらからも桁の離れた値（`pin_differences`）
ABSORBER_FISSION_LIMIT = 1e-6


def compared_openmc_names() -> set[str]:
    """比べる核種の OpenMC の名前（OpenMC の結果が数密度を書く核種）。"""
    return {openmc_name(n) for n in FUEL_NUCLIDES + GD_NUCLIDES + ABSORBER_NUCLIDES}


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


def orbit_of(i: int, j: int, n: int) -> tuple[int, int]:
    """n × n のセル（行 i、列 j）の対称操作の組。中央からの距離の対を大きい順に並べ、整数にする（n が偶数なら半分の距離を 2 倍にする）。"""
    a, b = sorted((abs(2 * i - (n - 1)), abs(2 * j - (n - 1))), reverse=True)
    return (a, b)


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
        # xlc のセルの番号は水のギャップの 1 周を含む（n + 2 列）。代表のセルの組が、OpenMC の物質の組と同じこと（体積が同じ別の組を
        # 取り違えない。軸上と対角線上の組は同じセル数）
        n = len(omc["pins"][0])
        row, column = divmod(r["cell"], n + 2)
        if orbit_of(row - 1, column - 1, n) != tuple(2 * x for x in m["orbit"]):
            raise SystemExit(
                f"region {i}: cell {r['cell']} is not in the orbit {m['orbit']}"
            )
        if abs(r["volume"] - m["volume"]) > VOLUME_TOLERANCE * m["volume"]:
            raise SystemExit(
                f"region {i} ({r['material']}): volume {r['volume']} vs {m['volume']}"
            )


def orbits(n: int) -> list[list[tuple[int, int]]]:
    """n × n のセルの対称操作の組（中央からの距離の対が同じセル）。組の中のセルは上の行から、各行は左から。"""
    found: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for i in range(n):
        for j in range(n):
            found.setdefault(orbit_of(i, j, n), []).append((i, j))
    return list(found.values())


def orbit_mean(values: list[list[float]], members: list[tuple[int, int]]) -> float:
    return sum(values[i][j] for i, j in members) / len(members)


def orbit_error(std: list[list[float]], members: list[tuple[int, int]]) -> float:
    """組の平均の統計誤差 sqrt(Σ σ²) / 個数（組のピンの誤差を独立とする）。"""
    return math.sqrt(sum(std[i][j] ** 2 for i, j in members)) / len(members)


def is_fissile(kind: str) -> bool:
    """燃焼する物質の種類が核分裂する物質（燃料と Gd 燃料のリング）か。IFBA の被覆は燃料のピンの中にあり、Pyrex と WABA は案内管の中"""
    return kind == "fuel" or re.fullmatch(r"gd\d+", kind) is not None


def cell_kinds(materials: list[dict]) -> tuple[set, set]:
    """OpenMC の燃焼する物質（`regions`）から、核分裂する物質のある組と、燃焼する吸収体（Pyrex、WABA）だけの組（`orbit_of` の形）。"""
    fissile = {
        tuple(2 * x for x in m["orbit"]) for m in materials if is_fissile(m["material"])
    }
    if not fissile:
        raise SystemExit(
            f"no fuel among the burnable materials {sorted({m['material'] for m in materials})}"
        )
    absorbers = {tuple(2 * x for x in m["orbit"]) for m in materials} - fissile
    return fissile, absorbers


def fuel_pin_mean(rows: list[list[float]], fissile: set) -> float:
    """核分裂する物質のあるセルのピン出力の平均。OpenMC の抽出（openmc_run.py の fuel_mean）は値が正のセルの平均で割るので、燃焼する
    吸収体のセルに丸めの誤差の核分裂（`pin_differences`）があると、そのセルも数えた平均で割られている。比べる前にこれで割り直す"""
    n = len(rows)
    values = [
        rows[i][j] for i in range(n) for j in range(n) if orbit_of(i, j, n) in fissile
    ]
    return sum(values) / len(values)


def orbit_std_rms(
    mean: list[list[float]], std: list[list[float]], fissile: set
) -> float:
    """ピン出力（[行][列]）を対称操作の組で平均したときの、核分裂する物質のある組の（統計誤差 / 平均）の二乗平均平方根。"""
    n = len(mean)
    relative_errors = [
        orbit_error(std, members) / orbit_mean(mean, members)
        for members in orbits(n)
        if orbit_of(*members[0], n) in fissile
    ]
    return math.sqrt(sum(e * e for e in relative_errors) / len(relative_errors))


def pin_differences(
    ours: list[list[float]],
    mean: list[list[float]],
    std: list[list[float]],
    fissile: set,
    absorbers: set,
) -> dict:
    """xlc のピン出力 `ours` と OpenMC の `mean`（統計誤差 `std`）を、それぞれ対称操作の組で平均して比べる（ADR-0181）。核分裂する
    物質のある組（`fissile`）の相対差の rms と、絶対値が最大の組の差、位置、統計誤差。核分裂のあるセルが物質の配置と食い違えば止める。

    燃焼する吸収体だけの組（`absorbers`）の OpenMC の核分裂は比べずに、その最大を `absorber_max` に残す。OpenMC の燃焼はチェーンの
    全核種を CRAM で解くので、アクチニドのない吸収体にも、物質の原子の数に対して丸めの誤差の大きさのアクチニド（負の値も含む）が
    でき、核分裂の tally が値を持つ（2F の Pyrex で、平均 1 のピン出力に対して 1e-24〜5e-16。負の値もありうる）。どのセルを比べるかは
    大きさで 0 とみなすしきい値でなく、物質の配置で決める。ただし吸収体のセルの値が `ABSORBER_FISSION_LIMIT` を超えれば、丸めの
    誤差でなく配置の取り違え（吸収体の位置に燃料がある、など）として止める"""
    n = len(mean)
    if (
        len(ours) != n
        or len(std) != n
        or any(len(row) != n for row in [*ours, *mean, *std])
    ):
        raise SystemExit(
            f"pin maps of other shapes: xlc {[len(r) for r in ours]}, OpenMC {[len(r) for r in mean]}, "
            f"OpenMC errors {[len(r) for r in std]}"
        )
    absorber_max = 0.0
    for i in range(n):
        for j in range(n):
            orbit = orbit_of(i, j, n)
            if orbit in fissile:
                if not (ours[i][j] > 0 and mean[i][j] > 0):
                    raise SystemExit(
                        f"pin ({i}, {j}): no fission in a fuel pin (xlc {ours[i][j]}, OpenMC {mean[i][j]})"
                    )
            elif ours[i][j] != 0 or (orbit not in absorbers and mean[i][j] != 0):
                raise SystemExit(
                    f"pin ({i}, {j}): fission in a cell without fuel (xlc {ours[i][j]}, OpenMC {mean[i][j]})"
                )
            elif abs(mean[i][j]) > ABSORBER_FISSION_LIMIT:
                raise SystemExit(
                    f"pin ({i}, {j}): OpenMC fission {mean[i][j]} in a burnt absorber is not a round-off"
                )
            else:
                absorber_max = max(absorber_max, abs(mean[i][j]))
    diffs = []
    for members in orbits(n):
        theirs = orbit_mean(mean, members)
        values = [ours[i][j] for i, j in members]
        # xlc は対称操作で 1/8 を解いて展開するので、組の中の値は同じ（丸めは 1e-9 よりずっと小さい）。違えば格子が対称でないか入力の誤りで、
        # 平均で隠さない
        if max(values) - min(values) > 1e-9 * max(abs(v) for v in values):
            raise SystemExit(
                f"xlc pin powers differ within the symmetry orbit of {members[0]}: {values}"
            )
        if orbit_of(*members[0], n) in fissile:
            diffs.append(
                (
                    orbit_mean(ours, members) / theirs - 1.0,
                    list(members[0]),
                    orbit_error(std, members) / theirs,
                )
            )
    worst, at, error = max(diffs, key=lambda d: abs(d[0]))
    return {
        "orbits": len(diffs),
        "rms": math.sqrt(sum(d * d for d, _, _ in diffs) / len(diffs)),
        "max": worst,
        "at": at,
        "max_std": error,
        "absorber_max": absorber_max,
    }


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
    if any("pin_powers" not in p for p in xlc["points"]):
        raise SystemExit(
            "the xlc history has no pin powers (ADR-0181); rerun the lattice stage"
        )
    if not len(xlc["points"]) == len(omc["pins"]) == len(omc["pins-std"]):
        raise SystemExit(
            f"{len(xlc['points'])} points with {len(omc['pins'])} OpenMC pin maps and "
            f"{len(omc['pins-std'])} error maps"
        )
    fissile, absorbers = cell_kinds(materials)
    openmc_pins = []
    for mean, std in zip(omc["pins"], omc["pins-std"]):
        scale = fuel_pin_mean(mean, fissile)
        openmc_pins.append(
            (
                [[v / scale for v in row] for row in mean],
                [[v / scale for v in row] for row in std],
            )
        )
    pins = [
        pin_differences(p["pin_powers"], mean, std, fissile, absorbers)
        for p, (mean, std) in zip(xlc["points"], openmc_pins)
    ]
    # 形はピン出力の比較で確かめてある
    pin_std_rms = [orbit_std_rms(mean, std, fissile) for mean, std in openmc_pins]
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
        "pins": pins,
    }


def percent(x: float | None) -> str:
    return "–" if x is None else f"{x * 100:+.2f}"


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the stage")
    if "pins" not in summary:
        raise SystemExit(
            "the JSON has no pin power differences (ADR-0181); rerun the stage"
        )
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
    if (
        not len(summary["burnup"])
        == len(summary["pins"])
        == len(summary["pin_std_rms"])
    ):
        raise SystemExit(
            f"{len(summary['burnup'])} burnup points with {len(summary['pins'])} pin rows and "
            f"{len(summary['pin_std_rms'])} error rows; rerun the stage"
        )
    lines.append(
        "ピン出力の差（対称操作の組、%）: 燃焼度 | 組 | rms | 最大 ± OpenMC の統計誤差 | 位置（行, 列。1 から） | OpenMC の統計誤差 rms"
    )
    for b, p, s in zip(summary["burnup"], summary["pins"], summary["pin_std_rms"]):
        row, column = p["at"]
        lines.append(
            f"{b:g} | {p['orbits']} | {p['rms'] * 100:.3f} | {p['max'] * 100:+.3f} ± {p['max_std'] * 100:.3f} | "
            f"({row + 1}, {column + 1}) | {s * 100:.3f}"
        )
    if any("absorber_max" not in p for p in summary["pins"]):
        raise SystemExit("the JSON has no fission of burnt absorbers; rerun the stage")
    lines.append(
        "燃焼する吸収体だけのセルの OpenMC の核分裂（ピン出力。比べない）の最大: "
        f"{max(p['absorber_max'] for p in summary['pins']):.1e}"
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
            # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
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
