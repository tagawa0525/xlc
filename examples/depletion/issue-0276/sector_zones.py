"""#0276 の試行: 燃料のピンを方位の扇形の組ごとにも燃焼させる（ADR-0244 の `[depletion]` の `sectors`）効果を測る、2A の燃料の
ピンを並べた最小の格子の入力を書く。

使い方:
  python3 sector_zones.py write   # この置き場に下の 8 本の入力と、その最初の 2 点だけの short-*.toml を書く
  python3 sector_zones.py compare --at <燃焼度>,... <要約.json>
      # この置き場の 8 本の履歴（<名前>.msgpack）を PAIRS の組で比べ、すべての燃焼点の k∞、反応度、ピン出力の差と、扇形の
      # 組の間の組成の広がりを JSON に書く。#0275 のピンセルの径方向の効き（PIN_CELL の組）も同じ燃焼点で添える。--at は
      # 記録の表に出す燃焼度（MWd/kg）。表は examples/verification_tables.py が JSON から描く（ADR-0176）
  python3 sector_zones.py report <要約.json>   # 組ごとの Δk、Δρ、ピン出力の差の極値（すべての燃焼点から）

格子（どちらも外周は鏡面反射）:
- guide: 中央が 2A の案内管の 3 × 3。案内管の隣のピンは、辺で接するピン（辺のピンの組）と角で接するピン（角のピンの組）の 2 つの
  組になる。外周が鏡面反射なので、案内管が 3 ピッチごとに並ぶ無限の格子（17 × 17 の集合体の案内管と計装管の割合 25/289 より
  密な 1/9）
- gap: 2A の燃料のピンだけの 3 × 3 の周りに、2A の入力と同じ幅 0.04 cm の水の余白を付けたもの。外周が鏡面反射なので、3 × 3 の
  小さな集合体が 0.08 cm の水ギャップを挟んで並ぶ格子で、角のピンは集合体の角（2 辺が水ギャップ）、辺のピンは集合体の辺、中央の
  ピンは水ギャップに接しないピンを表す。案内管と同じ格子に置かないのは、案内管と水ギャップの効きを分けるため

ケース（格子ごとに 4 本）:
- rings1: 本番に向けた 2A の入力（../p9-product/lattice-2A.toml）と同じ燃焼の単位（燃料 1 つ、MOC の 5 リング）。方位に分けない
- rings5: 燃料の 5 リングを 1 本ずつ 5 つの単位にしたもの（#0275 の zones-5-p9 と同じ分け方）。方位に分けない
- rings5-sectors4、rings5-sectors8: rings5 の各単位を、さらに 4 つ（象限）、8 つの扇形の組に分けたもの
MOC の扇形は 4 本とも 8 つ（2A の入力と同じ）で、方位に分けるかだけが違う。径方向の効きは rings5 − rings1、方位の効きは
rings5-sectors* − rings5 で、どちらも同じ格子で比べる。

compare の量:
- Δk = k_他 − k_基準、Δρ = 1/k_基準 − 1/k_他（pcm）
- ピン出力（履歴の `pin_powers`。核分裂のあるセルの平均が 1）の相対差 k_他/k_基準 − 1 の、燃料のピン（基準の出力が 0 でない
  位置）の rms と、絶対値が最大の差（%）。rms は位置ごと（移り合うピンを重ねて数える）で、ADR-0032 のピン出力 rms と同じ数え方
- 扇形の組の広がり: 扇形に分けた履歴の、組の代表のセル（ピンの組）とリング（物質）ごとの、扇形の組の数密度の
  (最大 − 最小) / 体積の平均（%）。体積は履歴の領域の `volume`（組のピンの和）。扇形がすべて 1 つの組に移り合うピン（gap の
  中央）は広がりを持たない（null）。最も外のリング（fuel-4）は、最後の燃焼点の扇形の組ごとの平均に対する相対差（%）も書く。扇形の組の番号は
  履歴の `sector`（x 軸から反時計回りの扇形の番号のうち組で最も小さいもの。ADR-0244）、ピンの位置は組の最初のピンの
  （行（上から）、列（左から））で、gap は水の余白も行と列に数える（燃料のピンは 1〜3 の行と列）。要約の行は --at の燃焼点
  だけで、最大（`max`）はすべての燃焼点から

形状、組成、温度、減速材、燃焼度の点、出力密度、UFG の許容差、トラックは ../p9-product/lattice-2A.toml の燃料のピンと案内管と
公称の状態（ホウ素 650 ppm、減速材 580 K、燃料 900 K）に合わせる。
- short-<名前>.toml: 上の入力の燃焼度を最初の 2 点（輸送 5 回）にしたもの。全部を流す前の時間の見積もり
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
PELLET = 0.4096
RINGS = 5
FUEL = "{ U-234 = 6.11864e-06, U-235 = 7.18132e-04, U-236 = 3.29861e-06, U-238 = 2.21546e-02, O-16 = 4.57642e-02 }"
ZIRCALOY = """\
Cr-50 = 3.30121e-06
Cr-52 = 6.36606e-05
Cr-53 = 7.21860e-06
Cr-54 = 1.79686e-06
Fe-54 = 8.68307e-06
Fe-56 = 1.36306e-04
Fe-57 = 3.14789e-06
Fe-58 = 4.18926e-07
Zr-90 = 2.18865e-02
Zr-91 = 4.77292e-03
Zr-92 = 7.29551e-03
Zr-94 = 7.39335e-03
Zr-96 = 1.19110e-03
Sn-112 = 4.68066e-06
Sn-114 = 3.18478e-06
Sn-115 = 1.64064e-06
Sn-116 = 7.01616e-05
Sn-117 = 3.70592e-05
Sn-118 = 1.16872e-04
Sn-119 = 4.14504e-05
Sn-120 = 1.57212e-04
Sn-122 = 2.23417e-05
Sn-124 = 2.79392e-05
"""
BURNUP = """\
  0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5,
  4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
  14, 15, 16, 18, 20, 22, 24, 26, 28, 30,
  32, 34, 36, 38, 40, 44, 48, 52, 56, 60,
"""
SHORT_BURNUP = "  0.1, 0.25,\n"
LATTICES = {
    # 名前: (行、余白)
    "guide": ('["pin pin pin", "pin gt pin", "pin pin pin"]', ""),
    "gap": (
        '["pin pin pin", "pin pin pin", "pin pin pin"]',
        'gap = { width = 0.04, material = "water" }\n',
    ),
}
CASES = {
    # 名前: (径方向の単位の数、方位の扇形の組の数。None は分けない)
    "rings1": (1, None),
    "rings5": (RINGS, None),
    "rings5-sectors4": (RINGS, 4),
    "rings5-sectors8": (RINGS, 8),
}
GUIDE_TUBE = """\
# 案内管（中は減速材）
[cell.gt]
zones = [
  { material = "water", radius = 0.561, rings = 3 },
  { material = "zircaloy", radius = 0.602, rings = 1 },
  { material = "water", radius = 0.75, rings = 2 },
  { material = "water", rings = 1 },
]
sectors = 8

"""


def fuel_names(units: int) -> list[str]:
    return ["fuel"] if units == 1 else [f"fuel-{k}" for k in range(units)]


def input_text(
    lattice: str, units: int, sectors: int | None, output: str, burnup: str = BURNUP
) -> str:
    rows, gap = LATTICES[lattice]
    names = fuel_names(units)
    materials = "".join(
        f"[material.{n}]\ntemperature = 900\nnuclides = {FUEL}\n\n" for n in names
    )
    # 単位は MOC の 5 リングを内から等体積に分けたもの（xlc は帯の中を等体積に分けるので、単位の数によらず MOC の領域が同じ）
    fuel_zones = "".join(
        f'  {{ material = "{n}", radius = {PELLET * math.sqrt((k + 1) / units)!r}, rings = {RINGS // units} }},\n'
        for k, n in enumerate(names)
    )
    listed = ", ".join(f'"{n}"' for n in names)
    split = "" if sectors is None else f"sectors = {sectors}\n"
    guide = GUIDE_TUBE if lattice == "guide" else ""
    radial = "燃料 1 つ" if units == 1 else f"等体積の {units} つの単位"
    azimuthal = (
        "方位に分けない"
        if sectors is None
        else f"各単位を {sectors} つの扇形の組に分ける"
    )
    return f"""\
# #0276 の試行（sector_zones.py が書く。手で直さない）: 格子 {lattice} の 2A の燃料のピンを、{radial}で、{azimuthal}。
# p9-product/lattice-2A.toml の公称の状態のまま燃焼させる。
library = "../p2/depletion.msgpack"
ufg-tolerance = 1e-4
pitch = 1.26
azimuthal = 30
spacing = 0.03
{gap}rows = {rows}

{materials}[material.zircaloy]
temperature = 580

[material.zircaloy.nuclides]
{ZIRCALOY}
[material.water]
temperature = 580
water = {{ density = 0.711901, boron = 650, b10 = 0.199 }}

# 燃料のピン（3.1%）
[cell.pin]
zones = [
{fuel_zones}  {{ material = "void", radius = 0.418, rings = 1 }},
  {{ material = "zircaloy", radius = 0.475, rings = 1 }},
  {{ material = "water", radius = 0.63, rings = 3 }},
  {{ material = "water", radius = 0.75, rings = 2 }},
  {{ material = "water", rings = 1 }},
]
sectors = 8

{guide}[depletion]
chain = "../chain.msgpack"
power = 38.41
burnup = [
{burnup}]
materials = [{listed}]
{split}output = "{output}"
"""


def write() -> None:
    for lattice in LATTICES:
        for case, (units, sectors) in CASES.items():
            name = f"{lattice}-{case}"
            text = input_text(lattice, units, sectors, f"{name}.msgpack")
            (HERE / f"{name}.toml").write_text(text)
            short = input_text(
                lattice, units, sectors, f"short-{name}.msgpack", SHORT_BURNUP
            )
            (HERE / f"short-{name}.toml").write_text(short)


def same_burnup(a: float, b: float) -> bool:
    return abs(a - b) <= 1e-9 * max(1.0, abs(b))


def signed(x: float, digits: int) -> str:
    """符号つきの小数（記録の書き方に合わせ、負の符号は U+2212。丸めて 0 になる値に符号を付けない）"""
    text = f"{x:+.{digits}f}"
    if float(text) == 0:
        text = f"{0:.{digits}f}"
    return text.replace("-", "−")


def pin_power_difference(base: list[list[float]], other: list[list[float]]) -> dict:
    """燃料のピンの出力の相対差の rms と、絶対値が最大の差（%）"""
    if [len(r) for r in base] != [len(r) for r in other]:
        raise SystemExit("pin powers of other shapes")
    diffs = []
    for p_row, q_row in zip(base, other, strict=True):
        for p, q in zip(p_row, q_row, strict=True):
            if p == 0:
                if q != 0:
                    raise SystemExit("power in a pin without fuel in the base")
                continue
            diffs.append((q / p - 1) * 100)
    return {
        "rms_percent": math.sqrt(sum(d * d for d in diffs) / len(diffs)),
        "max_percent": max(diffs, key=abs),
    }


def sector_spread(values: list[float], volumes: list[float]) -> float | None:
    """扇形の組の値の (最大 − 最小) / 体積の平均（%）。組が 1 つか、平均が 0 のときは None"""
    if len(values) < 2:
        return None
    mean = sum(v * w for v, w in zip(values, volumes, strict=True)) / sum(volumes)
    return None if mean == 0 else (max(values) - min(values)) / mean * 100


def sector_rows(h: dict, nuclides: list[str]) -> list[dict]:
    """燃焼点ごとの、ピンの組とリングごとの扇形の組の広がりと、最も外のリングの扇形の組ごとの相対差"""
    regions = h["regions"]
    if any("sector" not in r for r in regions):
        raise SystemExit("history has regions not split by sector")
    index = {n: i for i, n in enumerate(h["nuclides"])}
    cells = sorted({r["cell"] for r in regions})
    rings = sorted({r["material"] for r in regions})
    # ピンの組とリングごとの、扇形の組の領域（扇形の番号の順）
    members = {
        (c, m): sorted(
            (i for i, r in enumerate(regions) if r["cell"] == c and r["material"] == m),
            key=lambda i: regions[i]["sector"],
        )
        for c in cells
        for m in rings
    }
    outer = rings[-1]
    rows = []
    for p in h["points"]:
        d = p["densities"]
        out = []
        for c in cells:
            spreads = {
                n: [
                    sector_spread(
                        [d[i][index[n]] for i in members[(c, m)]],
                        [regions[i]["volume"] for i in members[(c, m)]],
                    )
                    for m in rings
                ]
                for n in nuclides
            }
            ms = members[(c, outer)]
            volumes = [regions[i]["volume"] for i in ms]
            by_sector = {"sectors": [regions[i]["sector"] for i in ms]}
            for n in nuclides:
                values = [d[i][index[n]] for i in ms]
                mean = sum(v * w for v, w in zip(values, volumes)) / sum(volumes)
                by_sector[n] = [
                    None if mean == 0 else (v / mean - 1) * 100 for v in values
                ]
            out.append(
                {
                    "position": regions[ms[0]]["cells"][0],
                    "ring_spread_percent": spreads,
                    "outer_by_sector_percent": by_sector,
                }
            )
        rows.append({"burnup": p["burnup"], "cells": out})
    return rows


def spread_maxima(rows: list[dict], nuclides: list[str]) -> dict:
    """核種ごとの、すべての燃焼点、ピンの組、リングの中で最大の扇形の組の広がり"""
    maxima = {}
    for n in nuclides:
        x, burnup, position, ring = max(
            (x, row["burnup"], cell["position"], ring)
            for row in rows
            for cell in row["cells"]
            for ring, x in enumerate(cell["ring_spread_percent"][n])
            if x is not None
        )
        maxima[n] = {
            "percent": x,
            "burnup": burnup,
            "position": position,
            "ring": ring,
        }
    return maxima


def spread_summary(case: str, rows: list[dict], at: list[float]) -> dict:
    """要約に書く扇形の組の広がり: 行は --at の燃焼点だけで、最も外のリングの扇形の組ごとの差は最後の燃焼点だけ（JSON を
    小さく保つ）。最大はすべての燃焼点から"""
    nuclides = list(rows[0]["cells"][0]["ring_spread_percent"])
    last = rows[-1]
    return {
        "case": case,
        "max": spread_maxima(rows, nuclides),
        "rows": [
            {
                "burnup": r["burnup"],
                "cells": [
                    {
                        "position": c["position"],
                        "ring_spread_percent": c["ring_spread_percent"],
                    }
                    for c in r["cells"]
                ],
            }
            for r in rows
            if any(same_burnup(r["burnup"], x) for x in at)
        ],
        "end": {
            "burnup": last["burnup"],
            "cells": [
                {
                    "position": c["position"],
                    "outer_by_sector_percent": c["outer_by_sector_percent"],
                }
                for c in last["cells"]
            ],
        },
    }


def summary(path: Path) -> dict:
    h = msgpack.unpackb(path.read_bytes())
    if h["format"] != "xlc-depletion-history":
        raise SystemExit(f"{path}: not an xlc depletion history")
    return {
        "path": str(path),
        "history": h,
        "points": [
            {"burnup": p["burnup"], "k": p["k_inf"], "pins": p["pin_powers"]}
            for p in h["points"]
        ],
    }


def pair_rows(base: dict, other: dict) -> list[dict]:
    """すべての燃焼点の、基準に対する差（記録の表の行と極値は、この行から描く）"""
    if len(other["points"]) != len(base["points"]):
        raise SystemExit(f"{other['path']}: burnup points differ from {base['path']}")
    rows = []
    for p, q in zip(base["points"], other["points"]):
        if not same_burnup(q["burnup"], p["burnup"]):
            raise SystemExit(
                f"{other['path']}: burnup {q['burnup']} != {p['burnup']} of {base['path']}"
            )
        pins = pin_power_difference(p["pins"], q["pins"])
        rows.append(
            {
                "burnup": p["burnup"],
                "k_base": p["k"],
                "dk_pcm": (q["k"] - p["k"]) * 1e5,
                "drho_pcm": (1 / p["k"] - 1 / q["k"]) * 1e5,
                "pin_rms_percent": pins["rms_percent"],
                "pin_max_percent": pins["max_percent"],
            }
        )
    return rows


def extremes_line(rows: list[dict]) -> str:
    hi = max(rows, key=lambda r: r["dk_pcm"])
    lo = min(rows, key=lambda r: r["dk_pcm"])
    far = max(rows, key=lambda r: abs(r["drho_pcm"]))
    rms = max(rows, key=lambda r: r["pin_rms_percent"])
    pin = max(rows, key=lambda r: abs(r["pin_max_percent"]))
    return (
        f"Δk の最大 {signed(hi['dk_pcm'], 2)} pcm（{hi['burnup']:g} MWd/kg）、"
        f"最小 {signed(lo['dk_pcm'], 2)} pcm（{lo['burnup']:g} MWd/kg）。"
        f"絶対値が最大の Δρ {signed(far['drho_pcm'], 2)} pcm（{far['burnup']:g} MWd/kg）。"
        f"ピン出力の差の rms の最大 {rms['pin_rms_percent']:.3f}%（{rms['burnup']:g} MWd/kg）、"
        f"絶対値が最大のピン出力の差 {signed(pin['pin_max_percent'], 3)}%（{pin['burnup']:g} MWd/kg）"
    )


# 記録の表の組: (組の名前, 基準, 比べる履歴)
PAIRS = [
    (f"{lattice}-{key}", f"{lattice}-{base}", f"{lattice}-{other}")
    for lattice in LATTICES
    for key, base, other in [
        ("radial", "rings1", "rings5"),
        ("sectors4", "rings5", "rings5-sectors4"),
        ("sectors8", "rings5", "rings5-sectors8"),
        ("8-4", "rings5-sectors4", "rings5-sectors8"),
    ]
]
# 扇形の組の広がりを出す履歴と核種
SPLIT = [
    f"{lattice}-{case}" for lattice in LATTICES for case in CASES if CASES[case][1]
]
SPREAD_NUCLIDES = ["U-235", "Pu-239"]
# #0275 のピンセルの径方向の効き（本番の離散化の 5 単位 − 1 単位。pellet_zones.py の PAIRS）
PIN_CELL = (HERE.parent / "issue-0133" / "pellet-zones" / "summary.json", "5p9-1p9")


def pin_cell_rows(at: list[float]) -> dict:
    path, key = PIN_CELL
    s = json.loads(path.read_text())
    if s.get("format") != "xlc-pellet-zones":
        raise SystemExit(f"{path}: not an xlc-pellet-zones summary")
    (pair,) = [p for p in s["pairs"] if p["key"] == key]
    rows = [r for r in pair["rows"] if any(same_burnup(r["burnup"], x) for x in at)]
    if len(rows) != len(at):
        raise SystemExit(f"{path}: pair {key} lacks some burnups of {at}")
    return {
        "source": str(path.relative_to(HERE.parent.parent.parent)),
        "md5": hashlib.md5(path.read_bytes()).hexdigest(),
        "key": key,
        "rows": [
            {"burnup": r["burnup"], "dk_pcm": r["dk_pcm"], "drho_pcm": r["drho_pcm"]}
            for r in rows
        ],
    }


def compare(at: list[float], out: Path) -> None:
    names = [f"{lattice}-{case}" for lattice in LATTICES for case in CASES]
    paths = {name: HERE / f"{name}.msgpack" for name in names}
    histories = {name: summary(path) for name, path in paths.items()}
    for name, h in histories.items():
        burnups = [p["burnup"] for p in h["points"]]
        missing = [x for x in at if not any(same_burnup(y, x) for y in burnups)]
        if missing:
            raise SystemExit(f"--at {missing}: not burnup points of {name}")
    result = {
        "format": "xlc-sector-zones",
        "nuclides": SPREAD_NUCLIDES,
        "at": at,
        "inputs": {
            f"{name}.msgpack": hashlib.md5(path.read_bytes()).hexdigest()
            for name, path in paths.items()
        },
        "pairs": [
            {
                "key": key,
                "base": base,
                "other": other,
                "rows": pair_rows(histories[base], histories[other]),
            }
            for key, base, other in PAIRS
        ],
        "pin_cell": pin_cell_rows(at),
        "spreads": [
            spread_summary(
                name, sector_rows(histories[name]["history"], SPREAD_NUCLIDES), at
            )
            for name in SPLIT
        ],
    }
    out.write_text(json.dumps(result, ensure_ascii=False, indent=1) + "\n")


def report(path: Path) -> None:
    s = json.loads(path.read_text())
    if s.get("format") != "xlc-sector-zones":
        raise SystemExit(f"{path}: not an xlc-sector-zones summary")
    for pair in s["pairs"]:
        print(f"{pair['other']} − {pair['base']}: {extremes_line(pair['rows'])}")
    for spread in s["spreads"]:
        for n, m in spread["max"].items():
            print(
                f"{spread['case']} {n}: 扇形の組の広がりの最大 {m['percent']:.3f}%"
                f"（{m['burnup']:g} MWd/kg、ピン {tuple(m['position'])}、リング {m['ring']}）"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("write")
    c = sub.add_parser("compare")
    c.add_argument(
        "--at", required=True, type=lambda t: [float(x) for x in t.split(",")]
    )
    c.add_argument("out", type=Path)
    r = sub.add_parser("report")
    r.add_argument("summary", type=Path)
    args = parser.parse_args()
    match args.command:
        case "write":
            write()
        case "compare":
            compare(args.at, args.out)
        case "report":
            report(args.summary)


if __name__ == "__main__":
    main()
