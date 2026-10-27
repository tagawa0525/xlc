"""検証の記録（examples/*/verification.md）の数値の表を、比較のステージが書いた JSON から描いて置き換える（ADR-0176）。

使い方: python3 examples/verification_tables.py [--check] [<verification.md>...]

- 記録の `<!-- table: <名前> -->` から `<!-- /table -->` までを、名前の表で置き換える。目印の中は丸ごと置き換わる
  （文章は目印の外に書く）。目印の外は触らない。名前と JSON の対応（TABLES）はこのファイルにある
- `--check`: 置き換えず、描いた表と記録の中身が違えば差分を出して 1 で終わる。JSON が無い（ステージを実行していない）、git で追跡されていない（`cache: false` にして `git add` する）、dvc.lock の md5 と
  違う（JSON を手で直した、ステージを再実行したのに relock していない）、目印が壊れているときも 1 で終わる
- ファイルを渡さないと examples/*/verification.md のすべて
- 標準ライブラリだけで動く（hook と `cargo test` から呼ぶ）
"""

import difflib
import hashlib
import json
import math
import re
import subprocess
import sys
import unicodedata
from functools import partial
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OPEN = re.compile(r"^<!-- table: ([\w-]+) -->$")
CLOSE = "<!-- /table -->"


class TableError(Exception):
    pass


def signed(x: float, digits: int) -> str:
    """符号つきの小数（記録の書き方に合わせ、負の符号は U+2212）"""
    text = f"{x:+.{digits}f}"
    if float(text) == 0:  # 丸めて 0 になる負の値に符号を付けない
        text = f"{0:.{digits}f}"
    return text.replace("-", "−")


def width(text: str) -> int:
    """表の列の幅（全角は 2）。prettier と同じ数え方"""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in text)


def markdown_table(header: list[str], rows: list[list[str]]) -> str:
    """列をそろえた Markdown の表"""
    widths = [max(width(r[i]) for r in [header, *rows]) for i in range(len(header))]

    def line(cells: list[str]) -> str:
        return (
            "| "
            + " | ".join(c + " " * (w - width(c)) for c, w in zip(cells, widths))
            + " |"
        )

    separator = "| " + " | ".join("-" * w for w in widths) + " |"
    return "\n".join([line(header), separator, *map(line, rows)])


def layer_name(level: int, levels: int) -> str:
    return {1: "下端の層", levels: "上端の層"}.get(level, f"{level} 層目")


def p3_summaries(sources: list[dict]) -> list[dict]:
    for s in sources:
        if s.get("format") != "xlc-compare-p3" or "openmc" not in s["k"]:
            raise TableError(
                f"problem {s.get('problem')}: not an xlc-compare-p3 summary with OpenMC"
            )
    return sources


def p3_reference_results(sources: list[dict]) -> str:
    """Problem 3 の参照解: OpenMC と KENO-VI の k と軸方向のオフセット"""
    rows = []
    for s in p3_summaries(sources):
        k, axial = s["k"], s["axial_openmc_vs_keno"]
        rows.append(
            [
                s["problem"],
                f"{k['openmc']:.6f} ± {k['openmc_std']:.6f}",
                f"{k['keno']:.6f} ± {k['keno_std']:.6f}",
                f"{signed(k['openmc_keno_diff_pcm'], 1)} pcm",
                f"{signed(100 * axial['ao_subject'], 3)}%",
                f"{signed(100 * axial['ao_reference'], 3)}%",
            ]
        )
    return markdown_table(
        [
            "問題",
            "OpenMC（ENDF/B-VIII.0）",
            "KENO-VI（ENDF/B-VII.0）",
            "差",
            "軸方向のオフセット（OpenMC）",
            "軸方向のオフセット（KENO-VI）",
        ],
        rows,
    )


def p3_openmc_comparison(sources: list[dict]) -> str:
    """Problem 3 の炉心の計算: xlc と OpenMC の k、軸方向の出力、ピン出力、層ごとのピン出力"""
    rows = []
    for s in p3_summaries(sources):
        k, axial, pins, layer = (
            s["k"],
            s["axial_vs_openmc"],
            s["pins_vs_openmc"],
            s["layer_pins"],
        )
        rows.append(
            [
                s["problem"],
                f"{k['xlc']:.6f}",
                f"{k['openmc']:.6f}",
                f"{signed(k['diff_openmc_pcm'], 0)} ± {k['diff_openmc_std_pcm']:.0f} pcm",
                f"{axial['rms']:.4f}（{100 * axial['relative_rms']:.2f}%）",
                f"{100 * pins['rms']:.3f}%、{signed(100 * pins['max'], 3)}%",
                f"{signed(100 * layer['max'], 1)}%（{layer_name(layer['max_level'], len(layer['levels']))}）",
            ]
        )
    return markdown_table(
        [
            "問題",
            "xlc の k",
            "OpenMC の k",
            "差",
            "軸方向の出力の差の RMS",
            "ピン出力の差の RMS、最大",
            "層ごとのピン出力の差の最大",
        ],
        rows,
    )


def p5_by_problem(sources: list[dict]) -> tuple[dict, dict]:
    """compare_p5.py の要約を、KENO-VI との比較と OpenMC との比較に分けて問題で引く"""
    keno: dict[str, dict] = {}
    openmc: dict[str, dict] = {}
    for s in sources:
        if s.get("format") != "xlc-compare-p5":
            raise TableError(
                f"problem {s.get('problem')}: not an xlc-compare-p5 summary"
            )
        target = {"keno": keno, "openmc": openmc}[s["kind"]]
        if s["problem"] in target:
            raise TableError(f"{s['problem']} has two {s['kind']} summaries")
        target[s["problem"]] = s
    return keno, openmc


def p5_rod_worth(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    """制御棒の価値 (1/k_b − 1/k_a) × 1e5 pcm と、k の標準偏差を伝えた誤差（引数は (k, 標準偏差) の組）"""
    return (
        1e5 * (1 / b[0] - 1 / a[0]),
        1e5 * math.hypot(a[1] / a[0] ** 2, b[1] / b[0] ** 2),
    )


def p5_worth_text(a: tuple[float, float], b: tuple[float, float]) -> str:
    worth, std = p5_rod_worth(a, b)
    return f"{worth:.0f} ± {std:.0f} pcm"


def p5_assembly(d: dict) -> str:
    """集合体の相対差の RMS と最大（位置の名前つき）"""
    r = d["relative"]
    return f"{100 * r['rms']:.2f}%、{signed(100 * r['max'], 1)}%（{r['name']}）"


def p5_reference_results(sources: list[dict]) -> str:
    """Problem 5-2D の参照解: OpenMC と KENO-VI の k と制御棒の価値（OpenMC が無い問題の列は空欄）"""
    keno, openmc = p5_by_problem(sources)
    rows = []
    for problem, s in keno.items():
        k = s["k"]
        o = openmc.get(problem)
        base = (keno.get("5A"), openmc.get("5A"))
        row = [
            f"{problem}-2D",
            f"{o['k']['reference']:.6f} ± {o['k']['reference_std']:.6f}" if o else "",
            f"{k['reference']:.6f} ± {k['reference_std']:.6f}",
            f"{signed(o['k']['diff_keno_pcm'], 0)} ± {o['k']['diff_keno_std_pcm']:.0f} pcm"
            if o
            else "",
            "",
            "",
        ]
        if problem != "5A":
            if base[1] and o:
                row[4] = p5_worth_text(
                    (base[1]["k"]["reference"], base[1]["k"]["reference_std"]),
                    (o["k"]["reference"], o["k"]["reference_std"]),
                )
            row[5] = p5_worth_text(
                (base[0]["k"]["reference"], base[0]["k"]["reference_std"]),
                (k["reference"], k["reference_std"]),
            )
        rows.append(row)
    return markdown_table(
        [
            "問題",
            "OpenMC（ENDF/B-VIII.0）",
            "KENO-VI（ENDF/B-VII.0）",
            "差",
            "制御棒の価値（OpenMC）",
            "制御棒の価値（KENO-VI）",
        ],
        rows,
    )


def p5_openmc_comparison(sources: list[dict]) -> str:
    """Problem 5-2D の炉心の計算: xlc と OpenMC の k、集合体出力、ピン出力"""
    _, openmc = p5_by_problem(sources)
    rows = []
    for problem, s in openmc.items():
        pins = s["pins"]
        rows.append(
            [
                f"{problem}-2D",
                signed(s["k"]["diff_pcm"], 1),
                p5_assembly(s["assemblies"]),
                f"{100 * pins['all']['rms']:.2f}%、{signed(100 * pins['all']['max'], 0)}%",
                f"{100 * pins['peripheral']['rms']:.2f}%",
                f"{pins['max_pin']['xlc']:.3f}、{pins['max_pin']['reference']:.3f}",
            ]
        )
    return markdown_table(
        [
            "問題",
            "Δk/k (pcm)",
            "集合体 RMS、最大",
            "ピン RMS、最大",
            "外周の集合体のピン RMS",
            "最大のピン出力（xlc、OpenMC）",
        ],
        rows,
    )


def p5_keno_comparison(sources: list[dict]) -> str:
    """Problem 5-2D の炉心の計算: xlc と KENO-VI の k、制御棒の価値、集合体出力、最大のピン出力"""
    keno, _ = p5_by_problem(sources)
    rows = []
    for problem, s in keno.items():
        k = s["k"]
        worth = ""
        if problem != "5A":
            a = keno["5A"]["k"]
            xlc, _ = p5_rod_worth((a["xlc"], 0.0), (k["xlc"], 0.0))
            ref, _ = p5_rod_worth(
                (a["reference"], a["reference_std"]),
                (k["reference"], k["reference_std"]),
            )
            worth = f"{xlc:.0f}、{ref:.0f} pcm"
        rows.append(
            [
                f"{problem}-2D",
                f"{k['xlc']:.6f}",
                f"{k['reference']:.6f}",
                signed(k["diff_pcm"], 1),
                worth,
                p5_assembly(s["assemblies"]),
                f"{s['max_pin']['xlc']:.3f}、{s['max_pin']['reference']:.3f}",
            ]
        )
    return markdown_table(
        [
            "問題",
            "xlc の k",
            "KENO-VI の k",
            "Δk/k (pcm)",
            "制御棒の価値（xlc、KENO-VI）",
            "集合体 RMS、最大の相対差",
            "最大のピン出力（xlc、KENO-VI）",
        ],
        rows,
    )


P5_SOURCES = [
    *(f"examples/vera/p5/compare-{p}.json" for p in ("5A", "5B", "5C")),
    *(f"examples/vera/p5/compare-{p}-openmc.json" for p in ("5A", "5B", "5C")),
]


# KENO-VI との比較を、ENDF/B-VII.0 の検証用ライブラリ（examples/vera/endfb70。ADR-0175）と本番の ENDF/B-VIII.0 で並べる表。
# 描く関数に渡す JSON の前半が ENDF/B-VII.0、後半が ENDF/B-VIII.0（同じ問題の並び）
LATTICE_PROBLEMS = [
    "1B",
    "1C",
    "1D",
    "2A",
    "2B",
    "2C",
    "2D",
    "2E",
    "2F",
    "2G",
    "2H",
    "2I",
    "2J",
    "2K",
    "2Q",
]


def halves(sources: list[dict], fmt: str) -> list[tuple[dict, dict]]:
    """前半（ENDF/B-VII.0）と後半（ENDF/B-VIII.0）を、同じ問題ごとの対にする"""
    if len(sources) % 2:
        raise TableError("expected the same problems for both libraries")
    n = len(sources) // 2
    for s in sources:
        if s.get("format") != fmt:
            raise TableError(f"problem {s.get('problem')}: not a {fmt} summary")
    pairs = list(zip(sources[:n], sources[n:]))
    for vii, viii in pairs:
        if vii["problem"] != viii["problem"]:
            raise TableError(
                f"problem {vii['problem']} (VII.0) is paired with {viii['problem']} (VIII.0)"
            )
    return pairs


def pin_cell(pins: dict | None) -> str:
    return (
        ""
        if pins is None
        else f"{100 * pins['rms']:.3f}%、{signed(100 * pins['max'], 3)}%"
    )


def endfb70_lattices(sources: list[dict]) -> str:
    """ピンセルと集合体の格子計算: xlc と KENO-VI の k の差とピン出力の差を、ENDF/B-VII.0 と VIII.0 で並べる"""
    rows = []
    for vii, viii in halves(sources, "xlc-compare-keno-lattice"):
        k = vii["k"]
        rows.append(
            [
                vii["problem"],
                f"{k['xlc']:.6f}",
                f"{k['keno']:.6f} ± {k['keno_std']:.6f}",
                f"{signed(k['diff_pcm'], 0)} pcm",
                f"{signed(viii['k']['diff_pcm'], 0)} pcm",
                pin_cell(vii["pins"]),
                pin_cell(viii["pins"]),
            ]
        )
    return markdown_table(
        [
            "問題",
            "xlc の k（ENDF/B-VII.0）",
            "KENO-VI の k",
            "xlc − KENO-VI（ENDF/B-VII.0）",
            "xlc − KENO-VI（ENDF/B-VIII.0）",
            "ピン出力の差の RMS、最大（ENDF/B-VII.0）",
            "ピン出力の差の RMS、最大（ENDF/B-VIII.0）",
        ],
        rows,
    )


def endfb70_p3_k(sources: list[dict]) -> str:
    """Problem 3 の炉心の k: xlc と KENO-VI の差を ENDF/B-VII.0 と VIII.0 で、VIII.0 の OpenMC との差と並べる"""
    rows = []
    for vii, viii in halves(sources, "xlc-compare-p3"):
        k = viii["k"]
        rows.append(
            [
                vii["problem"],
                f"{vii['k']['xlc']:.6f}",
                f"{vii['k']['keno']:.6f} ± {vii['k']['keno_std']:.6f}",
                f"{signed(vii['k']['diff_keno_pcm'], 0)} pcm",
                f"{signed(k['diff_keno_pcm'], 0)} pcm",
                f"{signed(k['diff_openmc_pcm'], 0)} pcm",
                f"{signed(k['openmc_keno_diff_pcm'], 0)} pcm",
            ]
        )
    return markdown_table(
        [
            "問題",
            "xlc の k（ENDF/B-VII.0）",
            "KENO-VI の k",
            "xlc − KENO-VI（ENDF/B-VII.0）",
            "xlc − KENO-VI（ENDF/B-VIII.0）",
            "xlc − OpenMC（ENDF/B-VIII.0）",
            "OpenMC − KENO-VI（ENDF/B-VIII.0）",
        ],
        rows,
    )


def endfb70_p3_distributions(sources: list[dict]) -> str:
    """Problem 3 の炉心の軸方向の出力とピン出力: xlc と KENO-VI の差を ENDF/B-VII.0 と VIII.0 で並べる"""
    rows = []
    for vii, viii in halves(sources, "xlc-compare-p3"):
        a7, a8 = vii["axial_vs_keno"], viii["axial_vs_keno"]
        rows.append(
            [
                vii["problem"],
                f"{a7['rms']:.4f}（{100 * a7['relative_rms']:.2f}%）",
                f"{a8['rms']:.4f}（{100 * a8['relative_rms']:.2f}%）",
                f"{signed(100 * a7['ao_subject'], 3)}%",
                f"{signed(100 * a8['ao_subject'], 3)}%",
                f"{signed(100 * a7['ao_reference'], 3)}%",
                pin_cell(vii["pins_vs_keno"]),
                pin_cell(viii["pins_vs_keno"]),
            ]
        )
    return markdown_table(
        [
            "問題",
            "軸方向の出力の差の RMS（ENDF/B-VII.0）",
            "軸方向の出力の差の RMS（ENDF/B-VIII.0）",
            "軸方向のオフセット（xlc、ENDF/B-VII.0）",
            "軸方向のオフセット（xlc、ENDF/B-VIII.0）",
            "軸方向のオフセット（KENO-VI）",
            "ピン出力の差の RMS、最大（ENDF/B-VII.0）",
            "ピン出力の差の RMS、最大（ENDF/B-VIII.0）",
        ],
        rows,
    )


DEPLETION_NUCLIDES = ["U-235", "U-238", "Pu-239", "Xe-135", "Sm-149"]


def depletion_summaries(sources: list[dict]) -> list[dict]:
    for s in sources:
        if s.get("format") != "xlc-compare-depletion":
            raise TableError("not an xlc-compare-depletion summary")
    return sources


def percent_or_dash(x) -> str:
    return "–" if x is None else signed(100 * x, 2)


def burnup_label(b: float) -> str:
    return f"{b:g}"


def depletion_k(sources: list[dict]) -> str:
    """格子計算の燃焼: 燃焼点ごとの k∞ の差（xlc − OpenMC）とその 0 の点からの変化を、燃料 1 リングと 4 リングで"""
    one, four = depletion_summaries(sources)
    if one["burnup"] != four["burnup"]:
        raise TableError(
            f"the cases have other burnup points: {one['burnup']} and {four['burnup']}"
        )
    rows = []
    for b, a, c in zip(one["burnup"], one["k"], four["k"]):
        rows.append(
            [
                burnup_label(b),
                f"{signed(a['diff_pcm'], 0)} ± {a['openmc_std_pcm']:.0f}",
                signed(a["change_pcm"], 0),
                f"{signed(c['diff_pcm'], 0)} ± {c['openmc_std_pcm']:.0f}",
                signed(c["change_pcm"], 0),
            ]
        )
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "1 リング: 差 (pcm)",
            "変化 (pcm)",
            "4 リング: 差 (pcm)",
            "変化 (pcm)",
        ],
        rows,
    )


def depletion_densities(sources: list[dict]) -> str:
    """格子計算の燃焼: 燃料のピン全体の数密度の相対差（%）"""
    (s,) = depletion_summaries(sources)
    rows = [
        [burnup_label(b)]
        + [percent_or_dash(s["pin"][n][i]) for n in DEPLETION_NUCLIDES]
        for i, b in enumerate(s["burnup"])
    ]
    return markdown_table(["燃焼度 (MWd/kg)", *DEPLETION_NUCLIDES], rows)


def depletion_rings(sources: list[dict]) -> str:
    """格子計算の燃焼: 最後の燃焼点のリングごとの数密度の相対差（%）。リングは内側から"""
    (s,) = depletion_summaries(sources)
    last = len(s["rings"])

    def label(r: int) -> str:
        if last == 1:
            return "1"
        return {1: "1（内）", last: f"{last}（外）"}.get(r, str(r))

    rows = [
        [label(r + 1)] + [percent_or_dash(ring[n]) for n in DEPLETION_NUCLIDES]
        for r, ring in enumerate(s["rings"])
    ]
    return markdown_table(["リング", *DEPLETION_NUCLIDES], rows)


def depletion_times(sources: list[dict]) -> str:
    """格子計算の燃焼: 時間の内訳（秒）。xlc の全体は内訳の和、OpenMC は全体、輸送、CRAM"""
    rows = []
    for s, rings in zip(depletion_summaries(sources), ["1 リング", "4 リング"]):
        x, o = s["seconds"]["xlc"], s["seconds"]["openmc"]
        parts = [x[k] for k in ("shielding", "moc", "rates", "matrix", "cram")]
        rows.append(
            [
                f"xlc、{rings}",
                str(x["solves"]),
                f"{sum(parts):.1f}",
                f"{x['shielding']:.1f}",
                f"{x['moc']:.1f}",
                f"{x['rates']:.2f}",
                f"{x['matrix']:.2f}",
                f"{x['cram']:.2f}",
            ]
        )
        rows.append(
            [
                f"OpenMC、{rings}",
                str(o["solves"]),
                f"{o['total']:,.0f}",
                f"{o['transport']:,.0f}（1 回 {o['fastest']:.0f}〜{o['slowest']:.0f}）",
                "",
                "",
                "",
                f"{o['cram']:.2f}",
            ]
        )
    return markdown_table(
        [
            "計算",
            "輸送の回数",
            "全体",
            "自己遮蔽（OpenMC は輸送）",
            "MOC",
            "反応率と規格化",
            "燃焼行列",
            "CRAM48",
        ],
        rows,
    )


def depletion_assembly_pins(sources: list[dict]) -> str:
    """集合体の燃焼: 燃焼点ごとのピン出力の相対差（xlc / OpenMC − 1。対称操作の組で平均した値の比、%。ADR-0181）。最大は OpenMC の
    その組の統計誤差を添え、位置は上と左から 1 で数える"""
    (s,) = sources
    if s.get("format") != "xlc-compare-depletion-assembly":
        raise TableError("not an xlc-compare-depletion-assembly summary")
    if not len(s["burnup"]) == len(s["pins"]) == len(s["pin_std_rms"]):
        raise TableError(
            f"{len(s['burnup'])} burnup points with {len(s['pins'])} pin comparisons and "
            f"{len(s['pin_std_rms'])} errors; rerun the stage"
        )
    rows = []
    for b, p, std in zip(s["burnup"], s["pins"], s["pin_std_rms"]):
        row, column = p["at"]
        rows.append(
            [
                burnup_label(b),
                f"{100 * p['rms']:.3f}",
                f"{signed(100 * p['max'], 3)} ± {100 * p['max_std']:.3f}",
                f"({row + 1}, {column + 1})",
                f"{100 * std:.3f}",
            ]
        )
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "rms (%)",
            "最大 (%)",
            "最大の位置（行, 列。1 から）",
            "OpenMC の統計誤差 rms (%)",
        ],
        rows,
    )


def assembly_summary(sources: list[dict]) -> dict:
    (s,) = sources
    return assembly_checked(s)


def assembly_checked(s: dict) -> dict:
    if s.get("format") != "xlc-compare-depletion-assembly":
        raise TableError("not an xlc-compare-depletion-assembly summary")
    return s


def depletion_assembly_k(sources: list[dict]) -> str:
    """集合体の燃焼: 燃焼点ごとの k∞ の差（xlc − OpenMC。± は OpenMC の標準偏差）と 0 の点からの変化"""
    s = assembly_summary(sources)
    rows = [
        [
            burnup_label(b),
            f"{signed(k['diff_pcm'], 0)} ± {k['openmc_std_pcm']:.0f}",
            signed(k["change_pcm"], 0),
        ]
        for b, k in zip(s["burnup"], s["k"], strict=True)
    ]
    return markdown_table(["燃焼度 (MWd/kg)", "差 (pcm)", "変化 (pcm)"], rows)


def depletion_assembly_free_gas_k(sources: list[dict]) -> str:
    """集合体の燃焼: 静止した標的の本番と自由ガスの試行の、燃焼点ごとの k∞ の差（xlc − OpenMC。± は OpenMC の標準偏差）と 0 の点からの変化"""
    base, trial = (assembly_checked(s) for s in sources)
    if base["problem"] != trial["problem"]:
        raise TableError(
            f"problems {base['problem']} and {trial['problem']} are not the same"
        )
    if base["burnup"] != trial["burnup"]:
        raise TableError("the two summaries have different burnup points")
    rows = []
    for b, kb, kt in zip(base["burnup"], base["k"], trial["k"], strict=True):
        if kb["openmc_std_pcm"] != kt["openmc_std_pcm"]:
            raise TableError(f"the OpenMC errors differ at {b} MWd/kg")
        std = f"{kb['openmc_std_pcm']:.0f}"
        rows.append(
            [
                burnup_label(b),
                f"{signed(kb['diff_pcm'], 0)} ± {std}",
                signed(kb["change_pcm"], 0),
                f"{signed(kt['diff_pcm'], 0)} ± {std}",
                signed(kt["change_pcm"], 0),
            ]
        )
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "静止した標的 差 (pcm)",
            "静止した標的 変化 (pcm)",
            "自由ガス 差 (pcm)",
            "自由ガス 変化 (pcm)",
        ],
        rows,
    )


def assembly_density_table(s: dict, columns: list[tuple[str, list]]) -> str:
    rows = [
        [burnup_label(b)] + [percent_or_dash(series[i]) for _, series in columns]
        for i, b in enumerate(s["burnup"])
    ]
    return markdown_table(["燃焼度 (MWd/kg)", *[name for name, _ in columns]], rows)


def depletion_assembly_fuel_densities(sources: list[dict]) -> str:
    """集合体の燃焼: 燃料の数密度の相対差（%。燃料の領域の体積の重みの平均。xlc / OpenMC − 1）"""
    s = assembly_summary(sources)
    return assembly_density_table(s, list(s["kinds"]["fuel"].items()))


def depletion_assembly_absorber_densities(sources: list[dict]) -> str:
    """集合体の燃焼: 燃料以外の物質（被覆、Pyrex、WABA の B-10、Gd 燃料のリングの Gd）の数密度の相対差（%）。列は「物質 核種」"""
    s = assembly_summary(sources)
    columns = [
        (f"{kind} {name}", series)
        for kind, table in s["kinds"].items()
        if kind != "fuel"
        for name, series in table.items()
    ]
    if not columns:
        raise TableError(f"problem {s['problem']} has no material except the fuel")
    return assembly_density_table(s, columns)


def depletion_assembly_times(sources: list[dict]) -> str:
    """集合体の燃焼: xlc の時間の内訳（秒）。全体は内訳の和。OpenMC の時間は AWS の機種によるので載せない"""
    rows = []
    for s in sources:
        x = assembly_checked(s)["seconds"]["xlc"]
        parts = [x[k] for k in ("shielding", "moc", "rates", "matrix", "cram")]
        rows.append(
            [
                s["problem"],
                str(x["solves"]),
                f"{sum(parts):.1f}",
                f"{x['shielding']:.1f}",
                f"{x['moc']:.1f}",
                f"{x['rates']:.2f}",
                f"{x['matrix']:.2f}",
                f"{x['cram']:.2f}",
            ]
        )
    return markdown_table(
        [
            "問題",
            "輸送の回数",
            "全体",
            "自己遮蔽",
            "MOC",
            "反応率と規格化",
            "燃焼行列",
            "CRAM48",
        ],
        rows,
    )


# 集合体の燃焼の問題（examples/depletion/p2/dvc.yaml の jobs）
ASSEMBLY_DEPLETION_JOBS = ["2A", "2F", "2M", "2N", "2O", "2P", "2A-density", "2A-rods"]


def aws_summaries(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-compare-aws-depletion":
        raise TableError("not an xlc-compare-aws-depletion summary")
    return s


def aws_depletion_k(sources: list[dict]) -> str:
    """区間に分けた AWS の燃焼と手元の燃焼の、燃焼点ごとの k と差（AWS − 手元）。倍率は合成した σ に対するもの"""
    s = aws_summaries(sources)
    rows = [
        [
            burnup_label(b),
            f"{k['local']:.6f} ± {k['local_std_pcm'] / 1e5:.6f}",
            f"{k['aws']:.6f} ± {k['aws_std_pcm'] / 1e5:.6f}",
            signed(k["diff_pcm"], 1),
            f"{k['sigma_pcm']:.1f}",
            signed(k["ratio"], 2),
        ]
        for b, k in zip(s["burnup"], s["k"], strict=True)
    ]
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "k 手元 ± σ",
            "k AWS ± σ",
            "差 (pcm)",
            "合成 σ (pcm)",
            "倍率",
        ],
        rows,
    )


def aws_depletion_densities(sources: list[dict]) -> str:
    """区間に分けた AWS の燃焼と手元の燃焼の、燃料の数密度の相対差（%。AWS − 手元）"""
    s = aws_summaries(sources)
    rows = [
        [burnup_label(b)]
        + [
            "–" if s["density"][n][i] is None else signed(100 * s["density"][n][i], 3)
            for n in DEPLETION_NUCLIDES
        ]
        for i, b in enumerate(s["burnup"])
    ]
    return markdown_table(["燃焼度 (MWd/kg)", *DEPLETION_NUCLIDES], rows)


def core_depletion_k(sources: list[dict]) -> str:
    """core の燃焼（集合体 1 体）: 燃焼点ごとの格子計算の履歴と core の k∞ と反応度の差（core − 履歴。ADR-0190）と、テーブルの点の
    分岐の格子計算との反応度の差（core − 分岐）"""
    (s,) = sources
    if s.get("format") != "xlc-compare-core-depletion":
        raise TableError("not an xlc-compare-core-depletion summary")
    rows = [
        [
            burnup_label(p["burnup"]),
            f"{p['efpd']:.2f}",
            "○" if p["on_table"] else "",
            f"{p['history_k']:.8f}",
            f"{p['k']:.8f}",
            signed(p["rho_pcm"], 3),
            "" if p["branch_rho_pcm"] is None else signed(p["branch_rho_pcm"], 4),
        ]
        for p in s["points"]
    ]
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "EFPD",
            "テーブルの点",
            "履歴の k∞",
            "core の k∞",
            "core − 履歴 (pcm)",
            "core − 分岐 (pcm)",
        ],
        rows,
    )


def burnup_interpolation_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-compare-burnup-interpolation":
        raise TableError("not an xlc-compare-burnup-interpolation summary")
    return s


def burnup_interpolation_methods(sources: list[dict]) -> str:
    """燃焼度の内挿の方式（ADR-0192）: 方式ごと、テーブルと点の選び方ごとの、間引いた点の |Δρ| の最大と rms（pcm）と、
    真値がデータの範囲にある項目のオーバーシュートの最大（相対）"""
    s = burnup_interpolation_summary(sources)
    selections = [(t, sel) for t in s["tables"] for sel in t["selections"]]
    rows = []
    for method, label in s["methods"].items():
        results = [sel["results"][method] for _, sel in selections]
        overshoot = max((r["overshoot"] for r in results), default=0.0)
        rows.append(
            [
                label,
                *(f"{r['max_abs_pcm']:.1f} / {r['rms_pcm']:.1f}" for r in results),
                "0" if overshoot == 0.0 else f"{overshoot:.1e}",
            ]
        )
    return markdown_table(
        [
            "方式",
            *(
                f"{t['name']} {sel['name']}（{len(sel['burnup'])} 点）"
                for t, sel in selections
            ),
            "オーバーシュートの最大",
        ],
        rows,
    )


def burnup_interpolation_points(sources: list[dict]) -> str:
    """燃焼度の点の選び方（ADR-0192）: 採った方式（3 点の放物線の傾きの Hermite）の、選び方ごとの間引いた点の誤差（pcm）"""
    s = burnup_interpolation_summary(sources)
    rows = []
    for t in s["tables"]:
        for sel in t["selections"]:
            r = sel["results"]["hermite"]
            rho = r["rho_pcm"]
            worst = max(range(len(rho)), key=lambda i: abs(rho[i]), default=None)
            rows.append(
                [
                    t["name"],
                    sel["name"],
                    str(len(sel["burnup"])),
                    str(len(sel["removed"])),
                    ""
                    if worst is None
                    else f"{signed(rho[worst], 1)}（{burnup_label(sel['removed'][worst])}）",
                    f"{r['rms_pcm']:.1f}",
                    f"{signed(min(rho), 1)}〜{signed(max(rho), 1)}" if rho else "",
                ]
            )
    return markdown_table(
        [
            "集合体",
            "点の選び方",
            "点の数",
            "比べた点の数",
            "最大の Δρ (pcm)（燃焼度）",
            "rms (pcm)",
            "範囲 (pcm)",
        ],
        rows,
    )


def signed_exponent(x: float) -> str:
    """符号つきの指数表記の 2 桁（`signed` と同じく、負の符号は U+2212、0 は符号なし。指数の符号は ASCII のまま）"""
    if x == 0:
        return "0.0e+00"
    mantissa, exponent = f"{x:+.1e}".split("e")
    return f"{mantissa.replace('-', '−')}e{exponent}"


def core_depletion_boron(sources: list[dict]) -> str:
    """core の燃焼（集合体 1 体）で点ごとに探した臨界ホウ素（ADR-0191）と、同じ燃焼度の一様な状態の探索（ADR-0145）"""
    (s,) = sources
    if s.get("format") != "xlc-compare-core-boron":
        raise TableError("not an xlc-compare-core-boron summary")
    rows = [
        [
            burnup_label(p["burnup"]),
            f"{p['efpd']:.2f}",
            f"{p['boron']:.4f}",
            str(p["solves"]),
            f"{p['k']:.8f}",
            signed_exponent(p["k_minus_target"]),
            f"{p['uniform_boron']:.4f}",
            str(p["uniform_solves"]),
            signed(p["boron_rounded_difference"], 4),
        ]
        for p in s["points"]
    ]
    return markdown_table(
        [
            "燃焼度 (MWd/kg)",
            "EFPD",
            "ホウ素 (ppm)",
            "回数",
            "k",
            "k − 目標",
            "一様な状態のホウ素 (ppm)",
            "一様な状態の回数",
            "差 (ppm)",
        ],
        rows,
    )


def branch_cross_terms_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-compare-branch-cross-terms":
        raise TableError("not an xlc-compare-branch-cross-terms summary")
    return s


def branch_cross_terms(sources: list[dict]) -> str:
    """分岐のテーブルの交差の項（ADR-0200）: モデルごと、テーブルごとの、組み立てた状態の |Δρ| の最大（pcm。燃焼度の点が
    複数あればその最大）と、モデルが解く状態の数"""
    s = branch_cross_terms_summary(sources)
    models = []
    for t in s["tables"]:
        for p in t["points"]:
            models += [m for m in p["models"] if m not in models]
    rows = []
    for m in models:
        cells = []
        for t in s["tables"]:
            results = [p["models"][m] for p in t["points"] if m in p["models"]]
            if results:
                worst = max(r["max_abs_pcm"] for r in results)
                cells.append(f"{worst:.2f}（{results[0]['states']}）")
            else:
                cells.append("")
        rows.append([m, *cells])
    return markdown_table(["モデル", *(t["name"] for t in s["tables"])], rows)


def branch_thinning(sources: list[dict]) -> str:
    """分岐を燃焼度の一部の点だけで解いたとき（ADR-0200）: 選び方ごとの、間引いた点の |Δρ| の最大（pcm）。Hermite は最大の燃焼度を添える"""
    s = branch_cross_terms_summary(sources)
    rows = []
    for t in s["thinning"]:
        for sel in t["selections"]:
            linear = sel["results"]["linear"]
            hermite = sel["results"]["hermite"]
            rows.append(
                [
                    t["name"],
                    sel["name"],
                    str(len(sel["burnup"])),
                    f"{linear['max_abs_pcm']:.2f}",
                    f"{hermite['max_abs_pcm']:.2f}（{hermite['worst']['burnup']:g}）",
                ]
            )
    return markdown_table(
        ["テーブル", "選び方", "点", "線形 (pcm)", "Hermite (pcm、燃焼度)"], rows
    )


def branch_shielding_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-compare-branch-shielding":
        raise TableError("not an xlc-compare-branch-shielding summary")
    return s


def branch_shielding_groups(s: dict) -> dict:
    """集合体（ケースの名前の `-` の前）と燃焼度ごとの状態（ケースの順、状態の順）"""
    groups: dict = {}
    for case in s["cases"]:
        assembly = case["name"].split("-")[0]
        for state in case["states"]:
            groups.setdefault((assembly, state["burnup"]), []).append(state)
    return groups


def branch_shielding_k(sources: list[dict]) -> str:
    """分岐の自己遮蔽を公称の減速材で解いたテーブルの、状態ごとに解き直したテーブルとの k∞ の差（pcm。ADR-0201）。行はホウ素、
    列は減速材の温度と密度"""
    s = branch_shielding_summary(sources)
    groups = branch_shielding_groups(s)
    moderators = sorted(
        {
            (state["moderator_temperature"], state["density"])
            for states in groups.values()
            for state in states
        }
    )
    rows = []
    for (assembly, burnup), states in groups.items():
        for boron in sorted({state["boron"] for state in states}):
            at: dict = {}
            for state in states:
                if state["boron"] != boron:
                    continue
                key = (state["moderator_temperature"], state["density"])
                # 燃料の温度の軸は列にないので、1 点のテーブルだけを描く（同じ減速材の点の状態を黙って上書きしない）
                if key in at:
                    raise TableError(
                        f"{assembly} {burnup} MWd/kg {boron} ppm: two states at moderator {key}"
                    )
                at[key] = state["dk_pcm"]
            rows.append(
                [
                    assembly,
                    burnup_label(burnup),
                    f"{boron:g}",
                    *(signed(at[m], 2) if m in at else "–" for m in moderators),
                ]
            )
    return markdown_table(
        [
            "集合体",
            "燃焼度 (MWd/kg)",
            "ホウ素 (ppm)",
            *(f"{t:g} K、{d:g} g/cc" for t, d in moderators),
        ],
        rows,
    )


def branch_shielding_constants(sources: list[dict]) -> str:
    """分岐の自己遮蔽を公称の減速材で解いたテーブルの、燃焼度の点ごとの状態の Δk∞ の絶対値の最大（pcm）と、少数群定数の項目ごとの
    相対差の最大（ADR-0201）"""
    s = branch_shielding_summary(sources)
    rows = []
    for (assembly, burnup), states in branch_shielding_groups(s).items():
        rows.append(
            [
                assembly,
                burnup_label(burnup),
                str(len(states)),
                f"{max(abs(state['dk_pcm']) for state in states):.2f}",
                *(
                    f"{max(state['relative'][item] for state in states):.1e}"
                    for item in s["items"]
                ),
            ]
        )
    labels = {
        "absorption": "Σa",
        "nu_fission": "νΣf",
        "diffusion": "D",
        "scatter": "Σs",
        "adf": "ADF",
    }
    return markdown_table(
        [
            "集合体",
            "燃焼度 (MWd/kg)",
            "状態の数",
            "Δk∞ の絶対値の最大 (pcm)",
            *(labels[item] for item in s["items"]),
        ],
        rows,
    )


def swap_candidates(sources: list[dict]) -> str:
    """#0133 の候補ごとの、断面積を OpenMC の値に置き換えて測る k への寄与（pcm）。状態ごとの値と、最初の状態からの変化"""
    (s,) = sources
    if s.get("format") != "xlc-swap-candidates":
        raise TableError("not an xlc-swap-candidates summary")
    states = s["states"]
    first, later = states[0], states[1:]

    def value(x) -> str:
        return "–" if x is None else signed(x, 1)

    def change(a, b) -> str:
        return "–" if a is None or b is None else signed(b - a, 1)

    rows = [
        [
            "k の差（xlc − OpenMC）",
            *(f"{signed(t['diff_pcm'], 0)} ± {t['diff_std_pcm']:.0f}" for t in states),
            *(signed(t["diff_pcm"] - first["diff_pcm"], 0) for t in later),
        ]
    ]
    for key in s["order"]:
        pcm = [t["items"].get(key, {}).get("pcm") for t in states]
        # 組成に核種が無くて置き換えの無い状態（新しい燃料の Pu など）の寄与は 0 なので、変化は 0 からのもの
        start = 0.0 if key not in first["items"] else pcm[0]
        rows.append(
            [
                s["labels"][key],
                *(value(x) for x in pcm),
                *(change(start, x) for x in pcm[1:]),
            ]
        )
    return markdown_table(
        [
            "項目",
            *(t["name"] for t in states),
            *(f"{t['name']} − {first['name']}" for t in later),
        ],
        rows,
    )


def ufg_free_gas_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-ufg-free-gas":
        raise TableError("not an xlc-ufg-free-gas summary")
    return s


def ufg_free_gas_k(sources: list[dict]) -> str:
    """ADR-0203 の UFG の散乱核ごとの k の差（xlc − OpenMC、pcm）と、最初の状態からの増え"""
    s = ufg_free_gas_summary(sources)
    states = s["states"]
    first, later = states[0], states[1:]
    rows = []
    for mode in s["modes"]:
        key = mode["key"]
        diffs = [t["results"].get(key, {}).get("diff_pcm") for t in states]
        rows.append(
            [
                mode["label"],
                *("–" if d is None else signed(d, 1) for d in diffs),
                *(
                    "–" if d is None or diffs[0] is None else signed(d - diffs[0], 1)
                    for d in diffs[1:]
                ),
            ]
        )
    return markdown_table(
        [
            "散乱核",
            *(t["name"] for t in states),
            *(f"{t['name']} − {first['name']}" for t in later),
        ],
        rows,
    )


def ufg_free_gas_swaps(
    sources: list[dict],
    specs: list[str] | None = None,
    modes: list[str] | None = None,
) -> str:
    """ADR-0203 の UFG の散乱核ごとの、UFG の範囲の吸収を OpenMC の値に置き換えて測る k への寄与（pcm）。`specs` と `modes` は
    残す置き換えと散乱核（省略するとすべて。要約に無いものは止める）。状態が流していない置き換え（状態の `swaps` に無いもの。燃料に無い
    核種）は –"""
    s = ufg_free_gas_summary(sources)
    for wanted, known in [
        (specs, [x["spec"] for x in s["swaps"]]),
        (modes, [m["key"] for m in s["modes"]]),
    ]:
        unknown = sorted(set(wanted or []) - set(known))
        if unknown:
            raise TableError(f"the summary has no {unknown}")
    swaps = [x for x in s["swaps"] if specs is None or x["spec"] in specs]
    kept = [m for m in s["modes"] if modes is None or m["key"] in modes]

    def cell(state: dict, values: dict, spec: str) -> str:
        # 状態が流していない置き換え（燃料に無い核種）だけを –、流したのに結果が無いものは止める
        if spec not in state["swaps"]:
            return "–"
        if spec not in values:
            raise TableError(f"{state['name']}: no result of the swap {spec}")
        return signed(values[spec], 1)

    rows = [
        [
            t["name"],
            mode["label"],
            *(cell(t, t["results"][mode["key"]]["swaps"], x["spec"]) for x in swaps),
        ]
        for t in s["states"]
        for mode in kept
        if mode["key"] in t["results"]
    ]
    return markdown_table(["状態", "散乱核", *(x["label"] for x in swaps)], rows)


# ADR-0204 の表の置き換えの列: UFG の範囲を 2 つに分けたもの（ADR-0203 と同じ）と、U-238 と Pu-239 の 4〜48 eV の群ごと
UFG_RANGE_SWAPS = ["U-238:27-27", "*:27-27", "U-238:14-26", "*:14-26"]
UFG_GROUP_SWAPS = [
    f"{n}:{g}-{g}" for n in ["U-238", "Pu-239"] for g in [27, 26, 25, 24]
]


def ufg_free_gas_tilt(
    sources: list[dict], energies: tuple[float, float] | None = None
) -> str:
    """ADR-0203 の UFG の散乱核ごとの、燃料の中性子束の区間ごとの比 − 1（%。群の中で規格化。fine_compare.py）。`energies` は残す
    区間（中点）のエネルギーの範囲 (eV)（省略するとすべての区間）"""
    s = ufg_free_gas_summary(sources)
    rows, bins = [], None
    for t in s["states"]:
        for mode in s["modes"]:
            tilt = t["results"].get(mode["key"], {}).get("tilt")
            if tilt is None:
                continue
            if energies is not None:
                low, high = energies
                # 区間の端は fine_compare.py が 4 桁に丸めて書くので、群の境界と比べずに区間の中点で選ぶ
                tilt = [b for b in tilt if low < 0.5 * (b["low"] + b["high"]) < high]
                if not tilt:
                    raise TableError(
                        f"{t['name']}, {mode['label']}: the fine log has no bins in {low}-{high} eV"
                    )
            edges = [(b["low"], b["high"]) for b in tilt]
            if bins is None:
                bins = edges
            elif edges != bins:
                raise TableError("the fine logs use other bins")
            rows.append(
                [t["name"], mode["label"], *(signed(100 * b["ratio"], 2) for b in tilt)]
            )
    if not bins:
        raise TableError("no fine logs")
    return markdown_table(
        ["状態", "散乱核", *(f"{lo:g}〜{hi:g} eV" for lo, hi in bins)], rows
    )


def ufg_free_gas_transfers(sources: list[dict]) -> str:
    """ADR-0205 の UFG の散乱核ごとの、物質の移行の行を OpenMC の値に置き換えて測る k への寄与（pcm）。移行の置き換えを流した
    散乱核だけを並べ、要約の列（`transfers`）の置き換えが無ければ止める"""
    s = ufg_free_gas_summary(sources)
    columns = s.get("transfers", [])
    if not columns:
        raise TableError("the summary lists no transfer swaps")
    rows = []
    for t in s["states"]:
        for mode in s["modes"]:
            values = t["results"].get(mode["key"], {}).get("transfers")
            if values is None:
                continue
            missing = [c["spec"] for c in columns if c["spec"] not in values]
            if missing:
                raise TableError(
                    f"{t['name']}, {mode['label']}: no transfer swaps {missing}"
                )
            rows.append(
                [
                    t["name"],
                    mode["label"],
                    *(signed(values[c["spec"]], 1) for c in columns),
                ]
            )
    if not rows:
        raise TableError("no transfer logs")
    return markdown_table(["状態", "散乱核", *(c["label"] for c in columns)], rows)


def group_rows_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-group-rows":
        raise TableError("not an xlc-group-rows summary")
    return s


def group_rows_swaps(sources: list[dict]) -> str:
    """#0152 の、物質の移行の行を、状態の OpenMC の行と、核の実行の細かい入射の核を xlc の群の中の中性子束か 1/E でまとめた行に
    置き換えて測る k への寄与（pcm）。群の中の重みの寄与は 2 つの重みの差、中性子束の差は OpenMC の行と xlc の重みの行の差"""
    s = group_rows_summary(sources)
    labels = {k["key"]: k["label"] for k in s["kernels"]}
    rows = []
    for t in s["states"]:
        for spec in s["specs"]:
            swap = t["swaps"].get(spec["spec"])
            if swap is None:
                continue
            for kernel in t["kernels"]:
                k = swap["kernels"][kernel]
                rows.append(
                    [
                        t["name"],
                        spec["label"],
                        labels[kernel],
                        signed(swap["openmc"], 1),
                        signed(k["xlc"], 1),
                        signed(k["lethargy"], 1),
                        signed(k["xlc"] - k["lethargy"], 1),
                        signed(swap["openmc"] - k["xlc"], 1),
                    ]
                )
    return markdown_table(
        [
            "状態",
            "置き換え",
            "核",
            "OpenMC の行",
            "核 × xlc の中性子束",
            "核 × 1/E",
            "群の中の重み",
            "中性子束の差",
        ],
        rows,
    )


def group_rows_out(sources: list[dict], key: str) -> str:
    """#0152 の、物質と群 `key`（`<物質>:<群>`）の群から出る移行（1/cm）。状態の OpenMC の値と、xlc、核の実行の細かい入射の核を
    核の実行、xlc、1/E の重みでまとめた値の、OpenMC との比 − 1"""
    s = group_rows_summary(sources)
    labels = {k["key"]: k["label"] for k in s["kernels"]}
    rows = []
    for t in s["states"]:
        for kernel in t["kernels"]:
            row = t["rows"][kernel].get(key)
            if row is None:
                continue
            reference = row["state"][0]

            def ratio(name: str, row=row, reference=reference) -> str:
                return f"{signed((row[name][0] / reference - 1) * 100, 2)}%"

            rows.append(
                [
                    t["name"],
                    labels[kernel],
                    f"{reference:.5f}",
                    ratio("xlc"),
                    ratio("kernel-openmc"),
                    ratio("kernel-xlc"),
                    ratio("kernel-lethargy"),
                ]
            )
    if not rows:
        raise TableError(f"no rows of {key}")
    return markdown_table(
        [
            "状態",
            "核",
            "OpenMC (1/cm)",
            "xlc",
            "核 × 核の実行の中性子束",
            "核 × xlc の中性子束",
            "核 × 1/E",
        ],
        rows,
    )


def group_rows_spectrum(sources: list[dict]) -> str:
    """#0152 の、群の中の区間ごとの中性子束の xlc / OpenMC − 1（%。核の実行のある状態だけ）"""
    s = group_rows_summary(sources)
    rows = []
    header = None
    for t in s["states"]:
        if t["spectrum"] is None:
            continue
        for material, bins in t["spectrum"].items():
            names = [f"{b['low']:g}〜{b['high']:g} eV" for b in bins]
            if header is None:
                header = names
            elif names != header:
                raise TableError(f"{t['name']}, {material}: other bins {names}")
            rows.append(
                [t["name"], material, *(signed(b["ratio"] * 100, 2) for b in bins)]
            )
    if header is None:
        raise TableError("no state has a fine run")
    return markdown_table(["状態", "物質", *header], rows)


# dvc.yaml のステージの出力でなく、git が追跡する JSON（dvc.lock の md5 の照合をしない）。AWS の結果は references/ にあり、DVC のステージにしない（ADR-0174）。
# 燃焼度の軸を持つテーブルの確認は、履歴のステージを deps にすると重いステージが走るのでステージにしない（examples/depletion/verification.md）
GIT_TRACKED = {
    "examples/depletion/aws-check/compare-aws-1C-rings1.json",
    "examples/depletion/branches/compare-core-2A.json",
    "examples/depletion/branches/compare-burnup-interpolation.json",
    "examples/depletion/branches/compare-core-2A-boron.json",
    "examples/depletion/branches/compare-branch-cross-terms.json",
    "examples/depletion/branches/compare-shielding.json",
    "examples/depletion/branches/compare-shielding-boron.json",
    "examples/depletion/issue-0133/swap-candidates.json",
    "examples/depletion/issue-0133/ufg-free-gas/summary.json",
    "examples/depletion/issue-0133/ufg-heavy-free-gas/summary.json",
    "examples/depletion/issue-0133/ufg-sab/summary.json",
    "examples/depletion/issue-0133/group-rows/summary.json",
    "examples/depletion/issue-0133/ufg-transfer/summary.json",
    "examples/depletion/issue-0133/ufg-transfer/summary-thin.json",
    "examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-2A-fg.json",
    "examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-2O-fg.json",
    "examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-2A-static.json",
    "examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-2O-static.json",
}


def ufg_transfer_summary(sources: list[dict]) -> dict:
    (s,) = sources
    if s.get("format") != "xlc-ufg-transfer":
        raise TableError("not an xlc-ufg-transfer summary")
    return s


def ufg_transfer_k(sources: list[dict]) -> str:
    """ADR-0207 の、UFG の範囲の群の移行の形ごとの k の差（xlc − OpenMC、pcm）と、最初の状態からの増え"""
    s = ufg_transfer_summary(sources)
    states = s["states"]
    first = states[0]
    header = (
        ["移行の形"]
        + [t["name"] for t in states]
        + [f"{t['name']} − {first['name']}" for t in states[1:]]
    )
    rows = [
        ["OpenMC の k"]
        + [f"{t['openmc_k']:.5f}" for t in states]
        + ["–" for _ in states[1:]]
    ]
    for mode in s["modes"]:
        diffs = [(t["modes"][mode["key"]]["k"] - t["openmc_k"]) * 1e5 for t in states]
        rows.append(
            [mode["label"]]
            + [signed(d, 1) for d in diffs]
            + [signed(d - diffs[0], 1) for d in diffs[1:]]
        )
    return markdown_table(header, rows)


def ufg_transfer_swaps(sources: list[dict]) -> str:
    """ADR-0207 の、物質の移行の行を状態の OpenMC の行に置き換えた k への寄与（pcm）を、移行の形ごとに。最後の行は、最初の状態からの
    寄与の増え"""
    s = ufg_transfer_summary(sources)
    modes = s["modes"]
    states = s["states"]
    first = states[0]

    def values(t: dict, spec: str) -> list[float]:
        return [t["modes"][m["key"]]["swaps"][spec] for m in modes]

    rows = []
    for t in states:
        for spec in s["specs"]:
            v = values(t, spec["spec"])
            rows.append(
                [
                    t["name"],
                    spec["label"],
                    *[signed(x, 1) for x in v],
                    signed(v[-1] - v[0], 1),
                ]
            )
    for t in states[1:]:
        for spec in s["specs"]:
            v = [
                a - b
                for a, b in zip(values(t, spec["spec"]), values(first, spec["spec"]))
            ]
            rows.append(
                [
                    f"{t['name']} − {first['name']}",
                    spec["label"],
                    *[signed(x, 1) for x in v],
                    signed(v[-1] - v[0], 1),
                ]
            )
    return markdown_table(
        [
            "状態",
            "置き換え",
            *[m["label"] for m in modes],
            f"{modes[-1]['label']} − {modes[0]['label']}",
        ],
        rows,
    )


def ufg_transfer_rows(sources: list[dict]) -> str:
    """ADR-0207 の、群から出る移行の、状態の OpenMC の行との比 − 1（%）。xlc は移行の形ごと、「核 × xlc の中性子束」は OpenMC の細かい
    入射の核を xlc のリングの中性子束でまとめた行（#0152）。14〜26 群は 14 群と、15〜26 群のうち絶対値の最も大きいもの"""
    s = ufg_transfer_summary(sources)
    modes = [m["key"] for m in s["modes"]]
    labels = {k["key"]: k["label"] for k in s["kernels"]}
    rows = []
    for t in s["states"]:
        for kernel in t["kernels"]:
            found = {mode: t["rows"][mode][kernel] for mode in modes}
            materials = sorted(
                {name.split(":")[0] for name in found[modes[0]]},
                key=["water", "fuel"].index,
            )
            for material in materials:
                spans = [("27", [27]), ("14", [14]), ("15〜26", list(range(15, 27)))]
                for label, numbers in spans:
                    names = [
                        f"{material}:{g}"
                        for g in numbers
                        if f"{material}:{g}" in found[modes[0]]
                    ]
                    if not names:
                        continue

                    def ratio(
                        mode: str,
                        column: str,
                        of: str = "state",
                        found=found,
                        names=names,
                    ) -> float:
                        values = [
                            found[mode][n][column][0] / found[mode][n][of][0] - 1
                            for n in names
                        ]
                        return max(values, key=abs) * 100

                    rows.append(
                        [
                            t["name"],
                            labels[kernel],
                            material,
                            label,
                            f"{found[modes[0]][names[0]]['state'][0]:.5f}"
                            if len(numbers) == 1
                            else "–",
                            *[f"{signed(ratio(mode, 'xlc'), 2)}%" for mode in modes],
                            f"{signed(ratio(modes[-1], 'kernel-xlc'), 2)}%",
                            f"{signed(ratio(modes[-1], 'xlc', 'kernel-xlc'), 2)}%",
                        ]
                    )
    if not rows:
        raise TableError("no rows")
    return markdown_table(
        [
            "状態",
            "核",
            "物質",
            "群",
            "OpenMC (1/cm)",
            *[m["label"] for m in s["modes"]],
            "核 × xlc の中性子束",
            f"{s['modes'][-1]['label']} / 核 × xlc の中性子束",
        ],
        rows,
    )


def counted_pin_cell(pins: dict | None) -> str:
    """比べた点の数を添える。OpenMC − KENO-VI は 1/8 の表、xlc − OpenMC は燃料ピンすべてなので、点の数が違う"""
    return "" if pins is None else f"{pin_cell(pins)}（{pins['positions']} 点）"


def endfb70_openmc(sources: list[dict]) -> str:
    """ENDF/B-VII.0 の OpenMC（NJOY）を基準に、xlc と KENO-VI の k の差とピン出力の差を分ける（#0131）"""
    rows = []
    for s in sources:
        if s.get("format") != "xlc-compare-openmc-keno":
            raise TableError(
                f"problem {s.get('problem')}: not an xlc-compare-openmc-keno summary"
            )
        k = s["k"]
        rows.append(
            [
                s["problem"],
                f"{k['xlc']:.6f}",
                f"{k['openmc']:.6f} ± {k['openmc_std']:.6f}",
                f"{k['keno']:.6f} ± {k['keno_std']:.6f}",
                f"{signed(k['diff_xlc_openmc_pcm'], 0)} ± {k['diff_xlc_openmc_std_pcm']:.0f} pcm",
                f"{signed(k['diff_openmc_keno_pcm'], 0)} ± {k['diff_openmc_keno_std_pcm']:.0f} pcm",
                f"{signed(k['diff_xlc_keno_pcm'], 0)} pcm",
                counted_pin_cell(s["pins_openmc_vs_keno"]),
                counted_pin_cell(s["pins_xlc_vs_openmc"]),
            ]
        )
    return markdown_table(
        [
            "問題",
            "xlc の k",
            "OpenMC の k",
            "KENO-VI の k",
            "xlc − OpenMC",
            "OpenMC − KENO-VI",
            "xlc − KENO-VI",
            "ピン出力の差の RMS、最大（OpenMC − KENO-VI）",
            "ピン出力の差の RMS、最大（xlc − OpenMC）",
        ],
        rows,
    )


# 表の名前 → (JSON（リポジトリの直下から）, 描く関数)。JSON は dvc.yaml のステージの outs で、同じ場所の dvc.lock が md5 を持つ
TABLES = {
    "p3-reference-results": (
        [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        p3_reference_results,
    ),
    "p5-reference-results": (
        P5_SOURCES,
        p5_reference_results,
    ),
    "p5-openmc-comparison": (
        [f"examples/vera/p5/compare-{p}-openmc.json" for p in ("5A", "5B", "5C")],
        p5_openmc_comparison,
    ),
    "p5-keno-comparison": (
        [f"examples/vera/p5/compare-{p}.json" for p in ("5A", "5B", "5C")],
        p5_keno_comparison,
    ),
    "p3-openmc-comparison": (
        [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        p3_openmc_comparison,
    ),
    "endfb70-lattices": (
        [f"examples/vera/endfb70/compare-{p}.json" for p in LATTICE_PROBLEMS]
        + [f"examples/vera/compare-keno-{p}.json" for p in LATTICE_PROBLEMS],
        endfb70_lattices,
    ),
    "endfb70-p3-k": (
        [
            "examples/vera/endfb70/p3/compare-3A.json",
            "examples/vera/endfb70/p3/compare-3B.json",
        ]
        + [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        endfb70_p3_k,
    ),
    "endfb70-p3-distributions": (
        [
            "examples/vera/endfb70/p3/compare-3A.json",
            "examples/vera/endfb70/p3/compare-3B.json",
        ]
        + [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        endfb70_p3_distributions,
    ),
    "endfb70-openmc": (
        [
            "examples/vera/endfb70/openmc/compare-1B.json",
            "examples/vera/endfb70/openmc/compare-2B.json",
        ],
        endfb70_openmc,
    ),
    "depletion-1c-k": (
        [
            "examples/depletion/compare-1C-rings1.json",
            "examples/depletion/compare-1C-rings4.json",
        ],
        depletion_k,
    ),
    "depletion-1c-rings1-densities": (
        ["examples/depletion/compare-1C-rings1.json"],
        depletion_densities,
    ),
    "depletion-1c-rings4-densities": (
        ["examples/depletion/compare-1C-rings4.json"],
        depletion_densities,
    ),
    "depletion-1c-times": (
        [
            "examples/depletion/compare-1C-rings1.json",
            "examples/depletion/compare-1C-rings4.json",
        ],
        depletion_times,
    ),
    "aws-depletion-1c-k": (
        ["examples/depletion/aws-check/compare-aws-1C-rings1.json"],
        aws_depletion_k,
    ),
    "aws-depletion-1c-densities": (
        ["examples/depletion/aws-check/compare-aws-1C-rings1.json"],
        aws_depletion_densities,
    ),
    "depletion-1c-swap-candidates": (
        ["examples/depletion/issue-0133/swap-candidates.json"],
        swap_candidates,
    ),
    **{
        f"depletion-ufg-free-gas-{name}": (
            ["examples/depletion/issue-0133/ufg-free-gas/summary.json"],
            draw,
        )
        for name, draw in [
            ("k", ufg_free_gas_k),
            ("swaps", ufg_free_gas_swaps),
            ("tilt", ufg_free_gas_tilt),
        ]
    },
    **{
        f"depletion-ufg-heavy-free-gas-{name}": (
            ["examples/depletion/issue-0133/ufg-heavy-free-gas/summary.json"],
            draw,
        )
        for name, draw in [
            ("k", ufg_free_gas_k),
            (
                "swaps",
                partial(ufg_free_gas_swaps, specs=UFG_RANGE_SWAPS),
            ),
            (
                "groups",
                partial(
                    ufg_free_gas_swaps,
                    specs=UFG_GROUP_SWAPS,
                    modes=["static", "light", "light-cap", "heavy", "heavy-cap"],
                ),
            ),
            ("tilt-27", partial(ufg_free_gas_tilt, energies=(4.0, 9.877))),
            ("tilt-25", partial(ufg_free_gas_tilt, energies=(15.968, 27.7))),
            ("tilt-24", partial(ufg_free_gas_tilt, energies=(27.7, 48.052))),
        ]
    },
    **{
        f"depletion-ufg-sab-{name}": (
            ["examples/depletion/issue-0133/ufg-sab/summary.json"],
            draw,
        )
        for name, draw in [
            ("k", ufg_free_gas_k),
            ("swaps", partial(ufg_free_gas_swaps, specs=UFG_RANGE_SWAPS)),
            (
                "groups",
                partial(
                    ufg_free_gas_swaps,
                    specs=UFG_GROUP_SWAPS,
                    modes=["heavy", "sab", "heavy-cap", "sab-cap"],
                ),
            ),
            ("transfers", ufg_free_gas_transfers),
            ("tilt-27", partial(ufg_free_gas_tilt, energies=(4.0, 9.877))),
            ("tilt-25", partial(ufg_free_gas_tilt, energies=(15.968, 27.7))),
        ]
    },
    "depletion-group-rows-swaps": (
        ["examples/depletion/issue-0133/group-rows/summary.json"],
        group_rows_swaps,
    ),
    "depletion-group-rows-out-27": (
        ["examples/depletion/issue-0133/group-rows/summary.json"],
        partial(group_rows_out, key="water:27"),
    ),
    "depletion-group-rows-spectrum": (
        ["examples/depletion/issue-0133/group-rows/summary.json"],
        group_rows_spectrum,
    ),
    **{
        f"depletion-ufg-transfer-{name}": (
            ["examples/depletion/issue-0133/ufg-transfer/summary.json"],
            draw,
        )
        for name, draw in [
            ("k", ufg_transfer_k),
            ("swaps", ufg_transfer_swaps),
            ("rows", ufg_transfer_rows),
        ]
    },
    # ADR-0208 の間引いた格子の測定。置き換えの寄与の表は最後の列が最初の形との差で、4 つの形では読めないので描かない
    **{
        f"depletion-ufg-transfer-thin-{name}": (
            ["examples/depletion/issue-0133/ufg-transfer/summary-thin.json"],
            draw,
        )
        for name, draw in [("k", ufg_transfer_k), ("rows", ufg_transfer_rows)]
    },
    "depletion-1c-rings4-rings": (
        ["examples/depletion/compare-1C-rings4.json"],
        depletion_rings,
    ),
    **{
        f"depletion-{job.lower()}-pins": (
            [f"examples/depletion/p2/compare-{job}.json"],
            depletion_assembly_pins,
        )
        for job in ASSEMBLY_DEPLETION_JOBS
    },
    **{
        f"depletion-{job.lower()}-{name}": (
            [f"examples/depletion/p2/compare-{job}.json"],
            draw,
        )
        for job in ASSEMBLY_DEPLETION_JOBS
        for name, draw in [
            ("k", depletion_assembly_k),
            ("densities", depletion_assembly_fuel_densities),
            ("absorbers", depletion_assembly_absorber_densities),
        ]
    },
    "depletion-2a-core-k": (
        ["examples/depletion/branches/compare-core-2A.json"],
        core_depletion_k,
    ),
    "depletion-burnup-interpolation-methods": (
        ["examples/depletion/branches/compare-burnup-interpolation.json"],
        burnup_interpolation_methods,
    ),
    "depletion-branch-cross-terms": (
        ["examples/depletion/branches/compare-branch-cross-terms.json"],
        branch_cross_terms,
    ),
    "depletion-branch-thinning": (
        ["examples/depletion/branches/compare-branch-cross-terms.json"],
        branch_thinning,
    ),
    "depletion-burnup-interpolation-points": (
        ["examples/depletion/branches/compare-burnup-interpolation.json"],
        burnup_interpolation_points,
    ),
    "depletion-2a-core-boron": (
        ["examples/depletion/branches/compare-core-2A-boron.json"],
        core_depletion_boron,
    ),
    "depletion-branch-shielding-k": (
        ["examples/depletion/branches/compare-shielding.json"],
        branch_shielding_k,
    ),
    "depletion-branch-shielding-constants": (
        ["examples/depletion/branches/compare-shielding.json"],
        branch_shielding_constants,
    ),
    "depletion-branch-shielding-boron-k": (
        ["examples/depletion/branches/compare-shielding-boron.json"],
        branch_shielding_k,
    ),
    "depletion-branch-shielding-boron-constants": (
        ["examples/depletion/branches/compare-shielding-boron.json"],
        branch_shielding_constants,
    ),
    **{
        f"depletion-{job.lower()}-free-gas-k": (
            [
                f"examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-{job}-static.json",
                f"examples/depletion/issue-0133/ufg-free-gas-assemblies/compare-{job}-fg.json",
            ],
            depletion_assembly_free_gas_k,
        )
        for job in ("2A", "2O")
    },
    "depletion-assembly-times": (
        [f"examples/depletion/p2/compare-{job}.json" for job in ["2A", "2F", "2N"]],
        depletion_assembly_times,
    ),
}


def locked_md5(root: Path, source: str) -> str:
    """JSON の出力の md5 を、同じ場所の dvc.lock から読む"""
    path = root / source
    lock = path.with_name("dvc.lock")
    if not lock.exists():
        raise TableError(f"{lock.relative_to(root)} does not exist")
    found = re.findall(
        rf"- path: {re.escape(path.name)}\n\s+hash: md5\n\s+md5: ([0-9a-f]{{32}})",
        lock.read_text(),
    )
    if len(found) != 1:
        raise TableError(
            f"{lock.relative_to(root)} has {len(found)} outputs named {path.name} (expected 1)"
        )
    return found[0]


def require_tracked(root: Path, source: str) -> None:
    """JSON が git で追跡されている（ADR-0176）。DVC のキャッシュのまま .gitignore に入ると、worktree や nix の sandbox で読めない。

    .git のない root（nix の sandbox のソースは追跡されたファイルだけから作る）は調べられないので飛ばす。
    """
    if not (root / ".git").exists():
        return
    try:
        found = subprocess.run(
            ["git", "-C", str(root), "ls-files", "--error-unmatch", "--", source],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as e:
        raise TableError(
            f"cannot run git to check that {source} is tracked: {e}"
        ) from e
    if found.returncode == 1 and "did not match" in found.stderr:
        raise TableError(
            f"{source} is not tracked by git (ADR-0176); set cache: false on its outs in the dvc.yaml and git add it"
        )
    if found.returncode != 0:
        raise TableError(
            f"git ls-files failed for {source} (exit {found.returncode}): {found.stderr.strip()}"
        )


def load(root: Path, source: str, verify: bool) -> dict:
    path = root / source
    if not path.exists():
        raise TableError(
            f"{source} does not exist (the stage has not run, or it is not committed); run dvc repro {path.parent.relative_to(root)}/dvc.yaml and commit it"
        )
    if verify:
        require_tracked(root, source)
    if (
        verify
        and source not in GIT_TRACKED
        and hashlib.md5(path.read_bytes()).hexdigest() != locked_md5(root, source)
    ):
        raise TableError(
            f"{source} differs from the md5 in dvc.lock; restore it with git checkout {source}, or run the stage and commit dvc.lock"
        )
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        raise TableError(f"{source} is not valid JSON: {e}") from e


def render(text: str, root: Path, verify: bool = True) -> str:
    """目印の間を、描いた表で置き換えた記録"""
    out: list[str] = []
    name = None
    seen: set[str] = set()
    for number, line in enumerate(text.split("\n"), 1):
        m = OPEN.match(line)
        if not m and line != CLOSE and re.search(r"<!--\s*/?table\b", line):
            raise TableError(f"line {number}: malformed table mark: {line!r}")
        if name is None:
            if line == CLOSE:
                raise TableError(f"line {number}: {CLOSE} without an open mark")
            out.append(line)
            if m:
                name = m.group(1)
                if name not in TABLES:
                    raise TableError(f"line {number}: unknown table {name}")
                if name in seen:
                    raise TableError(f"line {number}: table {name} appears twice")
                seen.add(name)
                sources, draw = TABLES[name]
                try:
                    table = draw([load(root, s, verify) for s in sources])
                except (KeyError, TypeError, IndexError) as e:
                    raise TableError(
                        f"line {number}: table {name}: the JSON lacks {e!r}; rerun the stage that writes it"
                    ) from e
                out.extend(["", table, ""])
        elif line == CLOSE:
            out.append(line)
            name = None
        elif m:
            raise TableError(f"line {number}: table {m.group(1)} opens inside {name}")
    if name is not None:
        raise TableError(f"table {name} is not closed with {CLOSE}")
    return "\n".join(out)


def main(argv: list[str]) -> int:
    check = "--check" in argv
    files = [Path(a) for a in argv if a != "--check"] or sorted(
        ROOT.glob("examples/*/verification*.md")
    )
    status = 0
    for file in files:
        text = file.read_text()
        try:
            new = render(text, ROOT)
        except TableError as e:
            print(f"FAILED: {file}: {e}", file=sys.stderr)
            status = 1
            continue
        if new == text:
            continue
        if check:
            sys.stderr.writelines(
                difflib.unified_diff(
                    text.splitlines(True),
                    new.splitlines(True),
                    f"{file} (record)",
                    f"{file} (drawn)",
                )
            )
            print(
                f"FAILED: {file}: the tables differ from the JSON. Run python3 examples/verification_tables.py {file} and commit it",
                file=sys.stderr,
            )
            status = 1
        else:
            file.write_text(new)
            print(f"updated {file}")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
