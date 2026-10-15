"""検証の記録（examples/*/verification.md）の数値の表を、比較の段階が書いた JSON から描いて置き換える（ADR-0176）。

使い方: python3 examples/verification_tables.py [--check] [<verification.md>...]

- 記録の `<!-- table: <名前> -->` から `<!-- /table -->` までを、名前の表で置き換える。目印の中は丸ごと置き換わる
  （文章は目印の外に書く）。目印の外は触らない。名前と JSON の対応（TABLES）はこのファイルにある
- `--check`: 置き換えず、描いた表と記録の中身が違えば差分を出して 1 で終わる。JSON が無い（段階を実行していない）、dvc.lock の md5 と
  違う（JSON を手で直した、段階を再実行したのに relock していない）、目印が壊れているときも 1 で終わる
- ファイルを渡さないと examples/*/verification.md のすべて
- 標準ライブラリだけで動く（hook と `cargo test` から呼ぶ）
"""

import difflib
import hashlib
import json
import math
import re
import sys
import unicodedata
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


# 表の名前 → (JSON（リポジトリの直下から）, 描く関数)。JSON は dvc.yaml の段階の outs で、同じ場所の dvc.lock が md5 を持つ
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
    "depletion-1c-rings4-rings": (
        ["examples/depletion/compare-1C-rings4.json"],
        depletion_rings,
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


def load(root: Path, source: str, verify: bool) -> dict:
    path = root / source
    if not path.exists():
        raise TableError(
            f"{source} does not exist (the stage has not run, or it is not committed); run dvc repro {path.parent.relative_to(root)}/dvc.yaml and commit it"
        )
    if verify and hashlib.md5(path.read_bytes()).hexdigest() != locked_md5(
        root, source
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
