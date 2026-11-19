"""入力の 1 項目（`ITEMS`）だけを変えた同じ core の燃焼を 2 本ずつ比べ、点と中点の k、ホウ素、燃焼度、出力の割合の差と、固有値計算の数、
外側の反復の回数、時間を出す。

使い方: python3 compare_core_pair.py <項目> --json <要約の出力先> <名前>=<基準の core の入力>,<比べる core の入力> ...
        python3 compare_core_pair.py --report <要約の JSON>

- <項目> は `ITEMS` の鍵（`initial-guess`、`xenon-boron`、`acceleration`）。基準と比べる燃焼の項目の値、要約の `format` と `runs` の鍵は
  その表の行
- 燃焼の履歴は入力の `[depletion]` の `history`（msgpack）、ログは入力の拡張子を `.log` にしたファイル（core の標準エラー。tracing）。
  入力の項目の値が表の基準の値と比べる値の順でないか、2 つの入力が項目と `history` のほかでも違えば止まる（履歴はこれらの項目を
  書かないので、入力で 2 本を見分ける）。同じ入力を別の版の bin/core で解いた燃焼や、出力の名前を付け替えた燃焼は、入力で見分けられ
  ないので比べない
- 2 つの燃焼は、同じ燃焼点（EFPD と日数）、状態、定格の出力密度、軸方向のメッシュ、集合体を持つこと。差は比べる燃焼 − 基準の燃焼
- `--json`: 名前ごとに次を書く
  - `points`: 燃焼点ごとの EFPD、日数、状態と、k とホウ素の差、ノードの燃焼度の相対差の絶対値の最大、出力の割合の相対差の絶対値の
    最大（核分裂のないノードを除く）
  - `midpoints`: ログの中点の行（`depletion midpoint`）ごとの、k とホウ素の差
  - `runs`: 基準と比べる燃焼の、固有値計算（span `nem`）の数、外側の反復の回数の和、固有値計算の時間の和（フィールド `seconds`）、
    core の全体の時間（`core seconds=`）。項目の `per_solve` が真なら、固有値計算 1 回の外側の反復の回数と時間の [最小, 最大]（`per_solve`）
  `inputs` は読んだファイル（入力、履歴、ログの順）のパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import msgpack
from core_history import core_run, differing_keys, flat, largest_relative, setting

# span `nem` の閉じる行のフィールド
NEM = re.compile(r"nem\{iterations=(\d+) seconds=([0-9.eE+-]+)\}")
TOTAL = re.compile(r"core seconds=([0-9.eE+-]+)")
MIDPOINT = re.compile(
    r"depletion midpoint k=([0-9.eE+-]+) boron=([0-9.eE+-]+) solves=(\d+)"
)
# 2 つの燃焼で同じであるべき履歴の項目
SHARED = ("specific_power", "states", "axial_mesh", "assemblies", "heavy_metal")


@dataclass(frozen=True)
class Item:
    """比べる入力の項目。`section` は項目を書く表（None はケースの直下）、`values` は基準と比べる燃焼の値、`default` は書かない
    ときの core の値（crates/core/src/input.rs）、`runs` は要約の `runs` の鍵（基準, 比べる燃焼）"""

    key: str
    section: str | None
    values: tuple[str, str]
    default: str
    format: str
    runs: tuple[str, str]
    per_solve: bool


ITEMS = {
    # 固有値計算の初期値（ADR-0231）: 平らな初期値 → 既定
    "initial-guess": Item(
        "initial-guess",
        "depletion",
        ("flat", "previous"),
        "previous",
        "xlc-compare-core-guess",
        ("flat", "default"),
        False,
    ),
    # 平衡 Xe と臨界ホウ素の収束のさせ方（ADR-0219）: 解き直すたびに探す → 同時の反復
    "xenon-boron": Item(
        "xenon-boron",
        "depletion",
        ("nested", "joint"),
        "joint",
        "xlc-compare-core-xenon-boron",
        ("nested", "joint"),
        False,
    ),
    # 固有値計算の加速（ADR-0232）: 加速なし → CMFD
    "acceleration": Item(
        "acceleration",
        None,
        ("none", "cmfd"),
        "cmfd",
        "xlc-compare-core-acceleration",
        ("none", "cmfd"),
        True,
    ),
}


def run_of(core: dict, log: str) -> dict:
    """燃焼の固有値計算の数、外側の反復の回数と時間の和、全体の時間、燃焼点と中点の固有値計算の回数の和"""
    nem = [(int(i), float(s)) for i, s in NEM.findall(log)]
    total = TOTAL.findall(log)
    if not nem or len(total) != 1:
        raise SystemExit(
            "the log needs the nem spans and one total time (core seconds=)"
        )
    midpoints = [int(m[2]) for m in MIDPOINT.findall(log)]
    solves = sum(p["solves"] for p in core["points"]) + sum(midpoints)
    if solves != len(nem):
        raise SystemExit(
            f"the history and the midpoint lines count {solves} solves, the log {len(nem)} nem spans"
        )
    return {
        "solves": len(nem),
        "iterations": sum(i for i, _ in nem),
        "eigenvalue_seconds": sum(s for _, s in nem),
        "seconds": float(total[0]),
    }


def per_solve(log: str) -> dict:
    """固有値計算 1 回の外側の反復の回数と時間の [最小, 最大]"""
    nem = [(int(i), float(s)) for i, s in NEM.findall(log)]
    iterations = [i for i, _ in nem]
    seconds = [s for _, s in nem]
    return {
        "iterations": [min(iterations), max(iterations)],
        "seconds": [min(seconds), max(seconds)],
    }


def differences(base_core: dict, base_log: str, core: dict, log: str) -> dict:
    """`core` − `base_core` の点と中点の差"""
    for key in SHARED:
        if base_core.get(key) != core.get(key):
            raise SystemExit(f"the histories differ in {key}")
    a, b = core["points"], base_core["points"]
    if [(p["efpd"], p["days"], p["state"]) for p in a] != [
        (p["efpd"], p["days"], p["state"]) for p in b
    ]:
        raise SystemExit("the histories have different points")
    points = [
        {
            "efpd": q["efpd"],
            "days": q["days"],
            "state": q["state"],
            "k": p["k"] - q["k"],
            "boron": p["boron"] - q["boron"],
            "burnup": largest_relative(flat(p["burnups"]), flat(q["burnups"])),
            "power": largest_relative(flat(p["powers"]), flat(q["powers"])),
        }
        for p, q in zip(a, b)
    ]
    mids, base_mids = MIDPOINT.findall(log), MIDPOINT.findall(base_log)
    if len(mids) != len(base_mids):
        raise SystemExit(f"{len(mids)} and {len(base_mids)} midpoints")
    midpoints = [
        {"k": float(m[0]) - float(f[0]), "boron": float(m[1]) - float(f[1])}
        for m, f in zip(mids, base_mids)
    ]
    return {"points": points, "midpoints": midpoints}


def compare(item: Item, base_core: dict, base_log: str, core: dict, log: str) -> dict:
    """`core` − `base_core` の差と、2 つの燃焼の固有値計算の数、反復の回数、時間"""
    runs = {}
    for key, (c, text) in zip(item.runs, ((base_core, base_log), (core, log))):
        runs[key] = run_of(c, text)
        if item.per_solve:
            runs[key]["per_solve"] = per_solve(text)
    return {**differences(base_core, base_log, core, log), "runs": runs}


def pair_runs(
    item: Item, base: Path, other: Path
) -> tuple[dict, str, dict, str, list[Path]]:
    """基準の入力 `base` と比べる入力 `other` の、履歴とログと、読んだファイル。入力の項目の値が `item.values` の順でなければ止まる"""
    loaded: list = []
    files = []
    runs = []
    for source, want, role in zip((base, other), item.values, ("first", "second")):
        run, history, log = core_run(source)
        runs.append(run)
        got = setting(run, item.section, item.key, item.default)
        if got != want:
            raise SystemExit(
                f"{source} sets {item.key} = {got!r}; the {role} input of a pair must set {want!r}"
            )
        loaded += [msgpack.unpackb(history.read_bytes()), log.read_text()]
        files += [source, history, log]
    if files[1].resolve() == files[4].resolve():
        raise SystemExit(
            f"{base} and {other} both write {files[1]}; the later run overwrote the history of the earlier"
        )
    others = differing_keys(*runs, {(item.section, item.key), ("depletion", "history")})
    if others:
        raise SystemExit(
            f"{base} and {other} also differ in {others}; compare runs that change only {item.key}"
        )
    base_core, base_log, core, log = loaded
    return base_core, base_log, core, log, files


def item_of(summary: dict) -> Item:
    for item in ITEMS.values():
        if summary.get("format") == item.format:
            return item
    raise SystemExit(f"not a summary of this tool: {summary.get('format')}")


def ratio(a: float, b: float) -> str:
    return f"{b / a:.3f} 倍"


def span(values: list, digits: str) -> str:
    low, high = values
    return f"{low:{digits}}〜{high:{digits}}"


def report(summary: dict) -> list[str]:
    item = item_of(summary)
    lines = []
    for name, c in summary["cases"].items():
        solves = c["points"] + c["midpoints"]
        k = max(abs(s["k"]) for s in solves)
        boron = max(abs(s["boron"]) for s in solves)
        burnup = max(p["burnup"] for p in c["points"])
        power = max(p["power"] for p in c["points"])
        a, b = (c["runs"][key] for key in item.runs)
        lines.append(
            f"{name}: 点と中点の k の差の絶対値の最大 {k:.2e}、ホウ素 {boron:.2e} ppm、燃焼度の相対差の最大 {burnup:.2e}、"
            f"出力の割合の相対差の最大 {power:.2e}"
        )
        lines.append(
            f"{name}: 固有値計算 {a['solves']} → {b['solves']} 回（{ratio(a['solves'], b['solves'])}）、外側の反復 "
            f"{a['iterations']:,} → {b['iterations']:,} 回（{ratio(a['iterations'], b['iterations'])}）、固有値計算の時間 "
            f"{a['eigenvalue_seconds']:.1f} → {b['eigenvalue_seconds']:.1f} s"
            f"（{ratio(a['eigenvalue_seconds'], b['eigenvalue_seconds'])}）、全体 {a['seconds']:.1f} → {b['seconds']:.1f} s"
        )
        if item.per_solve:
            lines.append(
                f"{name}: 固有値計算 1 回の外側の反復 {span(a['per_solve']['iterations'], ',')} → "
                f"{span(b['per_solve']['iterations'], ',')} 回、時間 {span(a['per_solve']['seconds'], '.3f')} → "
                f"{span(b['per_solve']['seconds'], '.3f')} s"
            )
    return lines


def parse(arg: str) -> tuple[str, list[Path]]:
    name, _, files = arg.partition("=")
    parts = files.split(",")
    paths = [Path(p) for p in parts]
    if not name or len(paths) != 2 or not all(parts):
        raise SystemExit(f"give <name>=<base core input>,<core input>, got {arg}")
    return name, paths


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 4 and args[0] in ITEMS and args[1] == "--json":
        item = ITEMS[args[0]]
        cases = {}
        files = []
        for name, paths in map(parse, args[3:]):
            if name in cases:
                raise SystemExit(f"a second case named {name}")
            base_core, base_log, core, log, read = pair_runs(item, *paths)
            cases[name] = compare(item, base_core, base_log, core, log)
            files += read
        summary = {
            "format": item.format,
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            "cases": cases,
        }
        Path(args[2]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
