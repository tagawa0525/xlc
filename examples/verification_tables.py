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


# 表の名前 → (JSON（リポジトリの直下から）, 描く関数)。JSON は dvc.yaml の段階の outs で、同じ場所の dvc.lock が md5 を持つ
TABLES = {
    "p3-reference-results": (
        [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        p3_reference_results,
    ),
    "p3-openmc-comparison": (
        [
            "examples/vera/p3/compare-3A-openmc.json",
            "examples/vera/p3/compare-3B-openmc.json",
        ],
        p3_openmc_comparison,
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
        ROOT.glob("examples/*/verification.md")
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
