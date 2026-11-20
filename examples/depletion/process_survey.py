"""版ごとの中性子の評価のすべてを `xs nuclide` に通す走査（#0292）の入力、報告、要約。

読み込みの走査（read_survey.py、`xs check`。ADR-0266）が確かめない処理（再構成、ドップラー拡がり、非分離共鳴、群定数、MF8 の経路の和）
で止まる評価を、版ごとに洗い出す。

使い方: python3 process_survey.py input <置き場> <走査の報告（TOML）> <チェーン（TOML）> <チェーン（msgpack）> <版の xs の入力> [--reuse <前の置き場>]
          [--declarations <宣言の TOML>]
        python3 process_survey.py report <置き場> [--reuse <前の置き場>]
        python3 process_survey.py summary --json <要約の出力先> <版>=<置き場>,<走査の報告>,<チェーン（msgpack）>...
        python3 process_survey.py --report <要約の JSON>

- `input`: 置き場に xs の入力 `xs.toml`、範囲の入力 `range.toml`（`range_entries`）、`xs nuclide` に通す核種の一覧 `nuclides.txt` を
  書く。走査の評価のすべてを、温度 `TEMPERATURE` の 1 点で輸送する核種にする（走査で止まる評価も、`xs nuclide` が止まる場所を
  確かめるために入れる）。名前はチェーンの名前、チェーンに無いものはファイルの名前（`input_entries`）。S(α, β) の核種は版の xs の
  入力の項目をそのまま使う。`--reuse` は前に流した置き場（`xs.toml`、`xs-range.toml`、混ざったログ `xs.log`。Slurm 999 の形）で、
  そこで同じ項目のまま終わった核種（`reusable`）を `nuclides.txt` から除く。`--declarations` は、評価の欠陥を読み替える核種の項目
  （`[nuclide]` の核種ごとの表。ADR-0272、ADR-0274 の宣言）で、`declare` がその核種の項目に足す
- `report`: 置き場の `log/<核種>.log`（`xs nuclide` の標準エラー）と `log/<核種>.time`（GNU time の `-f "%e s %U user %M kB exit %x"`）、
  無い核種は `--reuse` の置き場のログ（範囲が同じときだけ）から、評価ごとの成否、止まった文、時間、最大の常駐、警告の種類の数を
  `report.json` に書く
- `summary`: 版ごとに、チェーンの中と外の成否、止まった評価の段（走査でも止まれば reading、ほかは processing）と原因（`CAUSES`）、
  宣言（`DECLARATION_KEYS`）で読み替えた評価（走査で止まっても、`corrections` を宣言した `xs nuclide` が通れば通ったと数える）、
  警告の種類ごとの評価の数、時間と常駐の大きい核種、版の間で時間か常駐が 10 倍以上違う核種を JSON に書き、文章を標準出力に書く
- `--report`: 要約の JSON を文章で描く
"""

import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

import msgpack
import tomllib

sys.path.insert(0, str(Path(__file__).resolve().parent))
from read_survey import nuclide_key

# 温度 (K)。VIII.0 を流し直さずに Slurm 999（#0292 の全核種の試行のライブラリ、565 K の 1 点）の結果を使うので、ほかの版もそろえる
TEMPERATURE = 565
# 版の間の時間か常駐の比がこれ以上の核種を挙げる
ACROSS_RATIO = 10
# 版ごとに挙げる、時間と常駐の大きい核種の数
TOP = 5
TIMED_OUT = "timed out"
FREE_NEUTRON = "中性子（Z = 0）の評価で、輸送する標的の核種でない"

# 止まった文（行番号と核種の名前を除いたもの）の正規表現 → (原因の分類, 止まる場所（crates/xs/src から）, issue)。最初に合うもの。
# 場所は止まった文を出す関数（2026-11-20 の xs のコード）
CAUSES: list[tuple[str, str, str, str]] = [
    (
        r"neutron LAW=[67] is not supported",
        "MF6 の中性子の LAW=6、7 を読めない",
        "endf/mf6.rs の parse",
        "#0311",
    ),
    (
        r"last NBT must equal NP",
        "MF5 MT17 の TAB1 の補間の区切りが点の数と食い違う（評価の誤り。入力の `corrections` の宣言で読める。ADR-0272）",
        "endf/tape.rs の tab1",
        "#0316（closed）",
    ),
    (
        r"fission needs MF1 MT452 and MF5 MT18",
        "核分裂の中性子の分布が MF6 MT18 で、輸送する核種の群定数を作れない",
        "process.rs の at_temperature",
        "#0312",
    ),
    (
        r"a charged particle \(ZA=",
        "R 行列の制限形の荷電粒子の対を読めない",
        "resonance/r_matrix_limited.rs の role",
        "#0264",
    ),
    (
        r"the 0 K cross sections are negative",
        "0 K の断面積が負（ADR-0267 でエラー。入力の `negative-to-zero` の宣言で読める）",
        "pointwise.rs の reject_negative",
        "#0309（closed）",
    ),
    (
        r"the MF3 background of MT\d+ at the unresolved grid energies",
        (
            "非分離共鳴（LSSF=0）の MF3 の捕獲か核分裂のバックグラウンドが負（ADR-0274 でエラー。入力の "
            "`unresolved-negative-background` の宣言で読める）"
        ),
        "process.rs の unresolved_shielding",
        "#0319（closed）",
    ),
]
# 止まった段の表示: 読み込みの走査（`xs check`）でも止まるもの、処理で初めて止まるもの
STAGES = {"reading": "走査でも止まる", "processing": "処理で初めて止まる"}

FILE_NAME = re.compile(r"n-\d{3}_([A-Za-z]+)_(\d{3})(?:m(\d+))?\.endf")
LINE = re.compile(r"^line \d+: ")
TEMPERATURE_PREFIX = re.compile(r"^[\d.]+ K: ")
# GNU time の行。混ざったログでは先頭に核種の名前がある。`exit` は核種ごとの時間のファイルだけ
TIME = re.compile(
    r"(?:(?P<name>\S+) )?(?P<e>[\d.]+) s (?P<u>[\d.]+) user (?P<m>\d+) kB(?: exit (?P<x>\d+))?"
)
SIGNAL = re.compile(r"Command terminated by signal (\d+)")
# tracing の行: 時刻、水準、区間（nuclide{name=…}…:）、target、文
LOG = re.compile(
    r"\S+\s+(?P<level>[A-Z]+) (?:nuclide\{name=(?P<name>[^}]+)\}\S*: )?(?P<target>\S+): (?P<message>.*)"
)
# 数は数字から始める（符号と前の点は区切りとして残す。`a..b` の 2 つ目を `.b` と読まない）
NUMBER = re.compile(r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def name_from_file(file: str) -> str:
    """中性子のファイルの名前の核種の名前（最初の準安定状態は `m`。xs の名前の規則）"""
    m = FILE_NAME.fullmatch(Path(file).name)
    if m is None:
        raise SystemExit(f"{file}: not a neutron file name")
    state = int(m[3] or 0)
    suffix = "" if state == 0 else "m" if state == 1 else f"m{state}"
    return f"{m[1]}-{int(m[2])}{suffix}"


def input_entries(
    survey: dict, chain: dict, natural_isotopes: list[str], base: dict
) -> tuple[dict[str, dict], list[tuple[str, str]]]:
    """走査の評価の xs の核種の項目（走査の順）と、入れない評価（ファイル, 理由）。

    天然の元素の評価（A = 0）は、チェーンの `natural-isotopes` の同じ元素の最初の核種の名前で 1 度だけ入れる（ADR-0262 で
    同じ元素の核種は同じ評価を使うので、処理は 1 度で足りる）。"""
    if chain.get("format") != "xlc-chain":
        raise SystemExit("not an xlc chain")
    names = {(n["z"], n["a"], n["state"]): n["name"] for n in chain["nuclides"]}
    z_of = {n["name"]: n["z"] for n in chain["nuclides"]}
    entries, left_out = {}, []
    for e in survey["evaluation"]:
        file = Path(e["endf"]).name
        z, a, state = nuclide_key(file)
        if z == 0:
            left_out.append((file, FREE_NEUTRON))
            continue
        if a == 0:
            listed = [n for n in natural_isotopes if z_of.get(n) == z]
            if not listed:
                raise SystemExit(f"{file}: no natural isotope of Z = {z} is listed")
            name = listed[0]
        else:
            name = names.get((z, a, state)) or name_from_file(file)
        if name in entries:
            raise SystemExit(f"{name} has two evaluations in the survey")
        thermal = base.get(name, {})
        if "thermal" in thermal:
            entries[name] = thermal
        else:
            entries[name] = {"endf": file, "temperatures": [TEMPERATURE]}
    return entries, left_out


def declare(entries: dict[str, dict], declarations: dict[str, dict]) -> dict[str, dict]:
    """`entries` の核種に、宣言のファイルの核種ごとの項目（ADR-0272 の `corrections`、ADR-0274 の非分離共鳴の読み替えなど）を足す。

    入力に無い核種への宣言と、走査の項目を置き換える宣言は止める（宣言は評価の欠陥の読み替えで、温度などの走査の条件は変えない）"""
    unknown = sorted(set(declarations) - set(entries))
    if unknown:
        raise SystemExit(f"declarations for nuclides outside the input: {unknown}")
    out = {}
    for name, entry in entries.items():
        items = declarations.get(name, {})
        replaced = sorted(set(items) & set(entry))
        if replaced:
            raise SystemExit(f"{name}: a declaration would replace {replaced}")
        out[name] = {**entry, **items}
    return out


# 評価の欠陥を読み替える宣言の項目（ADR-0272、ADR-0267、ADR-0274）。報告と要約はこの順に並べる
DECLARATION_KEYS = ("corrections", "negative-to-zero", "unresolved-negative-background")


def declared_items(entry: dict) -> list[str]:
    """xs の核種の項目 `entry` の宣言の項目"""
    return [k for k in DECLARATION_KEYS if k in entry]


def value(v) -> str:
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, list):
        return "[" + ", ".join(value(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{ " + ", ".join(f"{k} = {value(x)}" for k, x in v.items()) + " }"
    return str(v)


def render_input(head: list[str], common: dict, entries: dict) -> str:
    rows = [
        f"{n} = {{ " + ", ".join(f"{k} = {value(v)}" for k, v in e.items()) + " }"
        for n, e in entries.items()
    ]
    return "\n".join(
        [
            *(f"# {h}" for h in head),
            *(f"{k} = {value(v)}" for k, v in common.items()),
            "",
            "[nuclide]",
            *rows,
            "",
        ]
    )


def parse_time(text: str, named: bool = False):
    """GNU time の出力。`named` なら混ざったログの 1 行で、(核種, 値) を返す"""
    lines = [line for line in text.splitlines() if line.strip()]
    m = TIME.fullmatch(lines[-1]) if lines else None
    # 混ざったログの行は名前があって終了コードが無く、核種ごとの時間のファイルはその逆
    if m is None or (m["name"] is not None, m["x"] is None) != (named, named):
        raise SystemExit(f"not a time record: {text!r}")
    t = {
        "seconds": float(m["e"]),
        "user_seconds": float(m["u"]),
        "max_rss_kb": int(m["m"]),
    }
    if named:
        return m["name"], t
    t["exit"] = int(m["x"])
    for line in lines[:-1]:
        s = SIGNAL.fullmatch(line)
        if s is not None:
            t["signal"] = int(s[1])
    return t


def warning_kind(line: str) -> str | None:
    """警告の行の種類（target と、数を N にした文。行番号は除く）。警告でなければ None"""
    m = LOG.fullmatch(line)
    if m is None or m["level"] != "WARN":
        return None
    message = NUMBER.sub("N", LINE.sub("", m["message"], count=1))
    return f"{m['target']}: {message}"


def count_warnings(lines) -> dict[str, int]:
    out: dict[str, int] = {}
    for line in lines:
        kind = warning_kind(line)
        if kind is not None:
            out[kind] = out.get(kind, 0) + 1
    return out


def error_reason(name: str, error: str) -> str:
    """`error: <核種>: line N: <文>` の文"""
    text = error.removeprefix("error: ").removeprefix(f"{name}: ")
    return LINE.sub("", text, count=1)


def evaluation(name: str, endf: str, log: str, time: str) -> dict:
    """1 核種の `xs nuclide` の標準エラーと GNU time の出力から、報告の項目を作る"""
    t = parse_time(time)
    lines = log.splitlines()
    e = {"nuclide": name, "endf": endf, **t, "warnings": count_warnings(lines)}
    errors = [line for line in lines if line.startswith("error: ")]
    if "signal" in t:
        # tracing の行でない最初の行が、止まった理由（Rust の割り当ての失敗の文など）
        first = next((line for line in lines if not LOG.fullmatch(line)), "")
        e["error"] = f"signal {t['signal']} ({first})"
    elif t["exit"] == 124:
        e["error"] = TIMED_OUT
    elif t["exit"] != 0:
        if len(errors) != 1:
            raise SystemExit(f"{name}: exit {t['exit']} with {len(errors)} error lines")
        e["error"] = error_reason(name, errors[0])
    return e


def combined_runs(text: str) -> dict[str, dict]:
    """前に流した結果の混ざったログ（Slurm 999 の xs.log）の、終わった核種ごとの項目（`endf` は無い）。

    核種ごとの終了コードは書かれていないので、エラーの行が無く、核種の区間の close の行と時間の行があるものを終わったとする"""
    lines = text.splitlines()
    if any(line.startswith("error: ") for line in lines):
        raise SystemExit("the combined log has an error line of an unknown nuclide")
    spans: dict[str, list[str]] = {}
    closed, times = set(), {}
    for line in lines:
        m = LOG.fullmatch(line)
        if m is not None and m["name"] is not None:
            spans.setdefault(m["name"], []).append(line)
            if m["target"] == "xs" and m["message"].startswith("close "):
                closed.add(m["name"])
        elif TIME.fullmatch(line):
            name, t = parse_time(line, named=True)
            times[name] = t
    if set(times) != closed:
        raise SystemExit(
            f"runs without both a close line and a time line: {sorted(set(times) ^ closed)}"
        )
    return {
        n: {"nuclide": n, **t, "exit": 0, "warnings": count_warnings(spans.get(n, []))}
        for n, t in times.items()
    }


def cause_of(reason: str, causes) -> tuple:
    for pattern, cause, where, issue in causes:
        if re.search(pattern, reason):
            return cause, where, issue
    return None, None, None


def top(evaluations: list[dict], key: str) -> list[dict]:
    """`key` の大きい順の `TOP` 個の評価の時間と常駐"""
    ranked = sorted(evaluations, key=lambda e: -e[key])[:TOP]
    return [
        {
            "nuclide": e["nuclide"],
            "seconds": e["seconds"],
            "max_rss_kb": e["max_rss_kb"],
        }
        for e in ranked
    ]


def summarize(versions: list, causes) -> dict:
    """`versions` は (版, 報告, 走査の報告, チェーン) の並び"""
    summary = {"format": "xlc-process-survey", "versions": {}, "warnings": {}}
    passed_by: dict[str, dict] = {}
    for label, report, survey, chain in versions:
        if report.get("format") != "xlc-process-survey-report":
            raise SystemExit(f"{label}: not an xlc-process-survey-report")
        if chain.get("format") != "xlc-chain":
            raise SystemExit(f"{label}: not an xlc chain")
        names = {(n["z"], n["a"], n["state"]): n["name"] for n in chain["nuclides"]}
        runs = {e["endf"]: e for e in report["evaluations"]}
        if len(runs) != len(report["evaluations"]):
            raise SystemExit(f"{label}: the report lists a file twice")
        chain_names = set(names.values())
        inside = {"evaluations": 0, "passed": 0, "stopped": 0, "stopped_reading": 0}
        outside = {
            "evaluations": 0,
            "processed": 0,
            "passed": 0,
            "stopped": 0,
            "left_out": 0,
        }
        stopped, declared = [], []
        for s in survey["evaluation"]:
            file = Path(s["endf"]).name
            key = nuclide_key(file)
            run = runs.pop(file, None)
            if run is None:
                if key[0] != 0:
                    raise SystemExit(f"{label}: {file} is not in the report")
                outside["evaluations"] += 1
                outside["left_out"] += 1
                continue
            # 天然の元素の評価は、それを使うチェーンの核種（ADR-0262）の名前で通したので、チェーンの評価に数える
            in_chain = key in names or (key[1] == 0 and run["nuclide"] in chain_names)
            counts = inside if in_chain else outside
            counts["evaluations"] += 1
            if not in_chain:
                outside["processed"] += 1
            if run.get("declarations"):
                declared.append(
                    {"nuclide": run["nuclide"], "items": run["declarations"]}
                )
            if "error" not in run:
                # 読み込みの走査は宣言を持たないので、読む前にテープを直す `corrections` を宣言した評価は、走査で止まっても
                # 通る（ほかの宣言は読んだ後に当てるので、走査の止まりを通せない）
                if "error" in s and "corrections" not in run.get("declarations", []):
                    raise SystemExit(
                        f"{label}: {file} stops in the survey but passes xs nuclide"
                    )
                counts["passed"] += 1
                passed_by.setdefault(run["nuclide"], {})[label] = run
            else:
                counts["stopped"] += 1
                stage = "reading" if "error" in s else "processing"
                # 走査でも止まる評価は、同じ文で止まること（温度ごとの段の文は `<温度> K: ` が前に付く）
                if "error" in s and TEMPERATURE_PREFIX.sub(
                    "", run["error"]
                ) != LINE.sub("", s["error"], count=1):
                    raise SystemExit(
                        f"{label}: {file} stops in the survey with another reason ({s['error']})"
                    )
                if in_chain and stage == "reading":
                    inside["stopped_reading"] += 1
                cause, where, issue = cause_of(run["error"], causes)
                stopped.append(
                    {
                        "nuclide": run["nuclide"],
                        "file": file,
                        "chain": in_chain,
                        "stage": stage,
                        "reason": run["error"],
                        "cause": cause,
                        "where": where,
                        "issue": issue,
                    }
                )
            for kind in run["warnings"]:
                by = summary["warnings"].setdefault(kind, {})
                by[label] = by.get(label, 0) + 1
        if runs:
            raise SystemExit(f"{label}: {sorted(runs)} are not in the survey")

        summary["versions"][label] = {
            "in_chain": inside,
            "outside_chain": outside,
            "stopped": stopped,
            "declared": declared,
            "slowest": top(report["evaluations"], "seconds"),
            "largest": top(report["evaluations"], "max_rss_kb"),
        }
    across = []
    for name, by in passed_by.items():
        if len(by) < 2:
            continue
        seconds = {v: e["seconds"] for v, e in by.items()}
        rss = {v: e["max_rss_kb"] for v, e in by.items()}

        def ratio(values):
            low, high = min(values.values()), max(values.values())
            return math.inf if low <= 0 else high / low

        # 0.1 秒に満たない時間は GNU time の分解能で比が決まらないので、時間の比は 1 秒以上の版を含む核種だけで見る
        slow = max(seconds.values()) >= 1 and ratio(seconds) >= ACROSS_RATIO
        if slow or ratio(rss) >= ACROSS_RATIO:
            across.append({"nuclide": name, "seconds": seconds, "max_rss_kb": rss})
    summary["across"] = across
    return summary


def report_text(summary: dict) -> list[str]:
    lines = []
    if summary.get("inputs"):
        lines.append(
            "inputs: "
            + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
        )
        lines.append("")
    lines.append(
        "| 版 | チェーンの評価 | 通った | 止まった | うち走査でも止まる | チェーンの外の評価 | 通した | 通った | 止まった | 入れない |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for label, v in summary["versions"].items():
        i, o = v["in_chain"], v["outside_chain"]
        lines.append(
            f"| {label} | {i['evaluations']} | {i['passed']} | {i['stopped']} | {i['stopped_reading']} | "
            f"{o['evaluations']} | {o['processed']} | {o['passed']} | {o['stopped']} | {o['left_out']} |"
        )
    lines.append("")
    lines.append("## 止まった評価")
    for label, v in summary["versions"].items():
        for s in v["stopped"]:
            lines.append(
                f"- {label} {s['nuclide']}（{s['file']}、{'チェーン' if s['chain'] else 'チェーンの外'}、{STAGES[s['stage']]}）: "
                f"{s['reason']} — {s['cause'] or '–'}、{s['where'] or '–'}、{s['issue'] or '–'}"
            )
    lines.append("")
    lines.append("## 宣言で読み替えた評価")
    for label, v in summary["versions"].items():
        lines.append(
            f"- {label}: "
            + (
                "、".join(
                    f"{d['nuclide']}（{'、'.join(d['items'])}）"
                    for d in v.get("declared", [])
                )
                or "なし"
            )
        )
    lines.append("")
    lines.append("## 警告の種類（版ごとの評価の数）")
    labels = list(summary["versions"])
    for kind, by in summary["warnings"].items():
        lines.append(f"- {kind}: " + "、".join(f"{v} {by.get(v, 0)}" for v in labels))
    lines.append("")
    lines.append("## 時間と常駐の大きい核種")
    for label, v in summary["versions"].items():
        for key, title in (("slowest", "時間"), ("largest", "常駐")):
            lines.append(
                f"- {label} {title}: "
                + "、".join(
                    f"{e['nuclide']} {e['seconds']} s {e['max_rss_kb']} kB"
                    for e in v[key]
                )
            )
    lines.append("")
    lines.append(f"## 版の間で時間か常駐が {ACROSS_RATIO} 倍以上違う核種")
    for a in summary["across"]:
        lines.append(
            f"- {a['nuclide']}: "
            + "、".join(
                f"{v} {a['seconds'][v]} s {a['max_rss_kb'][v]} kB" for v in a["seconds"]
            )
        )
    return lines


def reusable(entries: dict, earlier: dict, done: set[str]) -> set[str]:
    """前の入力 `earlier` で終わった核種 `done` のうち、項目が今の入力 `entries` と同じもの（同じ処理なので結果を使える）"""
    unknown = done - set(entries)
    if unknown:
        raise SystemExit(f"{sorted(unknown)} of the reused run are not in the input")
    return {n for n in done if earlier.get(n) == entries[n]}


def range_entries(entries: dict, survey: dict) -> dict:
    """範囲（`xs range`）の入力の核種: 走査で止まる評価を除く（`xs range` は MF2 を読み、読めない評価で止まる。ライブラリに入れられる
    評価の範囲にする）"""
    stopped = {Path(e["endf"]).name for e in survey["evaluation"] if "error" in e}
    return {n: e for n, e in entries.items() if e["endf"] not in stopped}


def read_reuse(place: Path, reuse: Path) -> tuple[dict, set[str]]:
    """前に流した置き場 `reuse`（xs.toml、xs.log）の終わった核種の項目と、今の入力で使える核種"""
    text = (reuse / "xs.toml").read_text()
    earlier, now = tomllib.loads(text), tomllib.loads((place / "xs.toml").read_text())
    for key in ("sigma0", "scattering-order"):
        if earlier.get(key) != now.get(key):
            raise SystemExit(f"the reused input has another {key}")
    if (reuse / earlier["endf-directory"]).resolve() != (
        place / now["endf-directory"]
    ).resolve():
        raise SystemExit("the reused input reads another ENDF directory")
    runs = combined_runs((reuse / "xs.log").read_text())
    return runs, reusable(now["nuclide"], earlier["nuclide"], set(runs))


def write_input(args: list[str]) -> None:
    reuse = declarations_path = None
    for flag in ("--reuse", "--declarations"):
        if flag in args[6:]:
            k = args.index(flag, 6)
            if k + 1 >= len(args) or args[k + 1].startswith("--"):
                raise SystemExit(__doc__)
            if flag == "--reuse":
                reuse = Path(args[k + 1])
            else:
                declarations_path = Path(args[k + 1])
            args = args[:k] + args[k + 2 :]
    if len(args) != 6:
        raise SystemExit(__doc__)
    place, survey_path, chain_toml, chain_path, base_path = map(Path, args[1:])
    survey = tomllib.loads(survey_path.read_text())
    chain = msgpack.unpackb(chain_path.read_bytes(), strict_map_key=False)
    natural = tomllib.loads(chain_toml.read_text()).get("natural-isotopes", [])
    base = tomllib.loads(base_path.read_text())
    entries, left_out = input_entries(survey, chain, natural, base["nuclide"])
    if declarations_path is not None:
        declarations = tomllib.loads(declarations_path.read_text())
        entries = declare(entries, declarations.get("nuclide", {}))
    directories = {Path(e["endf"]).parent for e in survey["evaluation"]}
    if len(directories) != 1:
        raise SystemExit(f"the survey reads files of {len(directories)} directories")
    (directory,) = directories
    common = {
        "output": "unused.msgpack",
        "endf-directory": os.path.relpath(survey_path.parent / directory, place),
        "thermal-directory": os.path.relpath(
            base_path.parent / base["thermal-directory"], place
        ),
        "sigma0": base["sigma0"],
        "scattering-order": base["scattering-order"],
    }
    head = [
        f"#0292: {survey_path.name} の評価のすべてを {TEMPERATURE} K の 1 点で輸送する核種にした、`xs nuclide` の走査の入力。",
        (
            "process_survey.py が書く。手で直さない。S(α, β) の核種の項目、σ0、散乱の次数は "
            f"{base_path.name}（{base_path.parent.name}）と同じ。"
        ),
        *(f"入れない評価: {f}（{why}）" for f, why in left_out),
    ]
    if declarations_path is not None:
        head.append(
            f"核種の宣言の項目は {os.path.relpath(declarations_path, place)} から足す（評価の欠陥の読み替え）"
        )
    if reuse is not None:
        head.append(
            f"nuclides.txt は、{reuse} の実行で同じ項目のまま終わった核種を除く（report の --reuse）"
        )
    range_head = [
        f"#0292: xs.toml のうち、走査（{survey_path.name}）で止まらない評価だけの、`xs range` の入力。process_survey.py が書く。",
    ]
    place.mkdir(parents=True, exist_ok=True)
    (place / "xs.toml").write_text(render_input(head, common, entries))
    (place / "range.toml").write_text(
        render_input(range_head, common, range_entries(entries, survey))
    )
    done = read_reuse(place, reuse)[1] if reuse is not None else set()
    (place / "nuclides.txt").write_text(
        "".join(f"{n}\n" for n in entries if n not in done)
    )


def write_report(args: list[str]) -> None:
    reuse = None
    if len(args) == 4 and args[2] == "--reuse":
        reuse, args = Path(args[3]), args[:2]
    if len(args) != 2:
        raise SystemExit(__doc__)
    place = Path(args[1])
    entries = tomllib.loads((place / "xs.toml").read_text())["nuclide"]
    runs, usable = read_reuse(place, reuse) if reuse is not None else ({}, set())
    if (
        reuse is not None
        and (place / "xs-range.toml").read_text()
        != (reuse / "xs-range.toml").read_text()
    ):
        raise SystemExit("the reused run has another pointwise range")
    evaluations = []
    for name, entry in entries.items():
        log, time = place / "log" / f"{name}.log", place / "log" / f"{name}.time"
        if log.exists():
            e = evaluation(name, entry["endf"], log.read_text(), time.read_text())
            e["source"] = "log"
        elif name in usable:
            e = {**runs[name], "endf": entry["endf"], "source": "reused"}
        else:
            raise SystemExit(f"{name} has no log")
        if declared_items(entry):
            e["declarations"] = declared_items(entry)
        evaluations.append(e)
    report = {"format": "xlc-process-survey-report", "evaluations": evaluations}
    if reuse is not None:
        report["reused"] = [
            {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
            for f in (reuse / "xs.toml", reuse / "xs-range.toml", reuse / "xs.log")
        ]
    (place / "report.json").write_text(
        json.dumps(report, indent=1, ensure_ascii=False) + "\n"
    )


def write_summary(args: list[str]) -> None:
    if len(args) < 4 or args[1] != "--json":
        raise SystemExit(__doc__)
    versions, inputs = [], []
    for arg in args[3:]:
        label, sep, paths = arg.partition("=")
        parts = [Path(p) for p in paths.split(",")]
        if not sep or len(parts) != 3:
            raise SystemExit(f"{arg}: expected <version>=<place>,<survey>,<chain>")
        files = [parts[0] / "report.json", parts[1], parts[2]]
        data = [f.read_bytes() for f in files]
        versions.append(
            (
                label,
                json.loads(data[0]),
                tomllib.loads(data[1].decode()),
                msgpack.unpackb(data[2], strict_map_key=False),
            )
        )
        inputs += [
            {"path": str(f), "md5": hashlib.md5(d).hexdigest()}
            for f, d in zip(files, data)
        ]
    summary = summarize(versions, CAUSES)
    summary["inputs"] = inputs
    Path(args[2]).write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")
    print("\n".join(report_text(summary)))


def main() -> None:
    args = sys.argv[1:]
    if args[:1] == ["input"]:
        write_input(args)
    elif args[:1] == ["report"]:
        write_report(args)
    elif args[:1] == ["summary"]:
        write_summary(args)
    elif len(args) == 2 and args[0] == "--report":
        print("\n".join(report_text(json.loads(Path(args[1]).read_text()))))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
