"""版ごとの `xs check` の報告（ADR-0266）を、燃焼チェーンの核種とその外に分けて数える（#0292 の全核種を輸送する案の材料）。

使い方: python3 read_survey.py --json <要約の出力先> <版>=<報告（TOML）>,<チェーン（msgpack）>...
        python3 read_survey.py --report <要約の JSON>

- `--json`: 版ごとに、報告のファイルの名前（`n-ZZZ_Sy_AAA[mN].endf`）を (Z, A, 状態) にしてチェーンの核種に対応させ、チェーンの中と外の
  評価の数、読めた数、止まった数（`in_chain`、`outside_chain`）、チェーンの核種の数とファイルの無い核種の数、止まった評価の一覧
  （`stopped`。ファイル、チェーンの名前、理由）を書く。理由は `xs check` の文から先頭の行番号（`line N: `）を除いたもの。
  `reasons` は理由ごと、版ごとの評価（チェーンの名前、なければファイル）の一覧、`watched` は `WATCHED` の核種の版ごとの結果
  （読めれば `readable`、止まれば理由、ファイルが無ければ null）。`inputs` は読んだファイルのパス（ステージの cwd から）と md5。
  要約の文章を標準出力に書く
- `--report`: 要約の JSON を文章で描く
"""

import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack
import tomllib

# 輸送に足すと高燃焼度の負の差が消えた 9 核種（#0292）
WATCHED = [
    "Cm-245",
    "Am-242m",
    "Np-238",
    "Cm-243",
    "Am-242",
    "U-237",
    "Pu-243",
    "U-239",
    "Am-244",
]

FILE_NAME = re.compile(r"n-(\d{3})_[A-Za-z]+_(\d{3})(?:m(\d+))?\.endf")
LINE = re.compile(r"line \d+: ")


def nuclide_key(path: str) -> tuple[int, int, int]:
    """中性子のファイルの名前の (Z, A, 状態)。天然の元素は A = 0"""
    name = Path(path).name
    m = FILE_NAME.fullmatch(name)
    if m is None:
        raise SystemExit(
            f"{name}: not a neutron file name of the form n-ZZZ_Sy_AAA[mN].endf"
        )
    return int(m[1]), int(m[2]), int(m[3] or 0)


def reason(error: str) -> str:
    """`xs check` の理由から、先頭の行番号を除いた文"""
    return LINE.sub("", error, count=1)


def counts(evaluations: list[dict]) -> dict:
    stopped = sum(1 for e in evaluations if "error" in e)
    return {
        "evaluations": len(evaluations),
        "readable": len(evaluations) - stopped,
        "stopped": stopped,
    }


def summarize(versions: list) -> dict:
    """`versions` は (版, 報告, チェーン) の並び"""
    summary = {"format": "xlc-read-survey", "versions": {}, "reasons": {}}
    watched = {n: {} for n in WATCHED}
    for label, report, chain in versions:
        if chain.get("format") != "xlc-chain":
            raise SystemExit(f"{label}: not an xlc chain")
        names = {(n["z"], n["a"], n["state"]): n["name"] for n in chain["nuclides"]}
        evaluations = report["evaluation"]
        keys = [nuclide_key(e["endf"]) for e in evaluations]
        if len(set(keys)) != len(keys):
            raise SystemExit(f"{label}: a nuclide has two files in the report")
        inside = [e for e, k in zip(evaluations, keys) if k in names]
        outside = [e for e, k in zip(evaluations, keys) if k not in names]
        stopped = []
        for e, k in zip(evaluations, keys):
            if "error" not in e:
                continue
            file = Path(e["endf"]).name
            why = reason(e["error"])
            stopped.append({"file": file, "nuclide": names.get(k), "reason": why})
            summary["reasons"].setdefault(why, {}).setdefault(label, []).append(
                names.get(k, file)
            )
        summary["versions"][label] = {
            "in_chain": counts(inside),
            "outside_chain": counts(outside),
            "chain_nuclides": len(names),
            "chain_without_file": len(set(names) - set(keys)),
            "stopped": stopped,
        }
        by_name = {names[k]: e for e, k in zip(evaluations, keys) if k in names}
        for n in WATCHED:
            e = by_name.get(n)
            if e is None:
                watched[n][label] = None
            elif "error" in e:
                watched[n][label] = reason(e["error"])
            else:
                watched[n][label] = "readable"
    summary["watched"] = [{"nuclide": n, "versions": v} for n, v in watched.items()]
    return summary


def report(summary: dict) -> list[str]:
    """要約の JSON を文章にする"""
    lines = []
    if summary.get("inputs"):
        lines.append(
            "inputs: "
            + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"])
        )
    lines.append("")
    lines.append(
        "| 版 | チェーンの評価 | 読めた | 止まった | チェーンの外の評価 | 読めた | 止まった | ファイルの無いチェーンの核種 |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for label, v in summary["versions"].items():
        i, o = v["in_chain"], v["outside_chain"]
        lines.append(
            f"| {label} | {i['evaluations']} | {i['readable']} | {i['stopped']} | {o['evaluations']} | "
            f"{o['readable']} | {o['stopped']} | {v['chain_without_file']} / {v['chain_nuclides']} |"
        )
    lines.append("")
    lines.append("## 止まった理由（版ごとの評価）")
    for why, by_version in summary["reasons"].items():
        lines.append(f"- {why}")
        lines.extend(
            f"  - {label}: {len(labels)}（{' '.join(labels)}）"
            for label, labels in by_version.items()
        )
    lines.append("")
    lines.append("## 注目する 9 核種（#0292）")
    for w in summary["watched"]:
        results = "、".join(
            f"{label} {'ファイルなし' if r is None else r}"
            for label, r in w["versions"].items()
        )
        lines.append(f"- {w['nuclide']}: {results}")
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) >= 3 and args[0] == "--json":
        versions, inputs = [], []
        for arg in args[2:]:
            label, sep, paths = arg.partition("=")
            files = [Path(p) for p in paths.split(",")]
            if not sep or len(files) != 2:
                raise SystemExit(f"{arg}: expected <version>=<report>,<chain>")
            data = [f.read_bytes() for f in files]
            versions.append(
                (
                    label,
                    tomllib.loads(data[0].decode()),
                    msgpack.unpackb(data[1], raw=False),
                )
            )
            # 読んだファイルのパス（ステージの cwd から）と md5。dvc.lock の deps と一致する（record_provenance.py が照合する）
            inputs += [
                {"path": str(f), "md5": hashlib.md5(d).hexdigest()}
                for f, d in zip(files, data)
            ]
        summary = summarize(versions)
        summary["inputs"] = inputs
        Path(args[1]).write_text(
            json.dumps(summary, indent=1, ensure_ascii=False) + "\n"
        )
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
