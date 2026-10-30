"""履歴の補正の core（ADR-0215）の確認: R を足す core の集合体 1 体の燃焼を、テーブルと比べる。

    python3 compare_core_history.py --json core-history.json     # このディレクトリの出力から JSON を書く
    python3 compare_core_history.py --report core-history.json   # JSON の要約の文章

記録の表は、この JSON から examples/verification_tables.py が描く。

読むのは、R を持つ公称のテーブル（`lattice histories residual.toml` の出力）、R を書き足す前の公称のテーブル、第 1 項のテーブル（履歴の
分岐の燃焼を条件 j で解いた分岐）、R の比較の要約（`residual.json`。条件 j の R を公称の状態で使った誤差）と、core の 4 つの燃焼の履歴:

- 公称（`core-2A-residual-nominal.toml`。R を持つテーブル、565 K で一生。w = 0）と、同じ燃焼の R を持たないテーブル
  （`core-2A-plain-nominal.toml`）。2 つのファイルがバイトで同じか（`identical`）
- 条件 j（`core-2A-residual-branch.toml`。600 K で一生。w = 1）: R の燃焼度の点ごとに、core の k の第 1 項のテーブルの k∞ からの差と、
  公称の燃焼の core の k の、公称のテーブルの公称の状態の k∞ からの差（ミクロ燃焼の core と格子計算の組成の違いの目安）
- 途中で条件を変えた燃焼（`core-2A-residual-switch.toml`。20 MWd/kg まで 565 K、その後 600 K）: 燃焼点ごとに、ノードの減速材の密度の
  履歴 H、w = (H − H⁰) / (Hᵇ − H⁰)、誤差の上限 w × max |公称の状態で使った誤差|。R の燃焼度の点では、R の効き
  ρ(f + w R) − ρ(f) と、R の点ごとの誤差の上限 w × |公称の状態で使った誤差(u)|

f はノードの定数のうち R の前まで（テーブルのその状態の格子点の定数 + Σᵢ cᵢ (Nᵢ − N̄ᵢ⁰)。Nᵢ は core の数密度。core と同じ係数）で、
R の燃焼度の点の、状態が格子点の燃焼点だけで組み立てる。`reconstruct_pcm` は ρ(f + w R) の 2 群の k と core の k の差（core が R を
足したことの検査。1 ノードの反射境界の k は 2 群の無限媒質の k）。ρ = 1 − 1/k で、差は pcm。
"""

import argparse
import json
from pathlib import Path

from compare_residual import md5, read

FILES = {
    "residual_table": "table-residual-nominal-histories.msgpack",
    "plain_table": "table-residual-nominal.msgpack",
    "branch_table": "table-residual-density-branch.msgpack",
    "nominal": "core-2A-residual-nominal.msgpack",
    "plain": "core-2A-plain-nominal.msgpack",
    "branch": "core-2A-residual-branch.msgpack",
    "switch": "core-2A-residual-switch.msgpack",
    "residual_summary": "residual.json",
}
VARIABLE = "moderator-density"


def compare(table, branch_table, nominal, branch, switch, residual_summary, identical):
    raise NotImplementedError


def report(summary):
    raise NotImplementedError


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--json", help="write the summary to this file")
    parser.add_argument("--report", help="print the report of this summary")
    args = parser.parse_args()
    if args.report:
        with open(args.report, encoding="utf-8") as f:
            print(report(json.load(f)))
        return
    here = Path(__file__).resolve().parent
    with open(here / FILES["residual_summary"], encoding="utf-8") as f:
        residual_summary = json.load(f)
    summary = compare(
        read(here / FILES["residual_table"]),
        read(here / FILES["branch_table"]),
        read(here / FILES["nominal"]),
        read(here / FILES["branch"]),
        read(here / FILES["switch"]),
        residual_summary,
        (here / FILES["nominal"]).read_bytes() == (here / FILES["plain"]).read_bytes(),
    )
    summary["inputs"] = {file: md5(here / file) for file in FILES.values()}
    text = json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
    if args.json:
        Path(args.json).write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
