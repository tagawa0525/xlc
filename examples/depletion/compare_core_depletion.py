"""集合体 1 体の core の燃焼（core の `[depletion]` の履歴、ADR-0190）と、格子計算の燃焼の履歴（ADR-0173）の k を比べる。

使い方: python3 compare_core_depletion.py --json <要約の出力先> <格子計算の履歴（msgpack）> <テーブル（msgpack）> <core の履歴（msgpack）>
        python3 compare_core_depletion.py --report <要約の JSON>

- `--json`: core の履歴が集合体 1 体・軸方向のノード 1 つ（反射境界で k∞ を解くケース）であることと、core の燃焼点の燃焼度が
  格子計算の履歴の燃焼点のどれかと相対 1e-9 で一致することを確かめ、燃焼点ごとに反応度の差 Δρ = (1/k_履歴 − 1/k_core) × 1e5 pcm と、
  その燃焼度がテーブルの燃焼度の点か（点の上は内挿の誤差がない）を書く。テーブルの点では、core の状態（履歴の `state`。テーブルの
  格子点であること）のテーブルの状態の k∞（分岐の格子計算の k∞）との差 (1/k_分岐 − 1/k_core) × 1e5 pcm も書く（ADR-0187 の確認と
  同じ比べ方。分岐と履歴の差を除いた、テーブルの読み書きと内挿の差）。`inputs` は読んだファイルのパスと md5
- `--report`: 要約の JSON を文章で描く（記録の本文の数値はここからコピーする。ADR-0176）
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack

# core の燃焼度と履歴の燃焼点の一致の許容（core は出力密度 × EFPD で、履歴は燃焼度の点から時間を求めて燃焼させる）
BURNUP_TOLERANCE = 1e-9


def branch_k(table: dict, state: dict, burnup: float) -> float:
    """燃焼度 `burnup`（テーブルの点）と core の状態 `state` の、テーブルの状態の k∞（軸の直積の順で、最初の軸がいちばん遅く変わる）"""
    if table["split"] != 1:
        raise SystemExit("the table must not split the assembly")
    wanted = {
        "burnup": burnup,
        "boron": state["boron"],
        "moderator": state["moderator"],
        "fuel": state["fuel"],
    }
    index = 0
    for axis in table["axes"]:
        name = axis["name"]
        if name not in wanted:
            raise SystemExit(
                f"the table has a {name} axis the comparison does not take"
            )
        if name == "boron" and axis["b10"] != state["b10"]:
            raise SystemExit(
                f"the table has b10 {axis['b10']}, the core {state['b10']}"
            )
        points = axis["points"]
        if wanted[name] not in points:
            raise SystemExit(
                f"{name} {wanted[name]} is not a point of the table {points}"
            )
        index = index * len(points) + points.index(wanted[name])
    return table["states"][index]["k_inf"]


def summarize(lattice: dict, table: dict, core: dict) -> dict:
    """3 つのファイルを比べた要約（モジュールの説明の `--json`）。"""
    formats = (lattice.get("format"), table.get("format"), core.get("format"))
    if formats != ("xlc-depletion-history", "xlc-table", "xlc-core-history"):
        raise SystemExit(f"unexpected formats {formats}")
    if "boron" not in core["state"]:
        raise SystemExit(
            "the core history searches the boron; compare it with compare_core_boron.py"
        )
    if len(core["assemblies"]) != 1 or len(core["axial_mesh"]) != 2:
        raise SystemExit(
            "the core must be one assembly of one axial node (k-infinity of a reflected assembly)"
        )
    axes = [a for a in table["axes"] if a["name"] == "burnup"]
    if len(axes) != 1:
        raise SystemExit("the table has no burnup axis")
    table_burnups = axes[0]["points"]
    history = lattice["points"]
    points = []
    for p in core["points"]:
        b = p["burnup"]
        matched = [
            h
            for h in history
            if abs(h["burnup"] - b) <= BURNUP_TOLERANCE * max(1.0, h["burnup"])
        ]
        if len(matched) != 1:
            raise SystemExit(f"core burnup {b} is not a point of the history")
        (h,) = matched
        on_table = h["burnup"] in table_burnups
        branch = branch_k(table, core["state"], h["burnup"]) if on_table else None
        points.append(
            {
                "efpd": p["efpd"],
                "burnup": h["burnup"],
                "core_burnup": b,
                "burnup_relative": 0.0 if h["burnup"] == 0.0 else b / h["burnup"] - 1.0,
                "on_table": on_table,
                "history_k": h["k_inf"],
                "k": p["k"],
                "rho_pcm": (1.0 / h["k_inf"] - 1.0 / p["k"]) * 1e5,
                "branch_k": branch,
                "branch_rho_pcm": None
                if branch is None
                else (1.0 / branch - 1.0 / p["k"]) * 1e5,
            }
        )
    on = [abs(p["rho_pcm"]) for p in points if p["on_table"]]
    between = [abs(p["rho_pcm"]) for p in points if not p["on_table"]]
    return {
        "specific_power": core["specific_power"],
        "table_burnups": table_burnups,
        "points": points,
        "on_table_max_abs_pcm": max(on, default=None),
        "between_max_abs_pcm": max(between, default=None),
        "on_table_branch_max_abs_pcm": max(
            (abs(p["branch_rho_pcm"]) for p in points if p["on_table"]), default=None
        ),
        "burnup_max_relative": max(abs(p["burnup_relative"]) for p in points),
    }


def report(summary: dict) -> list[str]:
    if "inputs" not in summary:
        raise SystemExit("the JSON has no inputs; rerun the comparison")
    points = summary["points"]
    on = [p for p in points if p["on_table"]]
    between = [p for p in points if not p["on_table"]]
    lines = [
        "inputs: "
        + ", ".join(f"{i['path']} {i['md5'][:8]}" for i in summary["inputs"]),
        f"出力密度 {summary['specific_power']:g} W/g、燃焼点 {len(points)}（テーブルの点 {len(on)}、点の間 {len(between)}）",
        f"燃焼度の相対差（core / 履歴 − 1）の最大 {summary['burnup_max_relative']:.1e}",
    ]
    if on:
        lines.append(
            f"テーブルの点の |Δρ| の最大 {summary['on_table_max_abs_pcm']:.4f} pcm（履歴と）、"
            f"{summary['on_table_branch_max_abs_pcm']:.4f} pcm（分岐の格子計算と）"
        )
    if between:
        worst = max(between, key=lambda p: abs(p["rho_pcm"]))
        lines.append(
            f"点の間の Δρ は {min(p['rho_pcm'] for p in between):+.1f}〜{max(p['rho_pcm'] for p in between):+.1f} pcm"
            f"（|Δρ| の最大は {worst['burnup']:g} MWd/kg の {worst['rho_pcm']:+.1f} pcm）"
        )
    lines.append(
        "燃焼度 (MWd/kg) | EFPD | テーブルの点 | k 履歴 | k core | Δρ 履歴 (pcm) | Δρ 分岐 (pcm)"
    )
    for p in points:
        lines.append(
            f"{p['burnup']:g} | {p['efpd']:.4f} | {'○' if p['on_table'] else ''} | {p['history_k']:.8f} | "
            f"{p['k']:.8f} | {p['rho_pcm']:+.4f} | "
            + ("" if p["branch_rho_pcm"] is None else f"{p['branch_rho_pcm']:+.4f}")
        )
    return lines


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        print("\n".join(report(json.loads(Path(args[1]).read_text()))))
    elif len(args) == 5 and args[0] == "--json":
        files = [Path(a) for a in args[2:]]
        lattice, table, core = (msgpack.unpackb(f.read_bytes()) for f in files)
        summary = {
            "format": "xlc-compare-core-depletion",
            # 読んだファイルのパス（実行した cwd から）と md5
            "inputs": [
                {"path": str(f), "md5": hashlib.md5(f.read_bytes()).hexdigest()}
                for f in files
            ],
            **summarize(lattice, table, core),
        }
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        print("\n".join(report(summary)))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
