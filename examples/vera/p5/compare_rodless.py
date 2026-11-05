"""制御棒の分岐を解いた格子（`[branches.rod]`。ADR-0224、ADR-0227）の値 out の状態を、同じ層の棒のない格子の状態と比べる（#0225）。

使い方（examples/vera/p5 で）:
    python3 compare_rodless.py --json <JSON> [rod:]<制御棒の分岐の出力>=<棒のない出力> ...
    python3 compare_rodless.py --report <JSON>

出力の組ごとに、`format` で種類を決める。`rod:` を付けた組の制御棒の分岐の出力は制御棒の軸を持ち（燃料のノードのテーブル）、付けない
組は持たない（反射体のノードのテーブルと媒質。ADR-0227）。食い違えばエラー。

- 分岐のテーブル（`xlc-table`）: 制御棒の軸を最後に持つなら、その最初の値が `out` で、ほかの軸が棒のないテーブルの軸と同じことを
  確かめ、棒のないテーブルの点ごとに、点に値 `out` の番号 0 を足した点の状態（split² の部分ノード）と比べる
  （crates/lattice/tests/history_rod_branch.rs と同じ照らし方）。軸を持たなければ（反射体のノード）、軸と点が同じことを確かめ、
  点ごとに比べる。値 out の点の集合と棒のない点の集合、点と状態の数（split²）が合わなければエラー。`bounds` と `split` は同じでなければならず、`heavy_metal`、`micro`、`histories` は持たない（燃焼のない格子）
- 媒質（`xlc-medium`）: `bounds` と軸が同じことを確かめ、状態を順に比べる

状態は、すべての数がビットまで同じとき（型と、浮動小数点は `struct` の倍精度の表現で比べる）に同じと数える。NaN と無限大は
エラーにする。0.0 と −0.0 は違う状態に数え、相対差は 0 とする。数でない項目の違い（nil と数、空の列と空の map など）はエラー。違う状態があれば、項目（状態の
最上位のキー）ごとの最大の相対差 |a − b| / max(|a|, |b|) のうち最大のものと、テーブルでは k∞ の差の最大 |Δk∞| / k∞ × 1e5（pcm）を出す。

JSON は `{"format": "xlc-compare-rodless", "inputs": [読んだファイルのパスと md5], "outputs": [組ごとの name（制御棒の分岐の出力の
パス）、reference（棒のない出力のパス）、kind（table か medium）、rod_axis、states（比べた状態の数）、bit_equal（ビットまで同じ状態の数）、
max_relative_difference、field（その項目。すべて同じなら null）、max_k_inf_pcm（媒質は null）], "all_bit_equal"}`。
`--json` は標準出力に何も出さず、`--report` が JSON から文章を描く。標準ライブラリと msgpack だけ。
"""

import hashlib
import json
import math
import struct
import sys
from pathlib import Path

import msgpack

FORMAT = "xlc-compare-rodless"
# 燃焼の履歴から作ったテーブルだけが持つキー。p5 の HZP の格子は燃焼が無いので、あれば比べる対象を取り違えている
DEPLETION_KEYS = ("heavy_metal", "micro", "histories")


def read(path: str) -> dict:
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read(), raw=False, strict_map_key=False)


def leaves(x, path: str = ""):
    """入れ子の map と列の末端の (パス, 値)。map のキーは名前の順"""
    if isinstance(x, dict):
        for k in sorted(x, key=str):
            yield from leaves(x[k], f"{path}/{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, f"{path}[{i}]")
    if isinstance(x, (dict, list)):
        # 空の入れ子にも末端を出し、[] と {} を同じ（末端のない）項目としない
        if not x:
            yield path, type(x)
    else:
        yield path, x


def bits(x):
    # == は 0.0 と −0.0 を同じとし、True と 1 を同じとするので、型と倍精度の表現で比べる
    if isinstance(x, float):
        return float, struct.pack("<d", x)
    return type(x), x


def differences(a: dict, b: dict) -> tuple[bool, float, str | None]:
    """2 つの状態がビットまで同じか、最大の相対差と、その最上位のキー"""
    la, lb = dict(leaves(a)), dict(leaves(b))
    if la.keys() != lb.keys():
        raise SystemExit(
            f"the states hold different items: {sorted(la.keys() ^ lb.keys())[:5]}"
        )
    same = True
    worst, field = 0.0, None
    for k, x in la.items():
        y = lb[k]
        for v in (x, y):
            # 定数に NaN や無限大は無いはずで、相対差の最大の比べ方も壊すので止める
            if isinstance(v, float) and not math.isfinite(v):
                raise SystemExit(f"a non-finite value {v} at {k}")
        if bits(x) == bits(y):
            continue
        same = False
        if not (isinstance(x, (int, float)) and isinstance(y, (int, float))):
            raise SystemExit(f"items of different kinds at {k}: {x!r} and {y!r}")
        # 0.0 と −0.0 は符号だけの違いで、差は 0
        scale = max(abs(x), abs(y))
        d = abs(x - y) / scale if scale else 0.0
        if field is None or d > worst:
            worst, field = d, k.split("/")[1].split("[")[0]
    return same, worst, field


def pairs_of_table(rod: dict, rodless: dict, expected: bool) -> list[tuple[dict, dict]]:
    """比べる (制御棒の分岐の状態, 棒のない状態) の組。`expected` は制御棒の軸を持つべきか"""
    for key in ("bounds", "split"):
        if rod[key] != rodless[key]:
            raise SystemExit(f"{key} differs: {rod[key]} and {rodless[key]}")
    for t in (rod, rodless):
        if any(key in t for key in DEPLETION_KEYS):
            raise SystemExit(
                f"a table with depletion data ({[k for k in DEPLETION_KEYS if k in t]})"
            )
    axes = rod["axes"]
    rod_axis = bool(axes) and axes[-1]["name"] == "rod"
    if rod_axis != expected:
        raise SystemExit(
            f"the table {'has' if rod_axis else 'lacks'} the rod axis, expected {'it' if expected else 'none'}"
        )
    if rod_axis:
        if axes[-1]["values"][0] != "out":
            raise SystemExit(
                f"the first rod value is {axes[-1]['values'][0]!r}, not 'out'"
            )
        axes = axes[:-1]
    if axes != rodless["axes"]:
        raise SystemExit(f"the continuous axes differ: {axes} and {rodless['axes']}")
    if not rod_axis and rod["points"] != rodless["points"]:
        raise SystemExit("the points differ")
    n = rod["split"] ** 2
    for t in (rod, rodless):
        if len(t["states"]) != len(t["points"]) * n:
            raise SystemExit(
                f"{len(t['states'])} states for {len(t['points'])} points of {n} nodes"
            )
    if rod_axis:
        outs = sorted(p[:-1] for p in rod["points"] if p[-1] == 0)
        if outs != sorted(rodless["points"]):
            raise SystemExit("the value out points differ from the rodless points")
    index = {tuple(p): i for i, p in enumerate(rod["points"])}
    pairs = []
    for j, point in enumerate(rodless["points"]):
        at = (*point, 0) if rod_axis else tuple(point)
        if at not in index:
            raise SystemExit(f"no point {list(at)} in the rod branch table")
        i = index[at]
        pairs += zip(
            rod["states"][i * n : (i + 1) * n], rodless["states"][j * n : (j + 1) * n]
        )
    return pairs


def compare(rod_path: str, rodless_path: str, rod_axis: bool) -> dict:
    """1 つの組の比較。`rod_axis` は、制御棒の分岐の出力が制御棒の軸を持つべきか（燃料のノードのテーブル。反射体のノードの
    テーブルと媒質は持たない。ADR-0227）"""
    rod, rodless = read(rod_path), read(rodless_path)
    if rod.get("format") != rodless.get("format"):
        raise SystemExit(
            f"{rod_path} is {rod.get('format')} but {rodless_path} is {rodless.get('format')}"
        )
    if rod["format"] == "xlc-table":
        pairs = pairs_of_table(rod, rodless, rod_axis)
    elif rod["format"] == "xlc-medium":
        for key in ("bounds", "axes"):
            if rod[key] != rodless[key]:
                raise SystemExit(f"the media have different {key}")
        if len(rod["states"]) != len(rodless["states"]):
            raise SystemExit("the media have different numbers of states")
        if rod_axis:
            raise SystemExit(f"{rod_path}: a medium has no rod axis")
        pairs = list(zip(rod["states"], rodless["states"]))
    else:
        raise SystemExit(f"{rod_path}: unknown format {rod.get('format')}")
    equal, worst, field = 0, 0.0, None
    k_pcm = 0.0 if rod["format"] == "xlc-table" else None
    for a, b in pairs:
        same, d, f = differences(a, b)
        equal += same
        if f is not None and (field is None or d > worst):
            worst, field = d, f
        if k_pcm is not None:
            k_pcm = max(k_pcm, abs(a["k_inf"] - b["k_inf"]) / b["k_inf"] * 1e5)
    return {
        "name": rod_path,
        "reference": rodless_path,
        "kind": "table" if rod["format"] == "xlc-table" else "medium",
        "rod_axis": rod_axis,
        "states": len(pairs),
        "bit_equal": equal,
        "max_relative_difference": worst,
        "field": field,
        "max_k_inf_pcm": k_pcm,
    }


def md5_inputs(paths: list[str]) -> list[dict]:
    """読んだファイルのパス（ステージの cwd から）と md5"""
    return [
        {"path": p, "md5": hashlib.md5(Path(p).read_bytes()).hexdigest()} for p in paths
    ]


def report(summary: dict) -> list[str]:
    if summary.get("format") != FORMAT:
        raise SystemExit(f"not an {FORMAT} summary")
    lines = []
    for o in summary["outputs"]:
        what = (
            "table with the rod axis"
            if o["rod_axis"]
            else ("table" if o["kind"] == "table" else "medium")
        )
        line = f"{o['name']}: {what}, {o['bit_equal']} of {o['states']} states bit-equal to {o['reference']}"
        if o["field"] is not None:
            line += f"; max relative difference {o['max_relative_difference']:.3e} in {o['field']}"
        if o["max_k_inf_pcm"]:
            line += f", k_inf {o['max_k_inf_pcm']:.3e} pcm"
        lines.append(line)
    lines.append(f"all bit-equal: {'yes' if summary['all_bit_equal'] else 'no'}")
    return lines


def main(argv: list[str]) -> None:
    if len(argv) == 2 and argv[0] == "--report":
        print("\n".join(report(json.loads(Path(argv[1]).read_text()))))
        return
    if len(argv) < 3 or argv[0] != "--json":
        raise SystemExit(__doc__)
    pairs = []
    for arg in argv[2:]:
        rod, sep, rodless = arg.partition("=")
        axis = rod.startswith("rod:")
        rod = rod.removeprefix("rod:")
        if not sep or not rod or not rodless:
            raise SystemExit(
                f"expected [rod:]<rod branch output>=<rodless output>, got {arg!r}"
            )
        pairs.append((rod, rodless, axis))
    outputs = [compare(a, b, axis) for a, b, axis in pairs]
    summary = {
        "format": FORMAT,
        "inputs": md5_inputs([p for a, b, _ in pairs for p in (a, b)]),
        "outputs": outputs,
        "all_bit_equal": all(o["bit_equal"] == o["states"] for o in outputs),
    }
    Path(argv[1]).write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
