"""#0293 の項目 8: H-1 の 573.6 K を評価済みの表で解いたブランチ計算（ライブラリ a）と、550 と 600 K から内挿して解いたもの（b）の
k∞ の差（b − a、pcm）を、燃焼点ごとに書く。どちらも同じ履歴の組成（steps/out/<燃料>-40.msgpack）から、減速材だけを 573.6 K に
したテーブル（sab/out/<燃料>-<a|b>.msgpack）。

使い方: python3 compare_sab.py <燃料> [--json <出力>]
"""

import argparse
import json
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent


def burnups(table: dict) -> list[float]:
    """テーブルの状態ごとの燃焼度"""
    axis = [a["name"] for a in table["axes"]].index("burnup")
    points = table["axes"][axis]["points"]
    return [points[p[axis]] for p in table["points"]]


def compare(evaluated: dict, interpolated: dict) -> list[dict]:
    if burnups(evaluated) != burnups(interpolated):
        sys.exit("the two tables solve different burnup points")
    return [
        {
            "burnup": b,
            "k_evaluated": a["k_inf"],
            "dk_pcm": (i["k_inf"] - a["k_inf"]) * 1e5,
        }
        for b, a, i in zip(
            burnups(evaluated), evaluated["states"], interpolated["states"]
        )
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fuel", choices=["u", "mox"])
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    read = lambda v: msgpack.unpackb(
        (HERE / "sab" / "out" / f"{args.fuel}-{v}.msgpack").read_bytes(),
        strict_map_key=False,
    )
    rows = compare(read("a"), read("b"))
    for r in rows:
        print(
            f"{r['burnup']:g} MWd/kg: k {r['k_evaluated']:.6f}, interpolated − evaluated {r['dk_pcm']:+.2f} pcm"
        )
    if args.json:
        args.json.write_text(
            json.dumps({"fuel": args.fuel, "rows": rows}, indent=1) + "\n"
        )


if __name__ == "__main__":
    main()
