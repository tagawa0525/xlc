"""分岐の自己遮蔽を公称の減速材で 1 度だけ解いたテーブル（ADR-0201 の `shielding = "nominal-moderator"`）を、状態ごとに解き直した
テーブル（`"every-state"`）と比べる。

使い方: python3 compare_branch_shielding.py --json <要約の出力先> <入力（TOML）>
        python3 compare_branch_shielding.py --report <要約の JSON>

入力の TOML（パスは TOML のあるディレクトリから）:

    [[case]]
    name = "2A"
    every = "table-shielding-2A-every.msgpack"
    nominal = "table-shielding-2A-nominal.msgpack"

- 2 つのテーブルは同じ軸（燃焼度、ホウ素、減速材、燃料）と状態の数を持つこと。状態ごとに軸の値を添える
- k∞ の差は (k_nominal − k_every) × 1e5（pcm）
- 少数群定数の項目（`ITEMS`）ごとに、群（と散乱の元と先、ADF の辺）のうち |every| か |nominal| が 0 でないものの相対差
  |nominal − every| / max(|every|, |nominal|) の最大
- `inputs` は読んだファイルのパス（実行した cwd から）と md5
"""

import hashlib
import json
import sys
from pathlib import Path

import msgpack
import tomllib

# 比べる少数群定数の項目
ITEMS = ["absorption", "nu_fission", "diffusion", "scatter", "adf"]


def flat(value) -> list[float]:
    if isinstance(value, list):
        return [x for v in value for x in flat(v)]
    return [float(value)]


def relative(a, b) -> float:
    """2 つの値（入れ子の list）の、0 でない要素の相対差の最大"""
    x, y = flat(a), flat(b)
    if len(x) != len(y):
        raise SystemExit(f"{len(x)} vs {len(y)} values")
    worst = 0.0
    for u, v in zip(x, y):
        scale = max(abs(u), abs(v))
        if scale > 0.0:
            worst = max(worst, abs(u - v) / scale)
    return worst


def axis_values(axes: list[dict], index: int) -> dict:
    """状態 `index`（軸の直積の順。最初の軸がいちばん遅く変わる）の軸の値"""
    counts = [len(a["points"]) for a in axes]
    at, rest = [], index
    for n in reversed(counts):
        at.append(rest % n)
        rest //= n
    at.reverse()
    by_name = {a["name"]: (a, i) for a, i in zip(axes, at)}
    moderator, m = by_name["moderator"]
    return {
        "burnup": by_name["burnup"][0]["points"][by_name["burnup"][1]],
        "boron": by_name["boron"][0]["points"][by_name["boron"][1]],
        "moderator_temperature": moderator["points"][m],
        "density": moderator["density"][m],
        "fuel": by_name["fuel"][0]["points"][by_name["fuel"][1]],
    }


def compare(name: str, every: dict, nominal: dict) -> dict:
    """1 つの集合体の 2 つのテーブルの、状態ごとの差"""
    if every["axes"] != nominal["axes"]:
        raise SystemExit(f"{name}: the tables have different axes")
    if len(every["states"]) != len(nominal["states"]):
        raise SystemExit(
            f"{name}: {len(every['states'])} vs {len(nominal['states'])} states"
        )
    states = []
    for index, (e, n) in enumerate(zip(every["states"], nominal["states"])):
        states.append(
            {
                **axis_values(every["axes"], index),
                "k_every": e["k_inf"],
                "k_nominal": n["k_inf"],
                "dk_pcm": (n["k_inf"] - e["k_inf"]) * 1e5,
                "relative": {item: relative(e[item], n[item]) for item in ITEMS},
            }
        )
    return {"name": name, "states": states}


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def summarize(spec_path: Path) -> dict:
    spec = tomllib.loads(spec_path.read_text())
    base = spec_path.parent
    cases, inputs = [], []
    for case in spec["case"]:
        tables = {}
        for mode in ["every", "nominal"]:
            path = base / case[mode]
            tables[mode] = msgpack.unpackb(path.read_bytes())
            inputs.append({"path": str(path), "md5": md5(path)})
        cases.append(compare(case["name"], tables["every"], tables["nominal"]))
    return {
        "format": "xlc-compare-branch-shielding",
        "items": ITEMS,
        "cases": cases,
        "inputs": inputs,
    }


def report(summary: dict) -> None:
    for case in summary["cases"]:
        dk = [s["dk_pcm"] for s in case["states"]]
        print(
            f"{case['name']}: {len(dk)} states, k-inf difference {min(dk):+.3f} to {max(dk):+.3f} pcm"
        )
        for s in case["states"]:
            worst = " ".join(
                f"{item} {s['relative'][item]:.2e}" for item in summary["items"]
            )
            print(
                f"  burnup {s['burnup']:g} boron {s['boron']:g} moderator {s['moderator_temperature']:g} K "
                f"{s['density']:g} g/cc fuel {s['fuel']:g} K: k {s['k_every']:.6f} {s['k_nominal']:.6f} "
                f"dk {s['dk_pcm']:+.3f} pcm, {worst}"
            )
    for i in summary["inputs"]:
        print(f"input {i['path']} md5 {i['md5']}")


def main() -> None:
    args = sys.argv[1:]
    if len(args) == 2 and args[0] == "--report":
        report(json.loads(Path(args[1]).read_text()))
    elif len(args) == 3 and args[0] == "--json":
        summary = summarize(Path(args[2]))
        Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
        report(summary)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
