"""#0286 の coverage.py の入力の結果を、変種 `full`（どの核種も自由ガス）と、核種を絞る変種（既定は `share-1e-4`。`--variant`）で
比べ、要約の JSON を書く。

ケースと変種ごとに coverage/ の出力（燃焼は履歴 `<ケース>-<変種>.msgpack`、分岐はテーブル `table-<ケース>-<変種>.msgpack`）と
ログ（`<ケース>-<変種>.log`。水準は `info,xlc_lattice::ufg=debug,xlc_lattice::shielding=debug`）から次を書く:

- 燃焼: 燃焼点ごとの k∞ の反応度の差 1/k(full) − 1/k(share)（pcm）と、主な核種の数密度の燃焼の領域ごとの相対の差（share − full）/ full
  の最大（核種ごと）。燃焼点ごとの少数群定数（履歴の `few_group`）とピン出力（履歴の `pin_powers`）の相対の差の最大
- 分岐: 状態ごとの k∞ の反応度の差、少数群定数（全断面積、吸収、νΣf、散乱行列）の相対の差の最大、燃料の温度と減速材の分岐の
  反応度の変化 1/k(公称) − 1/k(分岐) の差（share − full。pcm）、ミクロ燃焼の核種ごとの微視的断面積（数密度の微分を含む、
  テーブルの `micro` の数）の相対の差の最大
- どちらも: 燃料のリングの U-238 の 1 群の捕獲（ログの `UFG U-238 ring capture`）の相対の差の最大と、自己遮蔽ごとに静止した標的に
  した核種（ログの `UFG free gas rested by the scatter share` のうち割合が 0 でないもの）の一覧の移り変わり（絞る変種のログ）

使い方: python3 compare_coverage.py --json <要約の出力先> [--variant <変種>]（このディレクトリで）
"""

import argparse
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
OUT = HERE / "coverage"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


COVERAGE = load("coverage", HERE / "coverage.py")
CONSTANTS = ("total", "absorption", "nu_fission", "scatter")
# 燃焼の数密度を比べる核種（../compare_depletion.py の主要な核種と Gd）
NUCLIDES = (
    "U-235",
    "U-238",
    "Pu-239",
    "Pu-240",
    "Pu-241",
    "Pu-242",
    "Am-241",
    "Xe-135",
    "Sm-149",
    "Gd-155",
    "Gd-157",
)
FIELD = re.compile(r"(\w+)=(\S+)")
SHIELDING_CLOSE = " INFO shielding: xlc_lattice::shielding: close "


def flat(value) -> list[float]:
    """入れ子の数の並びを 1 列に（文字列と真偽値は除く）"""
    if isinstance(value, dict):
        return [x for k in sorted(value) for x in flat(value[k])]
    if isinstance(value, (list, tuple)):
        return [x for v in value for x in flat(v)]
    if value is None or isinstance(value, (str, bool)):
        return []
    return [float(value)]


def relative(a: list[float], b: list[float]) -> float:
    """`b` に対する `a` の相対の差の最大（`b` が 0 の要素は、`a` も 0 なら差 0、違えば無限大）"""
    if len(a) != len(b):
        sys.exit(f"cannot compare {len(a)} values with {len(b)}")
    if any(math.isnan(x) for x in a + b):
        sys.exit("cannot compare NaN values")
    worst = 0.0
    for x, y in zip(a, b):
        if y == 0.0:
            worst = max(worst, 0.0 if x == 0.0 else math.inf)
        else:
            worst = max(worst, abs(x - y) / abs(y))
    return worst


def reactivity(k_full: float, k_share: float) -> float:
    return (1.0 / k_full - 1.0 / k_share) * 1e5


def read(path: Path):
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def shieldings(path: Path) -> list[dict]:
    """ログの自己遮蔽ごとの、静止した標的にした核種（名前の並び）と、[セル, リング, U-238 の捕獲]。自己遮蔽の span の閉じた行で
    区切る（ピンの並列で 1 つの自己遮蔽の中の順は決まらないので、中はセルとリングで並べる）"""
    blocks, rested, captures = [], set(), []
    for line in path.read_text().splitlines():
        if "UFG free gas rested by the scatter share" in line:
            # 割合 0 は数密度 0 の輸送する核種（燃焼の領域の、まだできていない核種）で、外しても何も変わらないので数えない
            if float(dict(FIELD.findall(line))["share"]) > 0.0:
                rested.add(FIELD_NAME.search(line).group(1))
        elif "UFG U-238 ring capture" in line:
            f = dict(FIELD.findall(line))
            captures.append((int(f["cell"]), int(f["ring"]), float(f["capture"])))
        elif SHIELDING_CLOSE in line:
            blocks.append({"rested": sorted(rested), "captures": sorted(captures)})
            rested, captures = set(), []
    return blocks


FIELD_NAME = re.compile(r"name=\"?([^\" ]+)\"?")


def capture_difference(full: list[dict], share: list[dict], where: str) -> float:
    if len(full) != len(share):
        sys.exit(f"{where}: {len(share)} shieldings against {len(full)}")
    worst = 0.0
    for a, b in zip(share, full):
        # 同じピンを 2 度解く自己遮蔽（分岐の状態の中）はないので、セルとリングで組にする
        if [c[:2] for c in a["captures"]] != [c[:2] for c in b["captures"]]:
            sys.exit(f"{where}: the ring captures are not of the same pins and rings")
        if a["captures"]:
            worst = max(
                worst,
                relative([c[2] for c in a["captures"]], [c[2] for c in b["captures"]]),
            )
    return worst


def rested_changes(blocks: list[dict]) -> list[dict]:
    """自己遮蔽の番号ごとの、前の自己遮蔽から静止した標的に加わった核種と外れた核種（最初の自己遮蔽はすべてを加わったとする）"""
    changes, before = [], set()
    for i, block in enumerate(blocks):
        now = set(block["rested"])
        if now != before or i == 0:
            changes.append(
                {
                    "shielding": i,
                    "count": len(now),
                    "added": sorted(now - before),
                    "removed": sorted(before - now),
                }
            )
        before = now
    return changes


def few_group_difference(full: dict, share: dict, where: str) -> float:
    """燃焼点の少数群定数（履歴の `few_group` は、その点の状態を並べたテーブルと同じ形）の、状態ごとの CONSTANTS の相対の差の最大"""
    if full["points"] != share["points"]:
        sys.exit(f"{where}: the few-group states differ in points")
    return max(
        relative(flat(b[c]), flat(a[c]))
        for a, b in zip(full["states"], share["states"])
        for c in CONSTANTS
    )


def burn(case: str, variant: str) -> dict:
    full = read(OUT / f"{case}-full.msgpack")
    share = read(OUT / f"{case}-{variant}.msgpack")
    if full["nuclides"] != share["nuclides"] or len(full["points"]) != len(
        share["points"]
    ):
        sys.exit(f"{case}: the histories differ in nuclides or points")
    index = {n: i for i, n in enumerate(full["nuclides"])}
    points = []
    for a, b in zip(full["points"], share["points"]):
        if a["burnup"] != b["burnup"]:
            sys.exit(f"{case}: burnup {a['burnup']} against {b['burnup']}")
        densities = {}
        for name in NUCLIDES:
            i = index[name]
            densities[name] = max(
                relative([r[i] for r in b["densities"]], [r[i] for r in a["densities"]])
                if any(r[i] for r in a["densities"])
                else 0.0,
                0.0,
            )
        row = {
            "burnup": a["burnup"],
            "k_inf": a["k_inf"],
            "drho_pcm": reactivity(a["k_inf"], b["k_inf"]),
            "densities": densities,
        }
        if (a["few_group"] is None) != (b["few_group"] is None):
            sys.exit(
                f"{case}: burnup {a['burnup']} has few-group constants in one history only"
            )
        if a["few_group"] is not None:
            row["constants"] = few_group_difference(
                a["few_group"], b["few_group"], f"{case} burnup {a['burnup']}"
            )
        row["pin_powers"] = relative(flat(b["pin_powers"]), flat(a["pin_powers"]))
        points.append(row)
    blocks_full = shieldings(OUT / f"{case}-full.log")
    blocks_share = shieldings(OUT / f"{case}-{variant}.log")
    # CE/CM（ADR-0173）は点ごとに 1 回と、ステップの中点の組成で 1 回解くので、自己遮蔽 2k は点 k、2k + 1 は点 k と k + 1 の中点
    # （出力が一定なので燃焼度でも中点）
    burnups = [p["burnup"] for p in points]
    if len(blocks_share) != 2 * len(burnups) - 1:
        sys.exit(f"{case}: {len(blocks_share)} shieldings for {len(burnups)} points")
    rested = rested_changes(blocks_share)
    for change in rested:
        k, midpoint = divmod(change["shielding"], 2)
        change["burnup"] = (
            0.5 * (burnups[k] + burnups[k + 1]) if midpoint else burnups[k]
        )
    return {
        "points": points,
        "ring_capture": capture_difference(blocks_full, blocks_share, case),
        "rested": rested,
        "seconds": {"full": full["seconds"], variant: share["seconds"]},
    }


def coefficients(table: dict) -> dict:
    """燃焼点ごとの、燃料の温度と減速材の分岐の反応度の変化 1/k(公称) − 1/k(分岐)（pcm）。軸の点が 1 つなら空"""
    axes = {a["name"]: a for a in table["axes"]}
    names = [a["name"] for a in table["axes"]]
    states = {tuple(p): s for p, s in zip(table["points"], table["states"])}
    nominal = []
    for name in names[1:]:
        points = axes[name]["points"]
        nominal.append(points.index(axes[name]["nominal"]))
    out = {}
    for b, burnup in enumerate(axes["burnup"]["points"]):
        base = states.get((b, *nominal))
        if base is None:
            continue
        row = {}
        for j, name in enumerate(names[1:]):
            for p in range(len(axes[name]["points"])):
                if p == nominal[j]:
                    continue
                key = list(nominal)
                key[j] = p
                other = states.get((b, *key))
                if other is not None:
                    row[f"{name}={axes[name]['points'][p]:g}"] = (
                        1.0 / base["k_inf"] - 1.0 / other["k_inf"]
                    ) * 1e5
        out[f"{burnup:g}"] = row
    return out


def micro(full: dict, share: dict, states: list[int] | None = None) -> dict:
    """ミクロ燃焼の核種ごとの、状態（`states` の番号。None はすべて）の微視的断面積（数密度の微分を含む）の相対の差の最大"""
    names = full["micro"]["nuclides"]
    if names != share["micro"]["nuclides"]:
        sys.exit("the micro nuclides differ")
    picked = range(len(full["micro"]["states"])) if states is None else states
    out = {}
    for n, name in enumerate(names):
        a = [x for i in picked for x in flat(full["micro"]["states"][i][n])]
        b = [x for i in picked for x in flat(share["micro"]["states"][i][n])]
        out[name] = relative(b, a)
    rest = {k: v for k, v in full["micro"].items() if k not in ("nuclides", "states")}
    for key, value in rest.items():
        out[f"({key})"] = relative(flat(share["micro"][key]), flat(value))
    return out


def table(case: str, variant: str) -> dict:
    full = read(OUT / f"table-{case}-full.msgpack")
    share = read(OUT / f"table-{case}-{variant}.msgpack")
    if full["points"] != share["points"]:
        sys.exit(f"{case}: the tables differ in points")
    states = []
    for p, a, b in zip(full["points"], full["states"], share["states"]):
        states.append(
            {
                "point": p,
                "k_inf": a["k_inf"],
                "drho_pcm": reactivity(a["k_inf"], b["k_inf"]),
                "constants": max(relative(flat(b[c]), flat(a[c])) for c in CONSTANTS),
            }
        )
    coefficient_full, coefficient_share = coefficients(full), coefficients(share)
    coefficient = {
        burnup: {
            axis: {
                "full_pcm": value,
                "difference_pcm": coefficient_share[burnup][axis] - value,
            }
            for axis, value in row.items()
        }
        for burnup, row in coefficient_full.items()
        if row
    }
    blocks_full = shieldings(OUT / f"{case}-full.log")
    blocks_share = shieldings(OUT / f"{case}-{variant}.log")
    # 燃焼点ごと（テーブルの点の最初の軸）の微視的断面積の差
    burnups = full["axes"][0]["points"]
    micro_by_burnup = {
        f"{burnup:g}": micro(
            full, share, [i for i, p in enumerate(full["points"]) if p[0] == b]
        )
        for b, burnup in enumerate(burnups)
    }
    return {
        "states": states,
        "coefficients": coefficient,
        "micro": micro(full, share),
        "micro_by_burnup": micro_by_burnup,
        "ring_capture": capture_difference(blocks_full, blocks_share, case),
        "rested": rested_changes(blocks_share),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--variant", default="share-1e-4", choices=COVERAGE.VARIANTS)
    args = parser.parse_args()
    if args.variant == "full":
        sys.exit("--variant must name a variant other than full")
    # 変種の入力のあるケース（選び方の変種は COVERAGE.SELECTION_CASES だけ）
    names = sorted(
        n.removesuffix(f"-{args.variant}")
        for n in COVERAGE.inputs()
        if n.endswith(f"-{args.variant}")
    )
    summary = {}
    for case in names:
        summary[case] = (
            burn(case, args.variant)
            if case.startswith("burn-")
            else table(case, args.variant)
        )
    args.json.write_text(json.dumps(summary, indent=1, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
