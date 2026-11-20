"""#0292、ADR-0268: UFG の減速計算でポイントワイズの断面積で解く核種（`ufg-nuclides`）を絞ったときの k∞ の差と自己遮蔽の時間。

65 核種を足した試行のライブラリ（trial/library。輸送する核種 154）で燃焼させた U、MOX、Gd のセルの履歴（trial/<燃料>/nominal.msgpack）の
組成から、燃焼させ直さずに、`BURNUP` の燃焼点の基準条件（減速材 0.743 g/cc、565 K、ホウ素 1300 ppm）の 1 状態をブランチ計算で解く。
組（`SETS`）は、足した 65 核種（trial/added-65.txt）のうち UFG で解く核種で、ほかの輸送する核種（もとの 89 核種）と被覆管と減速材の
核種はどの組でも UFG で解く。`all` は `ufg-nuclides` を書かない（すべてを UFG で解く既定）。

全核種の系列（`--series all`）は、データのある全核種を輸送するライブラリ（all/library）で燃焼させた履歴（all/<燃料>/nominal.msgpack）の
組成から、同じ燃焼点を解く。足した核種は p2 のライブラリが輸送しない核種（all/added-to-p2.txt）で、組は `all`、`none`、上の `s17` と、
足した核種の自己遮蔽の目安（dilution_shielding.py の all/dilution-<燃料>.json。燃焼度を通じた最大）が閾値 `THRESHOLDS_PCM` を
どれかの燃料で超える核種の組（`t1` など。閾値が小さいほど核種が増え、前の組を含む）と、同じ閾値で全核種の燃焼の組成から出した目安
（all/dilution-all-<燃料>.json。added_effect.py の吸収の割合から）の組（`a1` など）。対照（control.toml）は、試行のライブラリと
同じ核種だけを輸送する対照のライブラリ（all/control。範囲と xs のビルドは全核種のライブラリと同じ）で、試行の履歴の組成を
すべて UFG で解き、上の試行の系列の `all` との差（ライブラリの範囲と xs の版の分）を出す。

- `python3 ufg_subset.py [--series all] inputs [--rankings t]` は <置き場>/<燃料>/<組>.toml を書く（置き場は `SERIES`、出力は
  <置き場>/<燃料>/out/<組>.msgpack）。`a` の組の目安は `all` を解いたテーブルから出すので、初めは `--rankings t` で書き、目安を出してから
  すべての順位で書き直す
- `python3 ufg_subset.py [--series all] compare --json <出力> <燃料>` は、組ごとの k∞（out/<組>.msgpack）と、自己遮蔽の区間の時間を
  `all` と比べる。差は 1/k(all) − 1/k(組)（pcm）。時間は、同じ入力を排他のジョブ（区分 timed）で解いたログ（out-timed/<組>.log）から
  読む（記録に書く時間は排他のジョブのものだけ。スキル verification-record）。`--untimed <組>...` に挙げた組は、まだ排他のジョブで
  解いていないものとして時間を書かない。全核種の系列は、対照の差 1/k(試行の all) − 1/k(対照)
  （`library_change_pcm`）も書く
- `python3 ufg_subset.py --report <JSON>...` は表を出す
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack
import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
from all_nuclides import transported as transported_names
from inputs import NOMINAL_DENSITY, ZIRCALOY, branch_input
from verification_tables import markdown_table, signed

FUELS = ("u", "mox", "gd")
# 記録の表の燃料の名前
FUEL_NAMES = {"u": "U", "mox": "MOX", "gd": "Gd"}
BURNUP = [0, 10, 30, 60]
# 足した 65 核種のうち UFG で解く核種。順位は examples/depletion/verification.md の「足す核種の自己遮蔽の目安」（燃焼度を通じた
# 最大の自己遮蔽の目安）: 7 pcm 以上の 4 核種と、1 pcm を超える 17 核種
SHIELDED_4 = ["Xe-131", "Pm-147", "Pd-108", "Rh-105"]
SHIELDED_17 = SHIELDED_4 + [
    "Pd-105", "Zr-93", "Pr-141", "Cs-135",
    "Pd-104", "La-139", "Mo-98", "Eu-155", "Mo-97", "Pd-107", "Mo-100", "Nd-148", "Tb-159",
]  # fmt: skip
SETS = {"all": None, "none": [], "s4": SHIELDED_4, "s17": SHIELDED_17}
# 全核種の系列の組の閾値 (pcm)。#0292 の 65 核種の測定（ADR-0268）で、1 pcm を超える 17 核種の組の差が 3.7 pcm 以下だった。
# それより小さい閾値で核種を増やし、差がどこで収束するかを見る
THRESHOLDS_PCM = (1.0, 0.3, 0.1)
# ホウ酸水の式（`water = { ... }`）が置く核種
WATER = ["H-1", "O-16", "B-10", "B-11"]
# 系列ごとの入力の置き場（issue-0292 から）、入力から見たライブラリと履歴、足した核種の一覧（issue-0292 から）、見出しの語
SERIES = {
    "trial": {
        "directory": "ufg-subset",
        "library": "../../trial/library/depletion.msgpack",
        "history": "../../trial/{fuel}/nominal.msgpack",
        "added": "trial/added-65.txt",
        "label": "",
        # 組は `SETS`
        "sets": "fixed",
    },
    "all": {
        "directory": "all-ufg",
        "library": "../../all/library/depletion.msgpack",
        "history": "../../all/{fuel}/nominal.msgpack",
        "added": "all/added-to-p2.txt",
        "label": "全核種の ",
        # 組は自己遮蔽の目安の閾値から（`sets`）。対照（control.toml）も書いて比べる
        "sets": "dilution",
        "control": True,
        # 輸送する核種をライブラリの入力から取る（燃焼の前に入力を書くため。履歴の輸送する核種は、チェーンの核種のうち
        # ライブラリが燃焼させる物質の温度で持つもので、入力のどの核種もチェーンの評価なので、その温度の点を持つ入力の核種と
        # 同じ。compare が履歴と食い違えば止める。`burnt_transported`）
        "library_input": "all/library/xs.toml",
    },
}
CONTROL_LIBRARY = "../../all/control/depletion.msgpack"


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def structural() -> list[str]:
    """被覆管と減速材の核種（中心でないリングの核種で、`ufg-nuclides` に挙げる）"""
    names = re.findall(r"^([A-Z][a-z]?-\d+m?) =", ZIRCALOY, re.MULTILINE)
    return names + [n for n in WATER if n not in names]


def listed(
    transported: list[str], added: list[str], kept: list[str] | None
) -> list[str] | None:
    """組の `ufg-nuclides`: 輸送する核種のうち、足した核種で `kept` に無いものを除き、被覆管と減速材の核種を足したもの"""
    if kept is None:
        return None
    unknown = [n for n in kept if n not in added]
    if unknown:
        raise SystemExit(f"{unknown} are not added nuclides")
    missing = [n for n in added if n not in transported]
    if missing:
        raise SystemExit(f"the added nuclides {missing} are not transported")
    names = [n for n in transported if n not in added or n in kept]
    return names + [n for n in structural() if n not in names]


def branch_text(fuel: str, case: str, head: str, library: str, history: str) -> str:
    """`history` の組成から `BURNUP` を解く基準条件の 1 状態のブランチ計算（ライブラリだけを `library` に替える）"""
    text = branch_input(fuel, f"out/{case}", history, BURNUP, None, NOMINAL_DENSITY)
    old = 'library = "../../../p2/depletion.msgpack"'
    if old not in text:
        raise SystemExit("the branch input names another library")
    body = text.replace(old, f'library = "{library}"')
    return head + body.split("\n", 1)[1]


def case_input(
    fuel: str, case: str, names: list[str] | None, kept: int, series: dict
) -> str:
    """組 `case` の入力（`names` は `ufg-nuclides`、`kept` は UFG で解く足した核種の数）"""
    head = (
        f"# #0292、ADR-0268: {series['label']}{fuel} の組 {case}（UFG で解く足した核種: "
        f"{'すべて' if names is None else kept}）。ufg_subset.py が書く。手で直さない\n"
    )
    text = branch_text(
        fuel, case, head, series["library"], series["history"].format(fuel=fuel)
    )
    if names is None:
        return text
    line = f'library = "{series["library"]}"'
    listing = "ufg-nuclides = [" + ", ".join(f'"{n}"' for n in names) + "]"
    return text.replace(line, f"{line}\n{listing}", 1)


def control_input(fuel: str) -> str:
    head = (
        f"# #0292: {fuel} の対照（試行のライブラリと同じ核種を輸送する ../../all/control で、試行の履歴の組成をすべて UFG で解く）。"
        "ufg_subset.py が書く。手で直さない\n"
    )
    return branch_text(
        fuel,
        "control",
        head,
        CONTROL_LIBRARY,
        SERIES["trial"]["history"].format(fuel=fuel),
    )


def threshold_sets(
    summaries: list[dict], thresholds, prefix: str = "t"
) -> dict[str, list[str]]:
    """閾値ごとの組（名前は `prefix` と閾値）: 自己遮蔽の目安（燃焼度を通じた最大）の絶対値が、どれかの燃料で閾値を超える核種を、
    燃料を通じた最大の絶対値の大きい順に"""
    largest = {}
    for d in summaries:
        for r in d["nuclides"]:
            largest[r["nuclide"]] = max(
                largest.get(r["nuclide"], 0.0), abs(r["max_pcm"])
            )
    ranked = sorted(largest, key=lambda n: -largest[n])
    return {f"{prefix}{t:g}": [n for n in ranked if largest[n] > t] for t in thresholds}


# 自己遮蔽の目安の順位ごとの JSON の名前と組の名前の頭（`t` は p2 の 89 核種の燃焼の組成から、`a` は全核種の燃焼の組成から。
# 89 核種の燃焼では足した核種が捕獲しないので、その捕獲の先の核種（Pm-148、Dy-160〜164 など）がほとんど無い）
RANKINGS = {"t": "dilution-{fuel}.json", "a": "dilution-all-{fuel}.json"}


def ranking_files(root: Path, rankings) -> dict[str, list[Path]]:
    """順位の頭 → 燃料の順の自己遮蔽の目安の JSON"""
    unknown = [r for r in rankings if r not in RANKINGS]
    if unknown:
        raise SystemExit(f"{unknown} are not rankings ({list(RANKINGS)})")
    return {
        r: [root / "all" / RANKINGS[r].format(fuel=fuel) for fuel in FUELS]
        for r in rankings
    }


def set_sources(series: str, root: Path, rankings=tuple(RANKINGS)) -> list[Path]:
    """系列の組を決めるファイル（自己遮蔽の目安の JSON。順位ごとに燃料の順）"""
    if SERIES[series]["sets"] == "fixed":
        return []
    return [p for paths in ranking_files(root, rankings).values() for p in paths]


def sets(series: str, root: Path, rankings=tuple(RANKINGS)) -> dict:
    """系列の組: 名前 → UFG で解く足した核種（None はすべて）。`rankings` は組を作る順位の頭（`a` の目安は `t` の系列の
    `all` を解いた後に出すので、その前の入力は `t` だけで書く）"""
    if SERIES[series]["sets"] == "fixed":
        return SETS
    out = {"all": None, "none": [], "s17": SHIELDED_17}
    for prefix, paths in ranking_files(root, rankings).items():
        summaries = [json.loads(p.read_text()) for p in paths]
        out.update(threshold_sets(summaries, THRESHOLDS_PCM, prefix))
    return out


def history_file(root: Path, config: dict, fuel: str) -> Path:
    """系列の燃料の履歴（入力の `[history]` の `file` は、入力の置き場 <置き場>/<燃料> から見たパス）"""
    return root / config["directory"] / fuel / config["history"].format(fuel=fuel)


def library_transported(path: Path) -> list[str]:
    """xs の入力の輸送する核種（放射化の核種でないもの）"""
    return transported_names(tomllib.loads(path.read_text())["nuclide"])


def burnt_transported(root: Path, config: dict, fuel: str) -> list[str]:
    """燃焼の履歴が輸送する核種として持つはずの核種: xs の入力の輸送する核種のうち、燃焼させるすべての物質の温度の点を持つもの
    （燃焼は、その温度のデータの無い核種を崩壊だけで追う。crates/lattice/src/depletion.rs の `has_data_at`。H-1 は S(α, β) の
    内挿の点だけを持ち、燃料の温度の点を持たない）。燃焼の入力は履歴の隣の同じ名前の .toml"""
    burnup = tomllib.loads(
        history_file(root, config, fuel).with_suffix(".toml").read_text()
    )
    temperatures = {
        burnup["material"][m]["temperature"] for m in burnup["depletion"]["materials"]
    }
    entries = tomllib.loads((root / config["library_input"]).read_text())["nuclide"]
    return [
        n
        for n in transported_names(entries)
        if temperatures <= set(entries[n]["temperatures"])
    ]


def write_inputs(series: str, rankings=tuple(RANKINGS)) -> None:
    config = SERIES[series]
    added = (HERE / config["added"]).read_text().split()
    cases = sets(series, HERE, rankings)
    for fuel in FUELS:
        if "library_input" in config:
            transported = library_transported(HERE / config["library_input"])
        else:
            transported = read(history_file(HERE, config, fuel))["transportable"]
        directory = HERE / config["directory"] / fuel
        (directory / "out").mkdir(parents=True, exist_ok=True)
        for case, kept in cases.items():
            names = listed(transported, added, kept)
            count = 0 if kept is None else len(kept)
            (directory / f"{case}.toml").write_text(
                case_input(fuel, case, names, count, config)
            )
        if config.get("control"):
            (directory / "control.toml").write_text(control_input(fuel))


SECONDS = {"s": 1.0, "ms": 1e-3, "µs": 1e-6, "us": 1e-6, "ns": 1e-9}


def shielding_seconds(log: str) -> list[float]:
    """ログの自己遮蔽の区間（`shielding` の span）の時間 (s)。並びは解いた順"""
    out = []
    for m in re.finditer(
        r"INFO shielding: \S+: close time\.busy=([0-9.]+)(s|ms|µs|us|ns) ", log
    ):
        out.append(float(m.group(1)) * SECONDS[m.group(2)])
    return out


def k_by_burnup(table: dict) -> dict[float, float]:
    axes = table["axes"]
    if axes[0]["name"] != "burnup" or any(len(a["points"]) != 1 for a in axes[1:]):
        raise SystemExit("not a table of one state per burnup point")
    burnup = axes[0]["points"]
    if len(table["points"]) != len(burnup):
        raise SystemExit("the table has another number of states than burnup points")
    return {burnup[p[0]]: s["k_inf"] for p, s in zip(table["points"], table["states"])}


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def kept_count(path: Path, added: set[str]) -> int | None:
    """解いた入力の `ufg-nuclides` のうち足した核種の数（一覧が無ければ None。すべてを UFG で解く）"""
    names = tomllib.loads(path.read_text()).get("ufg-nuclides")
    return None if names is None else sum(n in added for n in names)


def compare(
    fuel: str, root: Path, series: str, cases_of=None, untimed=frozenset()
) -> dict:
    """`cases_of` は系列の組の名前（省けば `sets`。組の核種の数は、組の名前の入力から数える）。`untimed` は排他のジョブで
    まだ解いていない組で、時間を書かない（基準の `all` は時間の比の分母なので挙げられない）"""
    if "all" in untimed:
        raise SystemExit("the reference set all must be timed")
    config = SERIES[series]
    inputs = []
    if cases_of is None:
        cases_of = sets(series, root)
        inputs += [
            {"path": str(p.relative_to(root.parent)), "md5": md5(p)}
            for p in set_sources(series, root)
        ]
    stray = sorted(set(untimed) - set(cases_of))
    if stray:
        raise SystemExit(f"the untimed {stray} are not sets of the series")
    added = set((root / config["added"]).read_text().split())
    directory = root / config["directory"] / fuel
    cases = {}
    for case in cases_of:
        table_path = directory / "out" / f"{case}.msgpack"
        log_path = directory / "out-timed" / f"{case}.log"
        k = k_by_burnup(read(table_path))
        input_path = directory / f"{case}.toml"
        paths = [input_path, table_path]
        if case in untimed:
            seconds = None
        else:
            seconds = shielding_seconds(log_path.read_text())
            if len(seconds) != len(k):
                raise SystemExit(
                    f"{log_path}: {len(seconds)} shieldings for {len(k)} states"
                )
            paths.append(log_path)
        cases[case] = {
            "k": k,
            "seconds": seconds,
            "kept": kept_count(input_path, added),
        }
        for path in paths:
            inputs.append(
                {"path": str(path.relative_to(root.parent)), "md5": md5(path)}
            )
    reference = cases["all"]
    burnup = sorted(reference["k"])
    rows = []
    for case, c in cases.items():
        if sorted(c["k"]) != burnup:
            raise SystemExit(f"{case}: other burnup points than all")
        rows.append(
            {
                "case": case,
                "ufg_added": c["kept"],
                "k_inf": [c["k"][b] for b in burnup],
                "pcm": [
                    1e5 * (1.0 / reference["k"][b] - 1.0 / c["k"][b]) for b in burnup
                ],
                "seconds": c["seconds"],
                "ratio": None
                if c["seconds"] is None
                else [s / r for s, r in zip(c["seconds"], reference["seconds"])],
            }
        )
    out = {
        "format": "xlc-compare-ufg-subset",
        "fuel": fuel,
        "burnup": burnup,
        "cases": rows,
        "inputs": inputs,
    }
    if "library_input" in config:
        history = read(history_file(root, config, fuel))
        if sorted(history["transportable"]) != sorted(
            burnt_transported(root, config, fuel)
        ):
            raise SystemExit(
                "the history transports other nuclides than the library input "
                "has at the burnt temperatures"
            )
    if config.get("control"):
        paths = [
            root / SERIES["trial"]["directory"] / fuel / "out" / "all.msgpack",
            directory / "out" / "control.msgpack",
        ]
        trial, control = (k_by_burnup(read(p)) for p in paths)
        if sorted(trial) != burnup or sorted(control) != burnup:
            raise SystemExit("the control has other burnup points than all")
        out["library_change_pcm"] = [
            1e5 * (1.0 / trial[b] - 1.0 / control[b]) for b in burnup
        ]
        # 2 つのテーブルと、それを出した入力
        paths += [
            root / SERIES["trial"]["directory"] / fuel / "all.toml",
            directory / "control.toml",
        ]
        inputs += [
            {"path": str(p.relative_to(root.parent)), "md5": md5(p)} for p in paths
        ]
    return out


def burnup_header(summaries: list[dict]) -> list[str]:
    burnup = summaries[0]["burnup"]
    if any(d["burnup"] != burnup for d in summaries):
        raise SystemExit("the fuels have other burnup points")
    return [f"{b:g}" for b in burnup]


def k_table(summaries: list[dict]) -> str:
    """k∞ の差（pcm）の表。行は燃料と組（基準の `all` を除く）、列は燃焼度 (MWd/kg)"""
    rows = [
        [f"{FUEL_NAMES[d['fuel']]} {c['case']}", *(signed(x, 1) for x in c["pcm"])]
        for d in summaries
        for c in d["cases"]
        if c["case"] != "all"
    ]
    return markdown_table(["燃料と組", *burnup_header(summaries)], rows)


def time_table(summaries: list[dict]) -> str:
    """自己遮蔽の時間 (s) と、括弧の中に `all` に対する比の表（排他のジョブで解いていない組は —）"""
    rows = []
    for d in summaries:
        width = len(d["burnup"])
        for c in d["cases"]:
            if c["seconds"] is None:
                cells = ["—"] * width
            else:
                cells = [
                    f"{s:.2f}（{r:.2f}）" for s, r in zip(c["seconds"], c["ratio"])
                ]
            rows.append([f"{FUEL_NAMES[d['fuel']]} {c['case']}", *cells])
    return markdown_table(["燃料と組", *burnup_header(summaries)], rows)


def library_table(summaries: list[dict]) -> str:
    """対照の差 1/k(試行の all) − 1/k(対照)（pcm。ライブラリの範囲と xs の版の分）の表。行は燃料、列は燃焼度 (MWd/kg)"""
    rows = [
        [FUEL_NAMES[d["fuel"]], *(signed(x, 1) for x in d["library_change_pcm"])]
        for d in summaries
    ]
    return markdown_table(["燃料", *burnup_header(summaries)], rows)


def report(paths: list[Path]) -> None:
    summaries = [json.loads(path.read_text()) for path in paths]
    print(k_table(summaries))
    print()
    print(time_table(summaries))
    controlled = [d for d in summaries if "library_change_pcm" in d]
    if controlled:
        print()
        print(library_table(controlled))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--report", nargs="+", type=Path)
    parser.add_argument("--series", choices=list(SERIES), default="trial")
    sub = parser.add_subparsers(dest="command")
    i = sub.add_parser("inputs")
    i.add_argument(
        "--rankings",
        nargs="+",
        default=list(RANKINGS),
        help="組を作る順位の頭（既定はすべて。a の目安を出す前は t だけ）",
    )
    c = sub.add_parser("compare")
    c.add_argument("--json", type=Path, required=True)
    c.add_argument("fuel", choices=FUELS)
    c.add_argument(
        "--untimed",
        nargs="+",
        default=[],
        help="排他のジョブでまだ解いていない組（時間を書かない）",
    )
    args = parser.parse_args()
    if args.report:
        report(args.report)
    elif args.command == "inputs":
        write_inputs(args.series, tuple(args.rankings))
    elif args.command == "compare":
        summary = compare(args.fuel, HERE, args.series, untimed=frozenset(args.untimed))
        args.json.write_text(json.dumps(summary, indent=1) + "\n")
    else:
        parser.error("give inputs, compare or --report")


if __name__ == "__main__":
    main()
