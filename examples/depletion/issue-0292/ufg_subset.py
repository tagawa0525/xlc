"""#0292、ADR-0268: UFG の減速計算でポイントワイズの断面積で解く核種（`ufg-nuclides`）を絞ったときの k∞ の差と自己遮蔽の時間。

65 核種を足した試行のライブラリ（trial/library。輸送する核種 154）で燃焼させた U、MOX、Gd のセルの履歴（trial/<燃料>/nominal.msgpack）の
組成から、燃焼させ直さずに、`BURNUP` の燃焼点の基準条件（減速材 0.743 g/cc、565 K、ホウ素 1300 ppm）の 1 状態をブランチ計算で解く。
組（`SETS`）は、足した 65 核種（trial/added-65.txt）のうち UFG で解く核種で、ほかの輸送する核種（もとの 89 核種）と被覆管と減速材の
核種はどの組でも UFG で解く。`all` は `ufg-nuclides` を書かない（すべてを UFG で解く既定）。

- `python3 ufg_subset.py inputs` は ufg-subset/<燃料>/<組>.toml を書く（出力は ufg-subset/<燃料>/out/<組>.msgpack）
- `python3 ufg_subset.py compare --json ufg-subset-<燃料>.json <燃料>` は、組ごとの k∞（out/<組>.msgpack）と、自己遮蔽の区間の時間を
  `all` と比べる。差は 1/k(all) − 1/k(組)（pcm）。時間は、同じ入力を排他のジョブ（区分 timed）で解いたログ（out-timed/<組>.log）から
  読む（記録に書く時間は排他のジョブのものだけ。スキル verification-record）
- `python3 ufg_subset.py --report <JSON>...` は表を出す
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent))
from inputs import NOMINAL_DENSITY, ZIRCALOY, branch_input
from verification_tables import markdown_table, signed

FUELS = ("u", "mox", "gd")
BURNUP = [0, 10, 30, 60]
# 足した 65 核種のうち UFG で解く核種。順位は examples/depletion/verification.md の「足す核種の自己遮蔽の目安」（燃焼度を通じた
# 最大の自己遮蔽の目安）: 7 pcm 以上の 4 核種と、1 pcm を超える 17 核種
SHIELDED_4 = ["Xe-131", "Pm-147", "Pd-108", "Rh-105"]
SHIELDED_17 = SHIELDED_4 + [
    "Pd-105", "Zr-93", "Pr-141", "Cs-135",
    "Pd-104", "La-139", "Mo-98", "Eu-155", "Mo-97", "Pd-107", "Mo-100", "Nd-148", "Tb-159",
]  # fmt: skip
SETS = {"all": None, "none": [], "s4": SHIELDED_4, "s17": SHIELDED_17}
# ホウ酸水の式（`water = { ... }`）が置く核種
WATER = ["H-1", "O-16", "B-10", "B-11"]
LIBRARY = "../../trial/library/depletion.msgpack"


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


def case_input(fuel: str, case: str, names: list[str] | None) -> str:
    text = branch_input(
        fuel,
        f"out/{case}",
        f"../../trial/{fuel}/nominal.msgpack",
        BURNUP,
        None,
        NOMINAL_DENSITY,
    )
    old = 'library = "../../../p2/depletion.msgpack"'
    if old not in text:
        raise SystemExit("the branch input names another library")
    head = (
        f"# #0292、ADR-0268: {fuel} の組 {case}（UFG で解く足した核種: "
        f"{'すべて' if names is None else len(SETS[case])}）。ufg_subset.py が書く。手で直さない\n"
    )
    line = f'library = "{LIBRARY}"'
    if names is not None:
        line += "\nufg-nuclides = [" + ", ".join(f'"{n}"' for n in names) + "]"
    body = text.replace(old, line)
    return head + body.split("\n", 1)[1]


def write_inputs() -> None:
    added = (HERE / "trial" / "added-65.txt").read_text().split()
    for fuel in FUELS:
        transported = read(HERE / "trial" / fuel / "nominal.msgpack")["transportable"]
        directory = HERE / "ufg-subset" / fuel
        (directory / "out").mkdir(parents=True, exist_ok=True)
        for case, kept in SETS.items():
            names = listed(transported, added, kept)
            (directory / f"{case}.toml").write_text(case_input(fuel, case, names))


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


def compare(fuel: str, root: Path) -> dict:
    directory = root / "ufg-subset" / fuel
    cases = {}
    inputs = []
    for case in SETS:
        table_path = directory / "out" / f"{case}.msgpack"
        log_path = directory / "out-timed" / f"{case}.log"
        k = k_by_burnup(read(table_path))
        seconds = shielding_seconds(log_path.read_text())
        if len(seconds) != len(k):
            raise SystemExit(
                f"{log_path}: {len(seconds)} shieldings for {len(k)} states"
            )
        cases[case] = {"k": k, "seconds": seconds}
        for path in (table_path, log_path):
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
                "ufg_added": None if SETS[case] is None else len(SETS[case]),
                "k_inf": [c["k"][b] for b in burnup],
                "pcm": [
                    1e5 * (1.0 / reference["k"][b] - 1.0 / c["k"][b]) for b in burnup
                ],
                "seconds": c["seconds"],
                "ratio": [s / r for s, r in zip(c["seconds"], reference["seconds"])],
            }
        )
    return {
        "format": "xlc-compare-ufg-subset",
        "fuel": fuel,
        "burnup": burnup,
        "cases": rows,
        "inputs": inputs,
    }


def burnup_header(summaries: list[dict]) -> list[str]:
    burnup = summaries[0]["burnup"]
    if any(d["burnup"] != burnup for d in summaries):
        raise SystemExit("the fuels have other burnup points")
    return [f"{b:g}" for b in burnup]


def k_table(summaries: list[dict]) -> str:
    """k∞ の差（pcm）の表。行は燃料と組（基準の `all` を除く）、列は燃焼度 (MWd/kg)"""
    rows = [
        [f"{d['fuel'].upper()} {c['case']}", *(signed(x, 1) for x in c["pcm"])]
        for d in summaries
        for c in d["cases"]
        if c["case"] != "all"
    ]
    return markdown_table(["燃料と組", *burnup_header(summaries)], rows)


def time_table(summaries: list[dict]) -> str:
    """自己遮蔽の時間 (s) と、括弧の中に `all` に対する比の表"""
    rows = [
        [
            f"{d['fuel'].upper()} {c['case']}",
            *(f"{s:.2f}（{r:.2f}）" for s, r in zip(c["seconds"], c["ratio"])),
        ]
        for d in summaries
        for c in d["cases"]
    ]
    return markdown_table(["燃料と組", *burnup_header(summaries)], rows)


def report(paths: list[Path]) -> None:
    summaries = [json.loads(path.read_text()) for path in paths]
    print(k_table(summaries))
    print()
    print(time_table(summaries))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--report", nargs="+", type=Path)
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("inputs")
    c = sub.add_parser("compare")
    c.add_argument("--json", type=Path, required=True)
    c.add_argument("fuel", choices=FUELS)
    args = parser.parse_args()
    if args.report:
        report(args.report)
    elif args.command == "inputs":
        write_inputs()
    elif args.command == "compare":
        args.json.write_text(json.dumps(compare(args.fuel, HERE), indent=1) + "\n")
    else:
        parser.error("give inputs, compare or --report")


if __name__ == "__main__":
    main()
