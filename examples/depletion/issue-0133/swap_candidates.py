"""#0133: 燃焼で増える k の差の候補ごとに、断面積を OpenMC の値に置き換えて測る k への寄与をまとめる。

使い方:
    python3 swap_candidates.py lists <格子の入力 (TOML)>       # swap_absorption の例の置き換えの一覧（1 行に 1 つ）
    python3 swap_candidates.py transfers                       # swap_transfer の例の引数（物質、最初の群、最後の群）
    python3 swap_candidates.py summarize <設定 (TOML)> <JSON>  # 置き換えのログから寄与の要約を書く

寄与は置き換えた k の変化の符号を変えたもの（pcm。xlc の断面積のその差が k の差（xlc − OpenMC）に足している分）。吸収は
`swap_absorption`（捕獲 + 核分裂、核分裂、νσf。SPH 因子は置き換える前のまま）、移行は `swap_transfer` の「SPH kept」（同じく
置き換える前の SPH 因子）。OpenMC の実行が複数あれば、実行ごとの寄与を平均する。群 1 に 10 MeV より上を数えない実行（`fast =
false`）は、15 keV より上の群を含む置き換えから外す。

設定は `[[state]]`（`name`、`lattice`（設定のディレクトリから））ごとに `[[state.run]]`（`openmc`（リポジトリの直下から。OpenMC の
k を読む msgpack）、`absorption`、`transfer`（どちらも設定のディレクトリからのログ。`transfer` は無くてよい）、`fast`）を並べる。
"""

import json
import math
import re
import sys
from pathlib import Path

import tomllib

ROOT = Path(__file__).resolve().parents[3]
ACTINIDES = {"Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf"}
# 群は WIMS 69 群の番号（1 始まり、高いエネルギーから）: 1〜13 は 15.03 keV より上、14〜27 は UFG の範囲（4 eV〜15.03 keV）で 27 は
# 4〜9.877 eV、28〜45 は 0.625〜4 eV、46〜56 は 0.1〜0.625 eV、57〜69 は 0.1 eV より下
LABELS = {
    "fast-u238": "U-238、15 keV より上（#0148）",
    "fast-other": "ほかのアクチニド、15 keV より上",
    "ufg27-u238": "U-238、4〜9.88 eV（#0152）",
    "ufg27-other": "U-238 と O-16 を除く燃料の核種、4〜9.88 eV（#0152）",
    "ufg-u238": "U-238、9.88 eV〜15 keV",
    "ufg-other": "U-238 と O-16 を除く燃料の核種、9.88 eV〜15 keV",
    "epithermal": "燃料の全核種、0.625〜4 eV",
    "thermal-u235": "U-235、0.625 eV より下（#0151）",
    "thermal-u238": "U-238、0.625 eV より下（#0151）",
    "thermal-pu": "Pu-239 と Pu-241、0.625 eV より下（#0151）",
    "thermal-other": "ほかの燃料の核種、0.625 eV より下（#0151）",
    "thermal-46-56": "燃料の全核種、0.1〜0.625 eV（#0151）",
    "thermal-57-69": "燃料の全核種、0.1 eV より下（#0151）",
    "fp": "核分裂生成物、4 eV〜15 keV と 4 eV より下（#0133）",
    "xe-sm": "Xe-135 と Sm-149、15 keV より下",
    "all": "吸収の置き換えの全部（燃料の全核種、15 keV より上はアクチニド）",
    "transfer-fuel-fast": "移行: 燃料、15 keV より上（#0148）",
    "transfer-water-fast": "移行: 減速材、15 keV より上（#0148）",
    "transfer-zircaloy-fast": "移行: 被覆管、15 keV より上（#0148）",
    "transfer-fuel-ufg": "移行: 燃料、4 eV〜15 keV",
    "transfer-water-ufg": "移行: 減速材、4 eV〜15 keV（#0152）",
    "transfer-zircaloy-low": "移行: 被覆管、15 keV より下",
    "transfer-fuel-thermal": "移行: 燃料、4 eV より下",
    "transfer-water-epithermal": "移行: 減速材、0.625〜4 eV",
    "transfer-water-thermal": "移行: 減速材、0.625 eV より下（#0151）",
    "rest": "残り（差 − 吸収の全部 − 移行の和。移行のある実行ごと）",
}
TRANSFERS = [
    ("transfer-fuel-fast", "fuel", 1, 13),
    ("transfer-water-fast", "water", 1, 13),
    ("transfer-zircaloy-fast", "zircaloy", 1, 13),
    ("transfer-fuel-ufg", "fuel", 14, 27),
    ("transfer-water-ufg", "water", 14, 27),
    ("transfer-zircaloy-low", "zircaloy", 14, 69),
    ("transfer-fuel-thermal", "fuel", 28, 69),
    ("transfer-water-epithermal", "water", 28, 45),
    ("transfer-water-thermal", "water", 46, 69),
]


class SummaryError(Exception):
    pass


def is_actinide(name: str) -> bool:
    return name.split("-")[0] in ACTINIDES


def absorption_items(nuclides: list[str]) -> list[tuple[str, str]]:
    """候補ごとの (名前, swap_absorption の置き換え)。組成に核種の無い置き換えは書かない"""
    actinides = [n for n in nuclides if is_actinide(n)]
    others = [n for n in nuclides if n not in ("U-238", "O-16")]
    products = [n for n in nuclides if not is_actinide(n) and n != "O-16"]
    heavy = ("U-235", "U-238", "Pu-239", "Pu-241", "O-16")

    def spec(names: list[str], groups: str) -> str:
        return ",".join(f"{n}:{groups}" for n in names)

    items = [
        ("fast-u238", spec(["U-238"], "1-13")),
        ("fast-other", spec([n for n in actinides if n != "U-238"], "1-13")),
        ("ufg27-u238", "U-238:27-27"),
        ("ufg27-other", spec(others, "27-27")),
        ("ufg-u238", "U-238:14-26"),
        ("ufg-other", spec(others, "14-26")),
        ("epithermal", "*:28-45"),
        ("thermal-u235", "U-235:46-69"),
        ("thermal-u238", "U-238:46-69"),
        (
            "thermal-pu",
            spec([n for n in ("Pu-239", "Pu-241") if n in nuclides], "46-69"),
        ),
        ("thermal-other", spec([n for n in nuclides if n not in heavy], "46-69")),
        ("thermal-46-56", "*:46-56"),
        ("thermal-57-69", "*:57-69"),
        ("fp", spec(products, "14-69")),
        ("xe-sm", spec([n for n in ("Xe-135", "Sm-149") if n in nuclides], "14-69")),
        ("all", "*:14-69," + spec(actinides, "1-13")),
    ]
    return [(key, s) for key, s in items if s]


def parse_absorption(text: str) -> tuple[float, dict[str, float]]:
    """swap_absorption のログから、置き換えない k と、置き換えごとの k の変化 (pcm)"""
    lines = [line for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith("k "):
        raise SummaryError("the absorption log does not start with the unswapped k")
    k0 = float(lines[0].split()[1])
    changes = {}
    for line in lines[1:]:
        m = re.fullmatch(r"(.+): k \d+\.\d+ \(([+-]\d+\.\d+) pcm\)", line)
        if m is None:
            raise SummaryError(f"cannot read the absorption log line: {line}")
        if m.group(1) in changes:
            raise SummaryError(f"the absorption log repeats `{m.group(1)}`")
        changes[m.group(1)] = float(m.group(2))
    return k0, changes


def k_of(line: str) -> float:
    m = re.search(r" k (\d+\.\d+)", line)
    if m is None:
        raise SummaryError(f"cannot read k in the transfer log line: {line}")
    return float(m.group(1))


def transfer_runs(text: str) -> dict[tuple[str, int, int], dict[str, float]]:
    """`material <名前>` を前に置いた swap_transfer のログから、(物質, 最初の群, 最後の群) ごとの、OpenMC の k（`openmc`）、置き換えない
    xlc の k（`base`）、置き換えた k の変化（`change`、pcm）"""
    found = {}
    material = None
    block: dict = {}
    for line in text.splitlines():
        if line.startswith("material "):
            material, block = line.split()[1], {}
        elif line.startswith("groups "):
            m = re.match(r"groups (\d+)\.\.=(\d+);", line)
            if m is None:
                raise SummaryError(
                    f"cannot read the groups in the transfer log line: {line}"
                )
            block = {"groups": (int(m.group(1)), int(m.group(2)))}
        elif line.startswith("OpenMC k "):
            block["openmc"] = k_of(line)
        elif line.startswith("xlc "):
            block["base"] = k_of(line)
        elif line.startswith("swapped, SPH kept"):
            if material is None or not {"groups", "openmc", "base"} <= block.keys():
                raise SummaryError(
                    f"a swapped k without its material, groups, OpenMC k or base: {line}"
                )
            key = (material, *block["groups"])
            if key in found:
                raise SummaryError(
                    f"the transfer log repeats {material} {key[1]}-{key[2]}"
                )
            found[key] = {
                "openmc": block["openmc"],
                "base": block["base"],
                "change": (k_of(line) - block["base"]) * 1e5,
            }
            block = {}
    return found


def mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def summarize_state(
    name: str,
    items: list[tuple[str, str]],
    transfers: list[tuple[str, str, int, int]],
    runs: list[dict],
) -> dict:
    """1 つの状態の、OpenMC の実行ごとの寄与と平均。残り（k の差 − 吸収の全部 − 移行の和）は、移行のある実行ごとに、その実行の
    OpenMC の k と置き換えで求めて平均する"""
    if not runs:
        raise SummaryError(f"{name}: no OpenMC runs")
    parsed = [parse_absorption(run["absorption"]) for run in runs]
    k0 = parsed[0][0]
    if any(abs(k - k0) > 5e-7 for k, _ in parsed):
        raise SummaryError(f"{name}: the runs were swapped on other xlc solutions")
    values: dict[str, list[float]] = {}
    for key, spec in items:
        values[key] = []
        for run, (_, changes) in zip(runs, parsed):
            if not run["fast"] and (key.startswith("fast-") or key == "all"):
                continue
            if spec not in changes:
                raise SummaryError(f"{name}: a log has no line `{spec}`")
            values[key].append(-changes[spec])
    for key, *_ in transfers:
        values[key] = []
    rests = []
    all_spec = dict(items).get("all")
    for run, (_, changes) in zip(runs, parsed):
        if run.get("transfer") is None:
            continue
        found = transfer_runs(run["transfer"])
        moved = []
        for key, material, first, last in transfers:
            if (material, first, last) not in found:
                raise SummaryError(
                    f"{name}: a transfer log has no {material} {first}-{last}"
                )
            block = found[(material, first, last)]
            # どちらのログも k を小数 6 桁で書くので、同じ解なら印字の丸めの半分（5e-7）の内で合う
            if abs(block["base"] - k0) > 5e-7:
                raise SummaryError(
                    f"{name}: the transfer log solved another xlc k ({block['base']} and {k0})"
                )
            if abs(block["openmc"] - run["k"]) > 5e-7:
                raise SummaryError(
                    f"{name}: the transfer log was swapped against another OpenMC run "
                    f"(k {block['openmc']} and {run['k']})"
                )
            values[key].append(-block["change"])
            moved.append(-block["change"])
        if all_spec is not None and run["fast"]:
            rests.append((k0 - run["k"]) * 1e5 + changes[all_spec] - sum(moved))
    k_openmc = mean([run["k"] for run in runs])
    out = {key: {"pcm": mean(v), "runs": v} for key, v in values.items()}
    out["rest"] = {"pcm": mean(rests), "runs": rests}
    return {
        "name": name,
        "k_xlc": k0,
        "openmc": [{"k": run["k"], "k_std": run["k_std"]} for run in runs],
        "diff_pcm": (k0 - k_openmc) * 1e5,
        "diff_std_pcm": math.sqrt(sum(run["k_std"] ** 2 for run in runs))
        / len(runs)
        * 1e5,
        "items": out,
    }


def fuel_nuclides(lattice: Path) -> list[str]:
    return list(tomllib.loads(lattice.read_text())["material"]["fuel"]["nuclides"])


def summarize(config: Path) -> dict:
    import msgpack

    base = config.parent
    states = []
    order: list[str] = []
    for state in tomllib.loads(config.read_text())["state"]:
        items = absorption_items(fuel_nuclides(base / state["lattice"]))
        runs = []
        for run in state["run"]:
            reference = msgpack.unpackb((ROOT / run["openmc"]).read_bytes())
            runs.append(
                {
                    "k": reference["k"],
                    "k_std": reference["k_std"],
                    "absorption": (base / run["absorption"]).read_text(),
                    "transfer": (base / run["transfer"]).read_text()
                    if "transfer" in run
                    else None,
                    "fast": run["fast"],
                }
            )
        s = summarize_state(state["name"], items, TRANSFERS, runs)
        s["sources"] = state["run"]
        states.append(s)
        order += [key for key in s["items"] if key not in order]
    order = [key for key in LABELS if key in order]
    return {
        "format": "xlc-swap-candidates",
        "order": order,
        "labels": {key: LABELS[key] for key in order},
        "states": states,
    }


def main(argv: list[str]) -> int:
    match argv:
        case ["lists", lattice]:
            for _, spec in absorption_items(fuel_nuclides(Path(lattice))):
                print(spec)
        case ["transfers"]:
            for _, material, first, last in TRANSFERS:
                print(material, first, last)
        case ["summarize", config, out]:
            try:
                summary = summarize(Path(config))
            except SummaryError as e:
                print(f"FAILED: {e}", file=sys.stderr)
                return 1
            Path(out).write_text(
                json.dumps(summary, ensure_ascii=False, indent=1) + "\n"
            )
        case _:
            print(__doc__, file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
