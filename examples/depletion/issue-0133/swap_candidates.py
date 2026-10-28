"""#0133: 燃焼で増える k の差の候補ごとに、断面積を OpenMC の値に置き換えて測る k への寄与をまとめる。

使い方:
    python3 swap_candidates.py lists <格子の入力 (TOML)>       # swap_absorption の例の燃料の核種の置き換えの一覧（1 行に 1 つ）
    python3 swap_candidates.py materials                       # swap_absorption の例の物質の吸収の置き換えの一覧
    python3 swap_candidates.py joint <格子の入力 (TOML)>       # 同時の置き換えと、それと同じ実行で別々に解く置き換えの一覧
    python3 swap_candidates.py transfers                       # swap_transfer の例の引数（物質、最初の群、最後の群）
    python3 swap_candidates.py summarize <設定 (TOML)> <JSON>  # 置き換えのログから寄与の要約を書く

寄与は置き換えた k の変化の符号を変えたもの（pcm。xlc の断面積のその差が k の差（xlc − OpenMC）に足している分）。燃料の核種の
吸収は `swap_absorption`（捕獲 + 核分裂、核分裂、νσf。SPH 因子は置き換える前のまま）、物質の吸収と χ も `swap_absorption`
（`--materials`、`--chi`）、移行は `swap_transfer` の「SPH kept」（同じく置き換える前の SPH 因子）。OpenMC の実行が複数あれば、
実行ごとの寄与を平均する。群 1 に 10 MeV より上を数えない実行（`fast = false`）は、群 1 を含む置き換えから外す。

置き換えどうしの相互作用は、吸収の全部（燃料の全核種）、物質の吸収、移行を 1 つの実行で同時に置き換えた寄与から、同じ実行で
別々に置き換えた寄与の和を引いたもの。残りは、同時の置き換えの OpenMC の実行の k の差から、同時の置き換えと χ の寄与を引いたもの
（差 − 吸収の全部 − 物質の吸収 − 移行の和 − χ − 相互作用と同じ）。

設定は `[[state]]`（`name`、`lattice`（設定のディレクトリから））ごとに `[[state.run]]`（`openmc`（リポジトリの直下から。OpenMC の
k を読む msgpack）、`absorption`、`materials`、`transfer`（どれも設定のディレクトリからのログ。`materials` と `transfer` は無くて
よい）、`fast`）を並べ、`[state.joint]`（`openmc`（燃料の吸収と移行の実行）、`materials`（物質の吸収の実行。どちらも `[[state.run]]`
の `openmc` のどれか）、`log`）と `[state.chi]`（`openmc`（openmc_chi.py の出力。記録のため）、`log`）を書く。
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
    "material-water-fast": "吸収: 減速材、15 keV より上",
    "material-water-ufg": "吸収: 減速材、4 eV〜15 keV",
    "material-water-epithermal": "吸収: 減速材、0.625〜4 eV",
    "material-water-thermal": "吸収: 減速材、0.625 eV より下",
    "material-zircaloy-fast": "吸収: 被覆管、15 keV より上",
    "material-zircaloy-ufg": "吸収: 被覆管、4 eV〜15 keV",
    "material-zircaloy-epithermal": "吸収: 被覆管、0.625〜4 eV",
    "material-zircaloy-thermal": "吸収: 被覆管、0.625 eV より下",
    "transfer-fuel-fast": "移行: 燃料、15 keV より上（#0148）",
    "transfer-water-fast": "移行: 減速材、15 keV より上（#0148）",
    "transfer-zircaloy-fast": "移行: 被覆管、15 keV より上（#0148）",
    "transfer-fuel-ufg": "移行: 燃料、4 eV〜15 keV",
    "transfer-water-ufg": "移行: 減速材、4 eV〜15 keV（#0152）",
    "transfer-zircaloy-low": "移行: 被覆管、15 keV より下",
    "transfer-fuel-thermal": "移行: 燃料、4 eV より下",
    "transfer-water-epithermal": "移行: 減速材、0.625〜4 eV",
    "transfer-water-thermal": "移行: 減速材、0.625 eV より下（#0151）",
    "chi": "χ（燃料のアクチニド）",
    "interaction": "置き換えどうしの相互作用（吸収の全部、物質の吸収、移行の同時 − 別々の和）",
    "rest": "残り（差 − 吸収の全部 − 物質の吸収 − 移行の和 − χ − 相互作用）",
}
# 物質の吸収を置き換える物質（xlc の入力の名前、openmc_pin_rates.py の OpenMC の物質の名前）とエネルギーの範囲
MATERIALS = [("water", "moderator"), ("zircaloy", "zircaloy-4")]
RANGES = [("fast", 1, 13), ("ufg", 14, 27), ("epithermal", 28, 45), ("thermal", 46, 69)]
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


def material_items() -> list[tuple[str, str]]:
    """物質と範囲ごとの (名前, swap_absorption の物質の吸収の置き換え)"""
    return [
        (f"material-{xlc}-{key}", f"@{xlc}={openmc}:{first}-{last}")
        for xlc, openmc in MATERIALS
        for key, first, last in RANGES
    ]


def joint_lines(
    items: list[tuple[str, str]],
    materials: list[tuple[str, str]],
    transfers: list[tuple[str, str, int, int]],
) -> list[str]:
    """同時の置き換えのログの一覧: 別々の置き換え（吸収の全部、物質の吸収、移行）と、それらをつないだ同時の置き換え（最後の行）"""
    all_spec = dict(items).get("all")
    if all_spec is None:
        raise SummaryError("no line swaps all the absorption")
    separate = [
        all_spec,
        *(spec for _, spec in materials),
        *(f"~{material}:{first}-{last}" for _, material, first, last in transfers),
    ]
    return [*separate, ",".join(separate)]


def touches_group_one(spec: str) -> bool:
    """置き換えが群 1（10 MeV より上を数えるかで変わる群）を含むか"""
    return any(part.split(":")[-1].startswith("1-") for part in spec.split(","))


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
    materials: list[tuple[str, str]],
    transfers: list[tuple[str, str, int, int]],
    runs: list[dict],
    joint: dict,
    chi: dict,
) -> dict:
    """1 つの状態の、OpenMC の実行ごとの寄与と平均。`joint` は同時の置き換えのログ（`log`）と、その燃料の吸収と移行の実行（`run`）、
    物質の吸収の実行（`materials_run`）の番号。`chi` は χ の置き換えのログ（`log`）"""
    if not runs:
        raise SummaryError(f"{name}: no OpenMC runs")
    parsed = [parse_absorption(run["absorption"]) for run in runs]
    k0 = parsed[0][0]
    if any(abs(k - k0) > 5e-7 for k, _ in parsed):
        raise SummaryError(f"{name}: the runs were swapped on other xlc solutions")

    def same_xlc(k: float, what: str) -> None:
        if abs(k - k0) > 5e-7:
            raise SummaryError(
                f"{name}: the {what} log solved another xlc k ({k} and {k0})"
            )

    values: dict[str, list[float]] = {}
    for key, spec in items:
        values[key] = []
        for run, (_, changes) in zip(runs, parsed):
            if not run["fast"] and touches_group_one(spec):
                continue
            if spec not in changes:
                raise SummaryError(f"{name}: a log has no line `{spec}`")
            values[key].append(-changes[spec])
    material_changes = {}
    for key, _ in materials:
        values[key] = []
    for index, run in enumerate(runs):
        if run.get("materials") is None:
            continue
        k, changes = parse_absorption(run["materials"])
        same_xlc(k, "materials")
        material_changes[index] = changes
        for key, spec in materials:
            if spec not in changes:
                raise SummaryError(f"{name}: a materials log has no line `{spec}`")
            if run["fast"] or not touches_group_one(spec):
                values[key].append(-changes[spec])
    for key, *_ in transfers:
        values[key] = []
    transfer_changes = {}
    for index, run in enumerate(runs):
        if run.get("transfer") is None:
            continue
        found = transfer_runs(run["transfer"])
        transfer_changes[index] = found
        for key, material, first, last in transfers:
            if (material, first, last) not in found:
                raise SummaryError(
                    f"{name}: a transfer log has no {material} {first}-{last}"
                )
            block = found[(material, first, last)]
            # 別々の実行が小数 6 桁で書いた 2 つの k なので、同じ解でも丸めの境をまたぐと最後の桁が 1 違いうる
            if abs(block["base"] - k0) > 1e-6 + 1e-12:
                raise SummaryError(
                    f"{name}: the transfer log solved another xlc k ({block['base']} and {k0})"
                )
            if abs(block["openmc"] - run["k"]) > 5e-7:
                raise SummaryError(
                    f"{name}: the transfer log was swapped against another OpenMC run "
                    f"(k {block['openmc']} and {run['k']})"
                )
            values[key].append(-block["change"])

    # 同時の置き換え。別々の行は、同じ OpenMC の値で解いたほかのログの行と合うこと（同じ置き換えなので、ログの丸めの中で同じ）
    i, j = joint["run"], joint["materials_run"]
    if i not in transfer_changes or not runs[i]["fast"]:
        raise SummaryError(
            f"{name}: the joint swap needs a run with transfers above 10 MeV"
        )
    if j not in material_changes or not runs[j]["fast"]:
        raise SummaryError(
            f"{name}: the joint swap needs a run with materials above 10 MeV"
        )
    k, joined = parse_absorption(joint["log"])
    same_xlc(k, "joint")
    lines = joint_lines(items, materials, transfers)
    missing = [line for line in lines if line not in joined]
    if missing:
        raise SummaryError(f"{name}: the joint log has no line `{missing[0]}`")
    all_spec, together = lines[0], lines[-1]
    expected = {all_spec: (parsed[i][1][all_spec], 0.1)}
    for _, spec in materials:
        expected[spec] = (material_changes[j][spec], 0.1)
    for _, material, first, last in transfers:
        # swap_transfer の変化は小数 6 桁の 2 つの k の差なので、丸めの分だけ広く
        change = transfer_changes[i][(material, first, last)]["change"]
        expected[f"~{material}:{first}-{last}"] = (change, 0.2)
    for spec, (change, tolerance) in expected.items():
        if abs(joined[spec] - change) > tolerance + 1e-9:
            raise SummaryError(
                f"{name}: the joint log has `{spec}` {joined[spec]:+.1f} pcm, the separate log {change:+.1f}"
            )
    joint_pcm = -joined[together]
    interaction = joint_pcm + sum(joined[spec] for spec in lines[:-1])
    k, chi_changes = parse_absorption(chi["log"])
    same_xlc(k, "chi")
    if "chi" not in chi_changes:
        raise SummaryError(f"{name}: the chi log has no line `chi`")
    chi_pcm = -chi_changes["chi"]
    rest = (k0 - runs[i]["k"]) * 1e5 - joint_pcm - chi_pcm

    k_openmc = mean([run["k"] for run in runs])
    out = {key: {"pcm": mean(v), "runs": v} for key, v in values.items()}
    out["chi"] = {"pcm": chi_pcm, "runs": [chi_pcm]}
    out["interaction"] = {"pcm": interaction, "runs": [interaction]}
    out["rest"] = {"pcm": rest, "runs": [rest]}
    return {
        "name": name,
        "k_xlc": k0,
        "openmc": [{"k": run["k"], "k_std": run["k_std"]} for run in runs],
        "diff_pcm": (k0 - k_openmc) * 1e5,
        "diff_std_pcm": math.sqrt(sum(run["k_std"] ** 2 for run in runs))
        / len(runs)
        * 1e5,
        "joint_pcm": joint_pcm,
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
                    **{
                        key: (base / run[key]).read_text() if key in run else None
                        for key in ("materials", "transfer")
                    },
                    "fast": run["fast"],
                }
            )
        sources = [run["openmc"] for run in state["run"]]
        joint = state["joint"]
        for key in ("openmc", "materials"):
            if sources.count(joint[key]) != 1:
                raise SummaryError(
                    f"{state['name']}: the joint {key} run {joint[key]} is not one of the runs"
                )
        s = summarize_state(
            state["name"],
            items,
            material_items(),
            TRANSFERS,
            runs,
            {
                "run": sources.index(joint["openmc"]),
                "materials_run": sources.index(joint["materials"]),
                "log": (base / joint["log"]).read_text(),
            },
            {"log": (base / state["chi"]["log"]).read_text()},
        )
        s["sources"] = state["run"]
        s["joint_sources"] = joint
        s["chi_sources"] = state["chi"]
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
        case ["materials"]:
            for _, spec in material_items():
                print(spec)
        case ["joint", lattice]:
            items = absorption_items(fuel_nuclides(Path(lattice)))
            for line in joint_lines(items, material_items(), TRANSFERS):
                print(line)
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
