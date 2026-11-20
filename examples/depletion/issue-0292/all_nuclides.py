"""#0292: データのある全核種を輸送する試行のライブラリの入力と、燃焼の入力を、核データの版（`VERSIONS`）ごとに書く。

版の燃焼チェーン（#0303 の版ごとのチェーン。ADR-0263）にあって、その版の中性子のファイルがある核種（読み込みの走査
../read-survey-endfb<版>.toml のチェーンの評価。天然の元素の評価は、チェーンの入力の `natural-isotopes` の同じ元素の核種それぞれに
当てる。ADR-0262）を、すべて輸送する核種にする。核種の並びと S(α, β) の核種の項目は 154 核種の試行のライブラリの入力
（trial/library/xs.toml。VIII.0）から取り、試行の入力の核種のうち版に評価の無いものは入れない（`missing_from_version`）。ファイルは版の
評価、温度は試行のセルの温度 `TEMPERATURE` の 1 点にする（S(α, β) の核種は内挿に使う温度で、版の熱散乱のファイルに無い温度を
除いた `thermal_temperatures`）。評価の欠陥を読み替える宣言は、版の走査の宣言（process-survey/<版>/declarations.toml）から足す
（`declare`）。走査で止まる評価は、`corrections` の宣言で読める（ADR-0272）もののほかは `EXCLUDED`（核種と issue）で、止まる評価と
食い違えば書かずに止める。外した核種のうち試行の入力で放射化の核種だったものは、放射化の核種のまま残す（冷却計算の断面積を
失わない。ADR-0238）。

- `python3 all_nuclides.py library [--version <版>]` は <置き場>/library/xs.toml（全核種）、ジョブが `xs nuclide` で処理する核種の一覧
  （<置き場>/library/nuclides.txt）、p2 のライブラリが輸送しない核種の一覧（<置き場>/added-to-p2.txt）を書く。VIII.0（既定）は
  さらに all/control/xs.toml（同じ核種の項目で、試行のライブラリが輸送する核種だけを輸送し、ほかを放射化の核種にした対照）と
  その nuclides.txt、all/added-to-trial.txt、all/unranked.txt を書く
- `python3 all_nuclides.py depletion [--version <版>]` は <置き場>/<燃料>/nominal.toml（trial/<燃料>/nominal.toml のライブラリを
  <置き場>/library に、チェーンを版のチェーンに替えたもの）を書く。VIII.0 はさらに確認のセル（inputs.py の `CHECK_FUELS`）の
  all/<セル>/nominal.toml（`check_depletion_input`）を書く
"""

import argparse
import sys
from pathlib import Path

import msgpack
import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import inputs
from read_survey import nuclide_key

FUELS = ("u", "mox", "gd")
# 試行のセル（inputs.py）の燃料、被覆管、減速材の温度 (K)。lattice は S(α, β) の核種のほかは表の温度と一致する温度だけを
# 引く（crates/format/src/library.rs の `lookup`）ので、この 1 点で試行のセルの計算は 13 温度のライブラリと同じになる
TEMPERATURE = 565
# 走査で xs が止まる、チェーンの評価（核種 → (issue, 理由)）。どの版も、宣言の後は走査の全評価が `xs nuclide` を通る
# （examples/depletion/process-survey.txt）
EXCLUDED: dict[str, tuple[str, str]] = {}
TRIAL = HERE / "trial" / "library" / "xs.toml"


def version_config(
    key: str, name: str, chain_place: str, thermal_temperatures: list | None
) -> dict:
    """版の設定: 置き場（issue-0292 から）、走査、チェーン（msgpack と入力）、評価の欠陥を読み替える核種の宣言（版の走査の宣言と
    同じファイル。process_survey.py の `declare`）、S(α, β) の温度（None は試行の入力のまま）。`trial_version` は試行の入力
    （trial/library/xs.toml）と同じ版（VIII.0）で、試行の核種はどれもその版のチェーンの評価のはずなので、食い違えば止める。対照、
    試行と比べる一覧、確認のセルも、この版だけが持つ"""
    chain = HERE.parent / chain_place
    trial_version = key == "endfb80"
    return {
        "key": key,
        "trial_version": trial_version,
        "directory": "all" if trial_version else f"all-{key}",
        "label": f"ENDF/B-{name}",
        "neutrons": f"ENDF-B-{name}",
        "survey": HERE.parent / f"read-survey-{key}.toml",
        "chain": chain / "chain.msgpack",
        "chain_input": chain / "chain.toml",
        "declarations": HERE.parent / "process-survey" / name / "declarations.toml",
        "thermal_temperatures": thermal_temperatures,
    }


# 版の置き場と名前は #0303 の版ごとの入力（../endfb70 など。VIII.0 は ../ の直下）と同じ。VII.0 と VII.1 の tsl-HinH2O.endf には
# 試行の 4 点のうち 573.6 K と 623.6 K が無い（docs/environment.md）ので、#0303 の版ごとのライブラリと同じく 550 K と 600 K にする
VERSIONS = {
    "endfb70": version_config("endfb70", "VII.0", "endfb70", [550, 600]),
    "endfb71": version_config("endfb71", "VII.1", "endfb71", [550, 600]),
    "endfb80": version_config("endfb80", "VIII.0", ".", None),
    "endfb81": version_config("endfb81", "VIII.1", "endfb81", None),
}
# 放射化の核種は遮蔽しないので取らない宣言（xs は書くと止まる。ADR-0274）
SHIELDING_ONLY = ("unresolved-negative-background",)
P2 = HERE.parent / "p2" / "xs.toml"


def chain_evaluations(
    survey: dict, chain: dict, natural: list[str] = ()
) -> dict[str, tuple[str, str | None]]:
    """走査の評価のうちチェーンの核種のもの: 核種の名前 → (ファイルの名前, 止まった理由か None)。走査の順。天然の元素の評価
    （A = 0）は、`natural`（チェーンの入力の `natural-isotopes`）の同じ元素の核種それぞれに、その並びで当てる（ADR-0262。チェーンは
    その核種に同じ評価の反応を持たせる）"""
    if chain.get("format") != "xlc-chain":
        raise SystemExit("not an xlc chain")
    names = {(n["z"], n["a"], n["state"]): n["name"] for n in chain["nuclides"]}
    z_of = {n["name"]: n["z"] for n in chain["nuclides"]}
    out, given = {}, set()
    for e in survey["evaluation"]:
        z, a, state = nuclide_key(e["endf"])
        if z > 0 and a == 0:
            isotopes = [n for n in natural if z_of.get(n) == z]
            given.update(isotopes)
            if not isotopes:
                raise SystemExit(
                    f"{e['endf']}: no natural isotope of Z = {z} is listed"
                )
        else:
            isotopes = [names[(z, a, state)]] if (z, a, state) in names else []
        for name in isotopes:
            if name in out:
                raise SystemExit(f"{name} has two evaluations in the survey")
            out[name] = (Path(e["endf"]).name, e.get("error"))
    ungiven = [n for n in natural if n not in given]
    if ungiven:
        raise SystemExit(
            f"{ungiven} are listed but no natural evaluation of their element is in the survey"
        )
    return out


def library_entries(
    base: dict,
    survey: dict,
    chain: dict,
    excluded: dict,
    natural: list[str] = (),
    resolved=frozenset(),
    thermal_temperatures: list | None = None,
    trial_version: bool = True,
) -> dict[str, dict]:
    """全核種のライブラリの核種の項目（試行の入力 `base` の核種のうち評価のあるものを `base` の順、続いて `base` に無い評価を
    走査の順）。ファイルは版の評価のもの。`resolved` は、走査で止まるが宣言（`corrections`）で読める核種。`thermal_temperatures` は
    S(α, β) の核種の温度（None は `base` のまま）。`trial_version` は試行の入力と同じ版で、`base` の核種のどれかがチェーンの評価で
    ないか、ファイルが評価と違えば止める（ほかの版は、評価の無い核種を入れず、ファイルを版のものにする）"""
    evaluations = chain_evaluations(survey, chain, natural)
    stopped = {
        n
        for n, (_, error) in evaluations.items()
        if error is not None and n not in resolved
    }
    if stopped != set(excluded):
        raise SystemExit(
            f"the exclusions differ from the evaluations xs stops on in {sorted(stopped ^ set(excluded))}"
        )
    if trial_version:
        unknown = [
            n for n, v in base.items() if evaluations.get(n, (None,))[0] != v["endf"]
        ]
        if unknown:
            raise SystemExit(f"{unknown} of the base input are not chain evaluations")
    out = {}
    for name, entry in base.items():
        if name not in evaluations:
            continue
        if name in excluded and not entry.get("activation"):
            raise SystemExit(f"{name} is excluded but the base input transports it")
        entry = {**entry, "endf": evaluations[name][0]}
        if "thermal" in entry:
            if thermal_temperatures is not None:
                entry["temperatures"] = list(thermal_temperatures)
        else:
            entry["temperatures"] = [TEMPERATURE]
        if name not in excluded:
            entry = {k: v for k, v in entry.items() if k != "activation"}
        out[name] = entry
    for name, (endf, _) in evaluations.items():
        if name not in out and name not in excluded:
            out[name] = {"endf": endf, "temperatures": [TEMPERATURE]}
    return out


def missing_from_version(base: dict, entries: dict) -> list[str]:
    """試行の入力 `base` の核種のうち、版の全核種のライブラリに無いもの（版に評価が無いか、外した核種）"""
    return [n for n in base if n not in entries]


def transported(entries: dict) -> list[str]:
    return [n for n, v in entries.items() if not v.get("activation")]


def declare(entries: dict, declarations: dict) -> dict[str, dict]:
    """`entries` の核種に宣言の項目を足す（放射化の核種には `SHIELDING_ONLY` を足さない）。入力に無い核種への宣言と、試行の入力の
    項目と値の違う宣言は止める（同じ値なら、試行の入力の宣言をそのまま使う）"""
    unknown = sorted(set(declarations) - set(entries))
    if unknown:
        raise SystemExit(f"declarations for nuclides outside the entries: {unknown}")
    out = {}
    for name, entry in entries.items():
        items = {
            k: v
            for k, v in declarations.get(name, {}).items()
            if not (entry.get("activation") and k in SHIELDING_ONLY)
        }
        conflicting = sorted(
            k for k, v in items.items() if k in entry and entry[k] != v
        )
        if conflicting:
            raise SystemExit(
                f"{name}: a declaration differs from the base entry in {conflicting}"
            )
        out[name] = {**entry, **items}
    return out


def control_entries(base: dict, entries: dict) -> dict[str, dict]:
    """`entries` の核種の項目のうち、試行の入力 `base` が輸送しない核種を放射化の核種にしたもの（放射化の核種は
    `SHIELDING_ONLY` の宣言を持たない）"""
    keep = set(transported(base))
    return {
        n: v
        if n in keep
        else {
            **{k: x for k, x in v.items() if k not in SHIELDING_ONLY},
            "activation": True,
        }
        for n, v in entries.items()
    }


def nuclide_lists(entries: dict, trial: dict, p2: dict) -> dict[str, list[str]]:
    """輸送する核種のうち、p2 のライブラリが輸送しないもの、試行のライブラリが輸送しないもの、p2 のライブラリに放射化の
    核種としても無いもの（89 核種の履歴の吸収の割合と冷却計算の断面積が無く、自己遮蔽の目安で順位を付けられない）"""
    names = transported(entries)
    p2_transported, trial_transported = set(transported(p2)), set(transported(trial))
    return {
        "added-to-p2": [n for n in names if n not in p2_transported],
        "added-to-trial": [n for n in names if n not in trial_transported],
        "unranked": [n for n in names if n not in p2],
    }


def value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    if isinstance(v, list):
        return "[" + ", ".join(value(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{ " + ", ".join(f"{k} = {value(x)}" for k, x in v.items()) + " }"
    return str(v)


def render(base_text: str, head: list[str], entries: dict) -> str:
    """注釈 `head`、試行の入力の全核種に共通の項目、核種の項目（[nuclide] のインラインの表）の入力"""
    common = []
    for line in base_text.splitlines():
        if line.startswith("["):
            break
        # ファイルの先頭の注釈（試行の入力の説明）は `head` に替え、共通の項目の間の注釈（書いておく理由など）は残す
        if line and (common or not line.startswith("#")):
            common.append(line)
    rows = [
        f"{n} = {{ " + ", ".join(f"{k} = {value(v)}" for k, v in e.items()) + " }"
        for n, e in entries.items()
    ]
    return "\n".join([*(f"# {h}" for h in head), *common, "", "[nuclide]", *rows, ""])


def exclusion_notes(base: dict, excluded: dict) -> list[str]:
    """外した核種ごとの注釈（試行の入力 `base` で放射化の核種なら、残すことを書く。`library_entries` と同じ判定）。外した核種が
    無ければ「なし」"""
    if not excluded:
        return ["なし"]
    return [
        f"{n}: {reason}（{issue}）"
        + ("。放射化の核種として残す" if base.get(n, {}).get("activation") else "")
        for n, (issue, reason) in excluded.items()
    ]


def job_lists(entries: dict, control: dict) -> dict[str, list[str]]:
    """ジョブが `xs nuclide` に渡す核種: 全核種のライブラリはすべて、対照は放射化の核種だけ（輸送する核種は、範囲とビルドの
    同じ全核種のライブラリの xs のファイルを使う）"""
    return {
        "library": list(entries),
        "control": [n for n, v in control.items() if v.get("activation")],
    }


CONTROL_HEAD = [
    "#0292 の試行の対照: ../library/xs.toml と同じ核種の項目で、154 核種の試行のライブラリ（../../trial/library）が",
    "輸送する核種だけを輸送し、ほかを放射化の核種にしたもの。範囲（xs-range）とビルドを全核種のライブラリとそろえ、",
    "全核種を輸送した差から、範囲と xs の版の違いの分を分ける。all_nuclides.py が書く。手で直さない。",
    "nuclides.txt の放射化の核種だけを ../library/xs-range.toml の範囲で xs-nuclides に作り、組み立ては、../library/xs-nuclides の上に",
    "xs-nuclides を重ねた置き場で ../library/xs-range.toml と行う（輸送する核種は全核種のライブラリの xs のファイル）。",
]


def library_head(
    version: dict, notes: list[str], extra: list[str], entries: dict
) -> list[str]:
    """全核種のライブラリの入力の注釈。`notes` は外した核種、`extra` は版のチェーンと評価の注記（VIII.0 は無い）。S(α, β) の点の数は
    書く核種の項目 `entries` から数える"""
    (thermal,) = [e for e in entries.values() if "thermal" in e]
    points = len(thermal["temperatures"])
    return [
        f"#0292 の試行: 燃焼チェーンにあって {version['label']} の中性子のファイルがある核種を、すべて輸送する核種にしたライブラリ。",
        f"all_nuclides.py が ../../trial/library/xs.toml と ../../../read-survey-{version['key']}.toml から書く。手で直さない。",
        f"温度は試行のセルの {TEMPERATURE} K の 1 点（H-1 は S(α, β) の内挿に使う {points} 点）。理由は all_nuclides.py の TEMPERATURE。",
        *extra,
        "xs が止まる評価は外す（外した核種と issue）:",
        *notes,
    ]


def version_entries(version: dict, base: dict) -> dict[str, dict]:
    """版の全核種のライブラリの核種の項目（宣言を足したもの）"""
    declarations = tomllib.loads(version["declarations"].read_text())["nuclide"]
    natural = tomllib.loads(version["chain_input"].read_text()).get(
        "natural-isotopes", []
    )
    return declare(
        library_entries(
            base,
            tomllib.loads(version["survey"].read_text()),
            msgpack.unpackb(version["chain"].read_bytes(), strict_map_key=False),
            EXCLUDED,
            natural=natural,
            resolved={n for n, d in declarations.items() if "corrections" in d},
            thermal_temperatures=version["thermal_temperatures"],
            trial_version=version["trial_version"],
        ),
        declarations,
    )


def directories(base_text: str, version: dict) -> str:
    """試行の入力の ENDF と熱散乱のディレクトリを版のものに替えた文"""
    old = "ENDF-B-VIII.0"
    lines = [x for x in base_text.splitlines() if "-directory = " in x]
    if len(lines) != 2 or not all(old in x for x in lines):
        raise SystemExit("the trial input has other ENDF directories")
    out = base_text
    for x in lines:
        out = out.replace(x, x.replace(old, version["neutrons"]))
    return out


def library_text(version: dict) -> tuple[str, dict]:
    """版の全核種のライブラリの入力の文と、核種の項目"""
    base_text = TRIAL.read_text()
    base = tomllib.loads(base_text)["nuclide"]
    entries = version_entries(version, base)
    extra = []
    if not version["trial_version"]:
        chain = Path("../../..") / version["chain_input"].relative_to(HERE.parent)
        extra.append(
            f"輸送する核種の組は {version['label']} の燃焼チェーン（{chain}）と評価で決める（#0303 の版ごとのチェーン）。"
        )
        missing = missing_from_version(base, entries)
        if missing:
            extra.append(
                f"試行の入力の核種のうち、この版に評価の無いもの（{len(missing)}）: "
                + ", ".join(missing)
            )
        natural = sorted(
            {e["endf"] for e in entries.values() if nuclide_key(e["endf"])[1] == 0}
        )
        for f in natural:
            names = [n for n, e in entries.items() if e["endf"] == f]
            extra.append(
                f"天然の元素の評価 {f} を当てる核種（ADR-0262）: " + ", ".join(names)
            )
    notes = [f"  {n}" for n in exclusion_notes(base, EXCLUDED)]
    text = render(
        directories(base_text, version),
        library_head(version, notes, extra, entries),
        entries,
    )
    return text, entries


def write_library(version: dict) -> None:
    text, entries = library_text(version)
    root = HERE / version["directory"]
    base = tomllib.loads(TRIAL.read_text())["nuclide"]
    p2 = tomllib.loads(P2.read_text())["nuclide"]
    lists = nuclide_lists(entries, base, p2)
    path = root / "library" / "xs.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    (path.parent / "nuclides.txt").write_text("\n".join(entries) + "\n")
    if not version["trial_version"]:
        (root / "added-to-p2.txt").write_text("\n".join(lists["added-to-p2"]) + "\n")
        return
    # 対照と、試行と比べる一覧は VIII.0 だけ（試行のライブラリが VIII.0）
    control = control_entries(base, entries)
    path = root / "control" / "xs.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render(TRIAL.read_text(), CONTROL_HEAD, control))
    (path.parent / "nuclides.txt").write_text(
        "\n".join(job_lists(entries, control)["control"]) + "\n"
    )
    for name, names in lists.items():
        (root / f"{name}.txt").write_text("\n".join(names) + "\n")


def depletion_input(text: str, version: dict) -> str:
    """試行の燃焼の入力 `text` のライブラリを版の置き場の library に、チェーンを版のチェーンに替えたもの"""
    old = 'library = "../library/depletion.msgpack"'
    lines = text.splitlines()
    if not lines[0].startswith("#") or lines[1] != old:
        raise SystemExit("the trial input has another head or library")
    if version["trial_version"]:
        head = "# #0292 の試行: ../../trial の同じ燃料の nominal.toml の、ライブラリだけを全核種を輸送する ../library に替えたもの。all_nuclides.py が書く。手で直さない"
        return "\n".join([head, old, *lines[2:]]) + "\n"
    chain_old = 'chain = "../../../chain.msgpack"'
    if chain_old not in lines:
        raise SystemExit("the trial input reads another chain")
    chain = Path("../../..") / version["chain"].relative_to(HERE.parent)
    head = (
        f"# #0292 の試行: ../../trial の同じ燃料の nominal.toml の、ライブラリを {version['label']} の全核種を輸送する ../library に、"
        "チェーンを同じ版のチェーンに替えたもの。all_nuclides.py が書く。手で直さない"
    )
    body = [f'chain = "{chain}"' if x == chain_old else x for x in lines[2:]]
    return "\n".join([head, old, *body]) + "\n"


def check_cell_nuclides(text: str, entries: dict) -> None:
    """燃焼の入力 `text` の物質の核種が、どれもライブラリの核種 `entries` にあること（無い核種は lattice が読めない）"""
    materials = tomllib.loads(text).get("material", {})
    missing = sorted(
        {
            n
            for m in materials.values()
            for n in m.get("nuclides", {})
            if n not in entries
        }
    )
    if missing:
        raise SystemExit(f"the cell has nuclides outside the library: {missing}")


def check_depletion_input(fuel: str) -> str:
    """確認のセル（inputs.py の `CHECK_FUELS`。試行の入力を持たない）の燃焼の入力: inputs.py の基準の密度の燃焼の入力を、
    試行の入力と同じ置き場の形（ライブラリは all/library、チェーンは 3 つ上、出力は同じ置き場）にしたもの"""
    lines = inputs.depletion_input(fuel, inputs.NOMINAL_DENSITY, "nominal").splitlines()
    places = {
        'library = "../../p2/depletion.msgpack"': 'library = "../library/depletion.msgpack"',
        'chain = "../../chain.msgpack"': 'chain = "../../../chain.msgpack"',
        'output = "out/nominal.msgpack"': 'output = "nominal.msgpack"',
    }
    missing = [old for old in places if old not in lines]
    if missing:
        raise SystemExit(f"the depletion input of {fuel} has no {missing}")
    head = f"# #0292 の確認のセル: {fuel} の燃焼（inputs.py の形状と組成）を全核種を輸送する ../library で。all_nuclides.py が書く。手で直さない"
    return "\n".join([head, *(places.get(x, x) for x in lines[1:])]) + "\n"


def write_depletion(version: dict) -> None:
    root = HERE / version["directory"]
    entries = tomllib.loads((root / "library" / "xs.toml").read_text())["nuclide"]
    for fuel in FUELS:
        text = depletion_input(
            (HERE / "trial" / fuel / "nominal.toml").read_text(), version
        )
        check_cell_nuclides(text, entries)
        directory = root / fuel
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "nominal.toml").write_text(text)
    if not version["trial_version"]:
        return
    for fuel in inputs.CHECK_FUELS:
        directory = root / fuel
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "nominal.toml").write_text(check_depletion_input(fuel))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["library", "depletion"])
    parser.add_argument("--version", choices=list(VERSIONS), default="endfb80")
    args = parser.parse_args()
    version = VERSIONS[args.version]
    if args.command == "library":
        write_library(version)
    else:
        write_depletion(version)


if __name__ == "__main__":
    main()
