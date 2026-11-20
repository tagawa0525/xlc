"""#0292: データのある全核種を輸送する試行のライブラリの入力と、燃焼の入力を書く。

燃焼チェーン（../chain.msgpack）にあって ENDF/B-VIII.0 の中性子のファイルがある核種（読み込みの走査 ../read-survey-endfb80.toml の
チェーンの評価）を、すべて輸送する核種にする。核種の項目は 154 核種の試行のライブラリの入力（trial/library/xs.toml）から取り、
温度だけを試行のセルの温度 `TEMPERATURE` の 1 点にする（S(α, β) の核種は、内挿に使う温度をそのまま持つ）。走査で止まる評価は
`EXCLUDED`（核種と issue）で、止まる評価と食い違えば書かずに止める。外した核種のうち試行の入力で放射化の核種だったものは、
放射化の核種のまま残す（冷却計算の断面積を失わない。ADR-0238）。

- `python3 all_nuclides.py library` は all/library/xs.toml（全核種）、all/control/xs.toml（同じ核種の項目で、試行のライブラリが
  輸送する核種だけを輸送し、ほかを放射化の核種にした対照）、核種の一覧（all/added-to-p2.txt、all/added-to-trial.txt、
  all/unranked.txt）と、ジョブが `xs nuclide` で処理する核種の一覧（all/<library か control>/nuclides.txt）を書く
- `python3 all_nuclides.py depletion` は all/<燃料>/nominal.toml（trial/<燃料>/nominal.toml のライブラリを all/library に替えたもの）を書く
"""

import argparse
import sys
from pathlib import Path

import msgpack
import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from read_survey import nuclide_key

FUELS = ("u", "mox", "gd")
# 試行のセル（inputs.py）の燃料、被覆管、減速材の温度 (K)。lattice は S(α, β) の核種のほかは表の温度と一致する温度だけを
# 引く（crates/format/src/library.rs の `lookup`）ので、この 1 点で試行のセルの計算は 13 温度のライブラリと同じになる
TEMPERATURE = 565
# 走査（../read-survey-endfb80.toml）で xs が止まる、チェーンの評価
EXCLUDED = {
    "H-2": ("#0311", "MF6 の中性子の LAW=6 を読めない"),
    "Be-9": ("#0311", "MF6 の中性子の LAW=7 を読めない"),
    "Th-232": (
        "#0312",
        "核分裂の中性子の分布が MF6 MT18 で、輸送する核種として処理できない",
    ),
    "Pa-231": (
        "#0312",
        "核分裂の中性子の分布が MF6 MT18 で、輸送する核種として処理できない",
    ),
    "Pa-233": (
        "#0312",
        "核分裂の中性子の分布が MF6 MT18 で、輸送する核種として処理できない",
    ),
    "Cl-35": ("#0264", "分離共鳴（R 行列の制限形）の荷電粒子の対を読めない"),
    "Ca-40": ("#0264", "分離共鳴（R 行列の制限形）の荷電粒子の対を読めない"),
}
SURVEY = HERE.parent / "read-survey-endfb80.toml"
CHAIN = HERE.parent / "chain.msgpack"
TRIAL = HERE / "trial" / "library" / "xs.toml"
P2 = HERE.parent / "p2" / "xs.toml"


def chain_evaluations(survey: dict, chain: dict) -> dict[str, tuple[str, str | None]]:
    """走査の評価のうちチェーンの核種のもの: 核種の名前 → (ファイルの名前, 止まった理由か None)。走査の順"""
    if chain.get("format") != "xlc-chain":
        raise SystemExit("not an xlc chain")
    names = {(n["z"], n["a"], n["state"]): n["name"] for n in chain["nuclides"]}
    out = {}
    for e in survey["evaluation"]:
        name = names.get(nuclide_key(e["endf"]))
        if name is None:
            continue
        if name in out:
            raise SystemExit(f"{name} has two evaluations in the survey")
        out[name] = (Path(e["endf"]).name, e.get("error"))
    return out


def library_entries(
    base: dict, survey: dict, chain: dict, excluded: dict
) -> dict[str, dict]:
    """全核種のライブラリの核種の項目（試行の入力 `base` の順、続いて `base` に無い評価を走査の順）"""
    evaluations = chain_evaluations(survey, chain)
    stopped = {n for n, (_, error) in evaluations.items() if error is not None}
    if stopped != set(excluded):
        raise SystemExit(
            f"the exclusions differ from the evaluations xs stops on in {sorted(stopped ^ set(excluded))}"
        )
    files = {f: n for n, (f, _) in evaluations.items()}
    unknown = [n for n, v in base.items() if files.get(v["endf"]) != n]
    if unknown:
        raise SystemExit(f"{unknown} of the base input are not chain evaluations")
    out = {}
    for name, entry in base.items():
        if name in excluded and not entry.get("activation"):
            raise SystemExit(f"{name} is excluded but the base input transports it")
        if name in excluded or "thermal" not in entry:
            entry = {**entry, "temperatures": [TEMPERATURE]}
        if name not in excluded:
            entry = {k: v for k, v in entry.items() if k != "activation"}
        out[name] = entry
    for name, (endf, _) in evaluations.items():
        if name not in out and name not in excluded:
            out[name] = {"endf": endf, "temperatures": [TEMPERATURE]}
    return out


def transported(entries: dict) -> list[str]:
    return [n for n, v in entries.items() if not v.get("activation")]


def control_entries(base: dict, entries: dict) -> dict[str, dict]:
    """`entries` の核種の項目のうち、試行の入力 `base` が輸送しない核種を放射化の核種にしたもの"""
    keep = set(transported(base))
    return {
        n: v if n in keep else {**v, "activation": True} for n, v in entries.items()
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
    return str(v)


def render(base_text: str, head: list[str], entries: dict) -> str:
    """注釈 `head`、試行の入力の全核種に共通の項目、核種の項目（[nuclide] のインラインの表）の入力"""
    common = []
    for line in base_text.splitlines():
        if line.startswith("["):
            break
        if line and not line.startswith("#"):
            common.append(line)
    rows = [
        f"{n} = {{ " + ", ".join(f"{k} = {value(v)}" for k, v in e.items()) + " }"
        for n, e in entries.items()
    ]
    return "\n".join([*(f"# {h}" for h in head), *common, "", "[nuclide]", *rows, ""])


def exclusion_notes(base: dict, excluded: dict) -> list[str]:
    """外した核種ごとの注釈（試行の入力 `base` で放射化の核種なら、残すことを書く。`library_entries` と同じ判定）"""
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


LIBRARY_HEAD = [
    "#0292 の試行: 燃焼チェーンにあって ENDF/B-VIII.0 の中性子のファイルがある核種を、すべて輸送する核種にしたライブラリ。",
    "all_nuclides.py が ../../trial/library/xs.toml と ../../../read-survey-endfb80.toml から書く。手で直さない。",
    f"温度は試行のセルの {TEMPERATURE} K の 1 点（H-1 は S(α, β) の内挿に使う 4 点）。理由は all_nuclides.py の TEMPERATURE。",
    "xs が止まる評価は外す（外した核種と issue）:",
]
CONTROL_HEAD = [
    "#0292 の試行の対照: ../library/xs.toml と同じ核種の項目で、154 核種の試行のライブラリ（../../trial/library）が",
    "輸送する核種だけを輸送し、ほかを放射化の核種にしたもの。範囲（xs-range）とビルドを全核種のライブラリとそろえ、",
    "全核種を輸送した差から、範囲と xs の版の違いの分を分ける。all_nuclides.py が書く。手で直さない。",
    "nuclides.txt の放射化の核種だけを ../library/xs-range.toml の範囲で xs-nuclides に作り、組み立ては、../library/xs-nuclides の上に",
    "xs-nuclides を重ねた置き場で ../library/xs-range.toml と行う（輸送する核種は全核種のライブラリの xs のファイル）。",
]


def write_library() -> None:
    base_text = TRIAL.read_text()
    base = tomllib.loads(base_text)["nuclide"]
    entries = library_entries(
        base,
        tomllib.loads(SURVEY.read_text()),
        msgpack.unpackb(CHAIN.read_bytes(), strict_map_key=False),
        EXCLUDED,
    )
    notes = [f"  {n}" for n in exclusion_notes(base, EXCLUDED)]
    control = control_entries(base, entries)
    jobs = job_lists(entries, control)
    for directory, head, items in [
        ("library", LIBRARY_HEAD + notes, entries),
        ("control", CONTROL_HEAD, control),
    ]:
        path = HERE / "all" / directory / "xs.toml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render(base_text, head, items))
        (path.parent / "nuclides.txt").write_text("\n".join(jobs[directory]) + "\n")
    p2 = tomllib.loads(P2.read_text())["nuclide"]
    for name, names in nuclide_lists(entries, base, p2).items():
        (HERE / "all" / f"{name}.txt").write_text("\n".join(names) + "\n")


def depletion_input(text: str) -> str:
    """試行の燃焼の入力 `text` のライブラリを all/library に替えたもの"""
    old = 'library = "../library/depletion.msgpack"'
    lines = text.splitlines()
    if not lines[0].startswith("#") or lines[1] != old:
        raise SystemExit("the trial input has another head or library")
    head = "# #0292 の試行: ../../trial の同じ燃料の nominal.toml の、ライブラリだけを全核種を輸送する ../library に替えたもの。all_nuclides.py が書く。手で直さない"
    return "\n".join([head, old, *lines[2:]]) + "\n"


def write_depletion() -> None:
    for fuel in FUELS:
        directory = HERE / "all" / fuel
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "nominal.toml").write_text(
            depletion_input((HERE / "trial" / fuel / "nominal.toml").read_text())
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["library", "depletion"])
    args = parser.parse_args()
    if args.command == "library":
        write_library()
    else:
        write_depletion()


if __name__ == "__main__":
    main()
