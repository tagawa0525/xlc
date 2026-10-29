"""小さな炉心の core の燃焼（examples/depletion/verification.md の「ステップの大きさの誤差と中点の臨界ホウ素（小さな炉心）」）の
ケース（branches/core-small-*.toml）を、基準のサイクルの出力のあるステップを等しい EFPD に分けて書く。

使い方: python3 small_core_cases.py [--check]

- 基準のサイクル（`CYCLE`）の出力のあるステップ（`efpd`）を、ケースごとの分割の数の等しい EFPD のステップに分ける（崩壊だけの
  ステップは分けない）。分けたステップは、もとのステップの出力の割合を持ち、状態はもとのステップの最初のステップにだけ書く
- restart と装荷替え（ADR-0220）: 基準のサイクルと 30 日の停止を切らずに燃焼させるケース（`core-small-restart-full.toml`）、それを 2 つの
  点で切った前と後の半分（`core-small-restart-<a|b><1|2>.toml`。後の半分は前の半分の履歴から restart するサイクル
  `core-small-restart-<a|b>2-cycle.toml`）、切らない燃焼の履歴から装荷を替えた 1 サイクル（`core-small-reload.toml` と
  `core-small-reload-cycle.toml`）
- `--check`: 書く内容がファイルと同じかを確かめ、違えば止める（ファイルを書かない）
"""

import sys
from pathlib import Path

BRANCHES = Path(__file__).resolve().parent / "branches"

# 定格の出力密度 (W/g)
SPECIFIC_POWER = 38.41

# 基準のサイクル: (種類（"burnup" か "decay"）, 炉心の平均の燃焼度 (MWd/kg) か停止の日数, 出力の割合（None は書かない）, 状態（None は前と
# 同じ）, コメント)。燃焼度は EFPD = 燃焼度 × 1000 / 出力密度
CYCLE = [
    ("burnup", 0.25, None, "aro", "0.25 MWd/kg"),
    ("burnup", 1.0, None, None, "1 MWd/kg"),
    ("burnup", 2.0, None, None, "2 MWd/kg"),
    ("burnup", 3.0, None, None, "3 MWd/kg"),
    ("burnup", 4.0, None, "rodded", "4 MWd/kg"),
    ("burnup", 5.0, None, None, "5 MWd/kg"),
    ("burnup", 6.0, 0.7, "aro", "6 MWd/kg"),
    ("burnup", 7.0, 0.7, None, "7 MWd/kg"),
    ("decay", 18.0, None, None, "停止"),
    ("burnup", 7.5, None, None, "7.5 MWd/kg"),
]

TRACKED = 'micro = { chain = "../chain.msgpack", xenon = "tracked" }'

# 燃焼の方式: 名前 → (ヘッダーの 4 行目、`micro` の行、`[depletion]` に足す行)
METHODS = {
    "search": ("中点のホウ素は ADR-0191 の探索。", TRACKED, []),
    "start": (
        "中点のホウ素は ステップの初めの点で探した濃度（ADR-0191 の案 S1）。",
        TRACKED,
        ['midpoint-boron = "start"'],
    ),
    "equilibrium": (
        "中点のホウ素は ADR-0191 の探索。Xe-135 と I-135 は出力のある解で平衡の値（#0200）。",
        'micro = { chain = "../chain.msgpack", xenon = "equilibrium" }',
        [],
    ),
}

# ケース: (方式, 分割)
CASES = (
    [("search", s) for s in (1, 2, 4, 8, 16, 32, 64, 128)]
    + [("start", s) for s in (1, 2, 4)]
    + [("equilibrium", s) for s in (1, 2, 4, 8, 16)]
)


def head(cycle: str) -> str:
    """ケースの初めの行（サイクルのファイル `cycle`、メッシュ、境界）"""
    return f"""\
cycle = "{cycle}"
axial-mesh = [0.0, 22.86, 45.72, 68.58, 91.44, 114.3, 137.16, 160.02, 182.88, 205.74, 228.6, 251.46, 274.32, 297.18, 320.04, 342.9, 365.76]
radial-nodes = 1
boundary = {{ x = ["reflective", "vacuum"], y = ["reflective", "vacuum"], z = ["marshak", "marshak"], outline = "marshak" }}

"""


# 状態の名前 → `[state.<名前>]` の表（コメントの行を含む）
STATES = {
    "aro": """\
# バンク D を抜いた状態（230 ステップで先端が炉心の上端）
[state.aro]
b10 = 0.199
moderator = 565
fuel = 565
search = { boron = [0.0, 1300.0], k = 1.071, tolerance = 1.0e-6, max-solves = 20 }
banks = { D = 230 }
""",
    "rodded": """\
# バンク D を 140 ステップ（先端 222.885 cm）まで入れた状態。下から 10 番目のノード（205.74〜228.6 cm）が部分挿入
[state.rodded]
b10 = 0.199
moderator = 565
fuel = 565
search = { boron = [0.0, 1300.0], k = 1.071, tolerance = 1.0e-6, max-solves = 20 }
banks = { D = 140 }
""",
}

DEPLETION = """\
[depletion]
specific-power = 38.41
"""


def states(names: list[str]) -> str:
    """状態 `names` の表を、空行で区切って並べたもの（燃焼の入力は、どのステップも指さない状態をエラーにする）"""
    return "\n".join(STATES[n] for n in names) + "\n"


HEAD = head("core-small-cycle.toml") + states(["aro", "rodded"]) + DEPLETION


def steps(split: int) -> list[str]:
    """基準のサイクルの出力のあるステップを `split` 個に分けたステップの行"""
    lines = []
    start = 0.0
    for kind, value, power, state, comment in CYCLE:
        if kind == "decay":
            lines.append(f"  {{ decay-days = {value!r} }},  # {comment}")
            continue
        end = value * 1000 / SPECIFIC_POWER
        for k in range(1, split + 1):
            efpd = end if k == split else start + (end - start) * k / split
            items = [f"efpd = {efpd!r}"]
            if power is not None:
                items.append(f"power = {power!r}")
            if state is not None and k == 1:
                items.append(f'state = "{state}"')
            line = "  { " + ", ".join(items) + " },"
            lines.append(line + (f"  # {comment}" if k == split else ""))
        start = end
    return lines


def case(method: str, split: int) -> str:
    sentence, micro, extra = METHODS[method]
    head = (
        "# 小さな炉心の core の燃焼（ステップの大きさの誤差と中点の探索の効き）: core-small-plant.toml の 2A の 4 × 4 の 1/4 と反射体を、\n"
        "# バンク D の出し入れ、出力の割合、18 日の停止を含むサイクルで、点と中点の臨界ホウ素を探して燃焼させる。基準のサイクル\n"
        f"# （core-small-search-1.toml）の出力のあるステップを等しい EFPD の {split} 個に分ける（崩壊だけのステップは分けない）。\n"
        f"# {sentence}\n"
    )
    depletion = [
        f'history = "core-small-{method}-{split}.msgpack"',
        "# 格子計算の燃焼（../p2/lattice-2A.toml）と同じチェーン",
        micro,
        *extra,
        "steps = [",
        *steps(split),
        "]",
    ]
    return head + HEAD + "\n".join(depletion) + "\n"


# restart と装荷替え（ADR-0220。examples/depletion/verification.md の「restart と装荷替え（小さな炉心）」）

# サイクルの間の停止（Problem 10 の 30 日）。前のサイクルのケースの最後の崩壊だけのステップに置く（ADR-0220 の案 S1）
OUTAGE = ("decay", 30.0, None, None, "サイクルの間の停止")

# 切らない燃焼: 基準のサイクル（分割 1）と停止
RESTART_CYCLE = [*CYCLE, OUTAGE]

# 切る点: 名前 → 前の半分のステップの数。a は 2 MWd/kg の点（前後のステップが同じ状態と出力）、b は 18 日の停止の後
CUTS = {"a": 3, "b": 9}

# 装荷替えのサイクル（燃焼度はこのサイクルの炉心の平均の増分。EFPD は 0 から数える）
RELOAD_CYCLE = [
    ("burnup", 0.25, None, "aro", "0.25 MWd/kg"),
    ("burnup", 1.0, None, None, "1 MWd/kg"),
    ("burnup", 2.0, None, None, "2 MWd/kg"),
    ("burnup", 3.0, None, None, "3 MWd/kg"),
]

EQUILIBRIUM = 'micro = { chain = "../chain.msgpack", xenon = "equilibrium" }'

LOADING = """\
loading = [
  ".    .    .",
  "A01  A02  .",
  "A03  A04  .",
]
"""

# A01 と A04 を入れ替え、A02 を 90° 回し、A03 を取り出して新しい A05 を中心側（バンク D の位置）に装荷する
RELOAD_LOADING = """\
loading = [
  ".    .       .",
  "A04  A02@90  .",
  "A05  A01     .",
]
"""


def efpd_of(burnup: float) -> float:
    """炉心の平均の燃焼度の増分 (MWd/kg) の累積の EFPD"""
    return burnup * 1000 / SPECIFIC_POWER


def cycle_steps(cycle: list) -> list[dict]:
    """サイクル `cycle` の、ステップごとの種類、終わりの累積の EFPD（崩壊だけのステップは前と同じ）、崩壊の日数、出力の割合、効いている状態、
    状態を書いたか、コメント"""
    out = []
    efpd = 0.0
    state = None
    for kind, value, power, named, comment in cycle:
        state = named or state
        if kind == "burnup":
            efpd = efpd_of(value)
        out.append(
            {
                "kind": kind,
                "efpd": efpd,
                "days": value if kind == "decay" else None,
                "power": power,
                "state": state,
                "named": named is not None,
                "comment": comment,
            }
        )
    return out


def restart_lines(steps: list[dict], start: float) -> list[str]:
    """ステップ `steps` の行。累積の EFPD は `start` から数え直し、最初のステップには効いている状態を書く"""
    lines = []
    for k, step in enumerate(steps):
        items = []
        if step["kind"] == "decay":
            items.append(f"decay-days = {step['days']!r}")
        else:
            items.append(f"efpd = {step['efpd'] - start!r}")
            if step["power"] is not None:
                items.append(f"power = {step['power']!r}")
        if k == 0 or step["named"]:
            items.append(f'state = "{step["state"]}"')
        lines.append("  { " + ", ".join(items) + " },  # " + step["comment"])
    return lines


def used_states(steps: list[dict]) -> list[str]:
    return sorted({s["state"] for s in steps})


def restart_case(
    name: str, comment: str, cycle: str, steps: list[dict], start: float
) -> str:
    """restart の確認のケース `core-small-<name>.toml`（平衡 Xe。中点のホウ素は ADR-0191 の探索）"""
    depletion = [
        f'history = "core-small-{name}.msgpack"',
        "# 格子計算の燃焼（../p2/lattice-2A.toml）と同じチェーン",
        EQUILIBRIUM,
        "steps = [",
        *restart_lines(steps, start),
        "]",
    ]
    return (
        comment
        + head(cycle)
        + states(used_states(steps))
        + DEPLETION
        + "\n".join(depletion)
        + "\n"
    )


def restart_cycle(comment: str, restart: str, loading: str) -> str:
    return (
        comment
        + 'plant = "core-small-plant.toml"\n'
        + 'inventory = "core-small-assemblies.toml"\n'
        + f"restart = {restart}\n"
        + loading
    )


def restart_files() -> dict[str, str]:
    """restart と装荷替えのケースとサイクル: ファイル名 → 内容"""
    steps = cycle_steps(RESTART_CYCLE)
    files = {
        "core-small-restart-full.toml": restart_case(
            "restart-full",
            "# 小さな炉心の restart の確認（ADR-0220）の、切らない燃焼: 基準のサイクル（core-small-equilibrium-1.toml）と、その後の\n"
            "# 30 日の停止（サイクルの間の停止。装荷替えの core-small-reload.toml の前のサイクル）。平衡 Xe\n",
            "core-small-cycle.toml",
            steps,
            0.0,
        )
    }
    for cut, n in CUTS.items():
        start = steps[n - 1]["efpd"]
        files[f"core-small-restart-{cut}1.toml"] = restart_case(
            f"restart-{cut}1",
            f"# 小さな炉心の restart の確認（ADR-0220）: core-small-restart-full.toml の燃焼を {steps[n - 1]['comment']} の点で切った前の半分\n",
            "core-small-cycle.toml",
            steps[:n],
            0.0,
        )
        files[f"core-small-restart-{cut}2-cycle.toml"] = restart_cycle(
            f"# 小さな炉心の restart の確認（ADR-0220）: core-small-cycle.toml と同じ装荷で、core-small-restart-{cut}1.toml の燃焼の履歴の\n"
            "# 最後の点から続ける\n",
            f'{{ history = "core-small-restart-{cut}1.msgpack" }}',
            LOADING,
        )
        files[f"core-small-restart-{cut}2.toml"] = restart_case(
            f"restart-{cut}2",
            f"# 小さな炉心の restart の確認（ADR-0220）: core-small-restart-full.toml の燃焼を {steps[n - 1]['comment']} の点で切った後の半分。\n"
            f"# 累積の EFPD は切った点（{start!r} EFPD）から数え直す\n",
            f"core-small-restart-{cut}2-cycle.toml",
            steps[n:],
            start,
        )
    files["core-small-reload-cycle.toml"] = restart_cycle(
        "# 小さな炉心の装荷替え（ADR-0220）: core-small-restart-full.toml（基準のサイクルと 30 日の停止）の燃焼の履歴の最後の点から、\n"
        "# A01 と A04 を入れ替え、A02 を 90° 回し、A03 を取り出して新しい A05 を中心側（バンク D の位置）に装荷する\n",
        '{ history = "core-small-restart-full.msgpack", fresh = ["A05"] }',
        RELOAD_LOADING,
    )
    files["core-small-reload.toml"] = restart_case(
        "reload",
        "# 小さな炉心の装荷替え（ADR-0220）の 1 サイクル: core-small-reload-cycle.toml の装荷を、バンク D を抜いて定格で、サイクルの\n"
        "# 炉心の平均の燃焼度の増分 3 MWd/kg まで燃焼させる。平衡 Xe\n",
        "core-small-reload-cycle.toml",
        cycle_steps(RELOAD_CYCLE),
        0.0,
    )
    return files


def main() -> None:
    args = sys.argv[1:]
    if args not in ([], ["--check"]):
        raise SystemExit(__doc__)
    stale = []
    files = {f"core-small-{m}-{s}.toml": case(m, s) for m, s in CASES} | restart_files()
    for name, text in files.items():
        path = BRANCHES / name
        if args:
            if not path.exists() or path.read_text() != text:
                stale.append(path.name)
        else:
            path.write_text(text)
    if stale:
        raise SystemExit(
            f"stale cases {stale}; run python3 examples/depletion/small_core_cases.py"
        )


if __name__ == "__main__":
    main()
