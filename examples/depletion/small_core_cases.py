"""小さな炉心の core の燃焼（examples/depletion/verification.md の「ステップの大きさの誤差と中点の臨界ホウ素（小さな炉心）」）の
ケース（branches/core-small-*.toml）を、基準のサイクルの出力のあるステップを等しい EFPD に分けて書く。

使い方: python3 small_core_cases.py [--check]

- 基準のサイクル（`CYCLE`）の出力のあるステップ（`efpd`）を、ケースごとの分割の数の等しい EFPD のステップに分ける（崩壊だけの
  ステップは分けない）。分けたステップは、もとのステップの出力の割合を持ち、状態はもとのステップの最初のステップにだけ書く
- `--check`: 書く内容がファイルと同じかを確かめ、違えば止める（ファイルを書かない）
"""

import sys
from pathlib import Path

BRANCHES = Path(__file__).resolve().parent / "branches"

# 定格の出力密度 (W/g)
SPECIFIC_POWER = 38.41

# 基準のサイクル: (炉心の平均の燃焼度 (MWd/kg) か停止の日数, 出力の割合, 状態, コメント)。燃焼度は EFPD = 燃焼度 × 1000 / 出力密度
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


HEAD = """\
cycle = "core-small-cycle.toml"
axial-mesh = [0.0, 22.86, 45.72, 68.58, 91.44, 114.3, 137.16, 160.02, 182.88, 205.74, 228.6, 251.46, 274.32, 297.18, 320.04, 342.9, 365.76]
radial-nodes = 1
boundary = { x = ["reflective", "vacuum"], y = ["reflective", "vacuum"], z = ["marshak", "marshak"], outline = "marshak" }

# バンク D を抜いた状態（230 ステップで先端が炉心の上端）
[state.aro]
b10 = 0.199
moderator = 565
fuel = 565
search = { boron = [0.0, 1300.0], k = 1.071, tolerance = 1.0e-6, max-solves = 20 }
banks = { D = 230 }

# バンク D を 140 ステップ（先端 222.885 cm）まで入れた状態。下から 10 番目のノード（205.74〜228.6 cm）が部分挿入
[state.rodded]
b10 = 0.199
moderator = 565
fuel = 565
search = { boron = [0.0, 1300.0], k = 1.071, tolerance = 1.0e-6, max-solves = 20 }
banks = { D = 140 }

[depletion]
specific-power = 38.41
"""


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


def main() -> None:
    args = sys.argv[1:]
    if args not in ([], ["--check"]):
        raise SystemExit(__doc__)
    stale = []
    for method, split in CASES:
        path = BRANCHES / f"core-small-{method}-{split}.toml"
        text = case(method, split)
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
