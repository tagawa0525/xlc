"""#0180 の自己遮蔽を省いたときの誤差の測定の入力を書く。

燃焼の軸（2A の分岐のテーブル。../branch/2A-60.toml の公称の 1 状態を、ADR-0199 の 27 燃焼点で解く）:

- 2A-every.toml: すべての燃焼点で自己遮蔽を解く（参照。既定）
- 2A-skip-2.toml: 27 点の 1 つおき（最初の点から）だけ自己遮蔽を解き、ほかの点はその前の点の組成で自己遮蔽する（`shielding-burnup`）
- 2A-skip-sparse.toml: SPARSE の 8 点だけで解く（間の点は、前の解く点からの燃焼度の差が 0.1〜12 MWd/kg）

微分の摂動の軸（ADR-0221 の Pu-240 の数密度の微分。燃焼点 DERIVATIVE_BURNUP）:

- 2A-derivative.toml: 摂動の組成で自己遮蔽を解き直す（既定）
- 2A-derivative-fixed.toml: すべての点を `shielding-burnup` に挙げる。自己遮蔽は状態の組成のまま、摂動は解き直さない

状態の軸の暗黙の効き（ADR-0247 の段 4 の `[sensitivity.few-group]` の `shielded`。../../issue-0279/inputs.py の 60 MWd/kg の組成の
ピンセルと 3×3）:

- pincell-moderator.toml、lattice-3x3-moderator.toml: 減速材の H-1 とホウ素の B-10 の数密度（減速材の密度とホウ素濃度の代わり）

燃焼の軸の 1 次の見込みに使う、ほかの燃焼度の燃料の暗黙の効き:

- pincell-few-group-<燃焼度>.toml: ../../issue-0279/inputs.py の pincell-few-group.toml の組成を、履歴の燃焼点 FUEL_BURNUP の
  ものにしたもの（自己遮蔽を解き直す核種は SHIELDED。60 MWd/kg で 1 割の変化に 1.3 pcm より大きい暗黙の効きを持つ核種と、
  燃焼の初めに変わる Xe-135、Sm-149）

使い方: python3 inputs.py（このディレクトリで）
"""

import importlib.util
import sys
from pathlib import Path

import msgpack

HERE = Path(__file__).resolve().parent
HISTORY = HERE / "../../p9-product/lattice-2A.msgpack"
BRANCH = HERE / "../branch/2A-60.toml"

# ADR-0199 の 27 点（../../p9-2A/branch-2A.toml の burnup）
BURNUP = [
    0, 0.1, 0.25, 0.5, 0.75, 1, 2, 4, 5, 6,
    7, 8, 9, 10, 11, 12, 14, 16, 20, 24,
    28, 32, 36, 40, 48, 56, 60,
]  # fmt: skip
SPARSE = [0, 0.5, 2, 6, 10, 16, 28, 48]
DERIVATIVE_BURNUP = [4, 20, 60]
FUEL_BURNUP = [4.0, 20.0]
SHIELDED = [
    "U-235", "U-236", "U-238", "Pu-239", "Pu-240", "Pu-241", "Pu-242", "Sm-150", "Xe-135", "Sm-149",
]  # fmt: skip


def stage4():
    """#0279 の段 4 の入力の部品（../../issue-0279/inputs.py）"""
    path = HERE / "../../issue-0279/inputs.py"
    spec = importlib.util.spec_from_file_location("issue_0279_inputs", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def listed(values) -> str:
    return (
        "["
        + ", ".join(
            f"{v:g}" if isinstance(v, float | int) else f'"{v}"' for v in values
        )
        + "]"
    )


def branch(text: str, output: str, burnup: list, extra: str, derivative: bool) -> str:
    """../branch/2A-60.toml の、先頭の注釈（その入力の説明）を除き、テーブルの出力、燃焼点、微分の有無を替え、`[history]` の
    最後に `extra` を足したもの"""
    start = text.find("library = ")
    if start < 0:
        sys.exit(f"{BRANCH} has no library line")
    replaced = text[start:].replace(
        'output = "table-2A-60.msgpack"', f'output = "{output}"'
    )
    replaced = replaced.replace("burnup = [60]\n", f"burnup = {listed(burnup)}\n")
    if not derivative:
        replaced = replaced.replace('derivative = ["Pu-240"]\n', "")
    if replaced.count(output) != 1 or f"burnup = {listed(burnup)}" not in replaced:
        sys.exit(f"{BRANCH} does not have the lines to replace")
    return replaced.rstrip("\n") + "\n" + extra


def fuels_at(history: dict, burnup: float) -> str:
    """燃焼点 `burnup` の、最初の燃焼の領域の組の 5 つのリングの物質（../../issue-0279/inputs.py の `fuels` と同じ形）"""
    point = next((p for p in history["points"] if p["burnup"] == burnup), None)
    if point is None:
        sys.exit(f"the history has no point at {burnup} MWd/kg")
    transportable = set(history["transportable"])
    first = history["regions"][0]["cell"]
    rings = [
        (region["material"], densities)
        for region, densities in zip(history["regions"], point["densities"])
        if region["cell"] == first
    ]
    if len(rings) != 5:
        sys.exit("the first region group must have the five fuel rings")
    kept = [
        i
        for i, name in enumerate(history["nuclides"])
        if name in transportable and any(d[i] > 0.0 for _, d in rings)
    ]
    lines = []
    for material, densities in rings:
        lines.append(
            f"[material.{material}]\ntemperature = 900\n\n[material.{material}.nuclides]"
        )
        for i in kept:
            lines.append(f"{history['nuclides'][i]} = {max(densities[i], 0.0)!r}")
        lines.append("")
    return "\n".join(lines)


def few_group(output: str, nuclides: list[str], shielded: list[str]) -> str:
    """段 4 と同じ直接の摂動と判定の `[sensitivity.few-group]`（../../issue-0279/inputs.py の `few_group`）"""
    return (
        "[sensitivity]\n"
        f"nuclides = {listed(nuclides)}\n"
        f'output = "{output}"\n'
        "direct = 0.1\n"
        "k-tolerance = 1e-12\n"
        "flux-tolerance = 1e-11\n\n"
        "[sensitivity.few-group]\n"
        "bounds = [1e7, 0.625, 1e-5]\n"
        "tolerance = 1e-10\n"
        "restart = 100\n"
        f"shielded = {listed(shielded)}\n"
    )


def main() -> None:
    parts = stage4()
    header = parts.HEADER.replace(
        '"../p2/depletion.msgpack"', '"../../p2/depletion.msgpack"'
    )
    if header == parts.HEADER:
        sys.exit("the stage-4 header has no library line to replace")
    with HISTORY.open("rb") as f:
        history = msgpack.unpackb(f.read(), strict_map_key=False)
    text = BRANCH.read_text()
    comment = "# このファイルは inputs.py が書く（手で直さない）。#0180 の自己遮蔽を省いたときの誤差\n"
    every_other = BURNUP[::2]
    cases = {
        "2A-every.toml": branch(text, "table-2A-every.msgpack", BURNUP, "", False),
        "2A-skip-2.toml": branch(
            text,
            "table-2A-skip-2.msgpack",
            BURNUP,
            f"shielding-burnup = {listed(every_other)}\n",
            False,
        ),
        "2A-skip-sparse.toml": branch(
            text,
            "table-2A-skip-sparse.msgpack",
            BURNUP,
            f"shielding-burnup = {listed(SPARSE)}\n",
            False,
        ),
        "2A-derivative.toml": branch(
            text, "table-2A-derivative.msgpack", DERIVATIVE_BURNUP, "", True
        ),
        "2A-derivative-fixed.toml": branch(
            text,
            "table-2A-derivative-fixed.msgpack",
            DERIVATIVE_BURNUP,
            f"shielding-burnup = {listed(DERIVATIVE_BURNUP)}\n",
            True,
        ),
    }
    for name, body in cases.items():
        (HERE / name).write_text(comment + body)
    at_60 = parts.fuels(history) + "\n" + parts.ZIRCALOY
    moderator = ["H-1", "B-10"]
    (HERE / "pincell-moderator.toml").write_text(
        comment
        + header
        + 'rows = ["pin"]\n\n'
        + at_60
        + "\n"
        + parts.PIN
        + "\n"
        + few_group("pincell-moderator-sensitivity.toml", moderator, moderator)
    )
    (HERE / "lattice-3x3-moderator.toml").write_text(
        comment
        + header
        + 'rows = ["pin pin pin", "pin gt  pin", "pin pin pin"]\n\n'
        + at_60
        + "\n"
        + parts.PIN
        + "\n"
        + parts.GUIDE_TUBE
        + "\n"
        + few_group("lattice-3x3-moderator-sensitivity.toml", moderator, moderator)
    )
    for burnup in FUEL_BURNUP:
        name = f"pincell-few-group-{burnup:g}"
        (HERE / f"{name}.toml").write_text(
            comment
            + header
            + 'rows = ["pin"]\n\n'
            + fuels_at(history, burnup)
            + "\n"
            + parts.ZIRCALOY
            + "\n"
            + parts.PIN
            + "\n"
            + few_group(f"{name}-sensitivity.toml", parts.NUCLIDES, SHIELDED)
        )


if __name__ == "__main__":
    main()
