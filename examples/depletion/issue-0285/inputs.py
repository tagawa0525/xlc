"""#0285 の段 1（ADR-0250）の入力を書く: 分岐の状態の自己遮蔽を基準の UFG と等価原理の比で求めたテーブルと、状態ごとに UFG で
解いたテーブル（`every-state`）を、燃料の種類と燃焼度を変えて比べる。

燃焼の履歴の入力（`[depletion]` より前）の格子、物質、離散化に、テーブルの出力（ピン出力の形状関数つき）と分岐と `[history]` を
足す。状態は ADR-0239 の減速材の軸と燃料の温度の軸の、両端と公称（減速材 560、580、615 K、燃料 560、900、1350 K）の 3 × 3 で、
ホウ素は公称の 650 ppm の 1 点。物質の状態は、どの履歴でも ADR-0239 の公称（燃料 900 K、被覆管と減速材 580 K、密度 0.711901 g/cc、
ホウ素 650 ppm）に書き直す（比の基準は入力の物質の減速材の点で、本番の基準と同じ点にするため。p2 の 2O の履歴は 600 K で燃焼させた
もので、組成だけを使う）。

入力の名前は <種類>-<燃焼度>-<解き方>.toml:

- every: 状態ごとに UFG を解く（`every-state`。比べる基準）
- f1、f2: ADR-0250 の (F1)（`moderator-ratio`）と (F2)（`moderator-fuel-ratio`）で、比が決めない量はどれも基準の値。f2 は
  CANDIDATES の燃焼点だけ
- f1-transfer、f2-transfer: 移行を弾性散乱の比で縮めたもの（基準の移行の行のままでは、行の和が比で変えた弾性散乱と食い違い、
  その差が吸収に入る。2A の 60 MWd/kg の Slurm 863）
- f1-factor、f1-dilute、f1-sph: f1-transfer の、燃料でないリングを自己遮蔽因子の比にする、無限希釈にする、SPH 因子を状態の UFG の
  ものにしたもの。CANDIDATES の燃焼点だけ
- f1-average、f1-state: f1-transfer の燃料のリングを、基準のリングに状態と基準の UFG のピンの平均の比を掛けたもの（リングの間の
  分布は基準のまま）、状態の UFG のリングの値にしたもの（比の残る差を、平均の変化とリングの間の分布と燃料でないリングに分ける
  診断）。CANDIDATES の燃焼点だけ
- f1-state-sph: f1-state の SPH 因子も状態の UFG のもの（残る差を燃料でないリングだけにする診断）。CANDIDATES の燃焼点だけ
- f1-exact、f1-state-sections: f1-state-sph の燃料でないリングも状態の UFG のもの（f1-exact は every と同じになるはずの確認）と、
  その群定数だけを状態の UFG にし、移行は基準の行を弾性散乱の比で縮めたもの。CANDIDATES の燃焼点だけ
- f1-static-fuel、f1-static: f1-transfer の燃料のリング（f1-static は燃料でないリングと SPH 因子も）を、基準の値に、状態と基準を
  静止した標的で解いた減速計算の比を掛けたもの。CANDIDATES の燃焼点だけ

時間の内訳の入力 timing-<種類>-<燃焼度>-<解き方>.toml（TIMING）は、状態を燃料 900 K の減速材 3 点に絞る（基準の UFG 1 回と、同じ
基準の比の状態 2 回。ログの水準を上げて段ごとの時間を読む）。解き方には SHIELDINGS のほかに、時間だけを測る TIMING_SHIELDINGS:

- f1-static-cp6: f1-static の静止した標的の減速計算の衝突確率の求積を、小区間あたり 6 点にしたもの（`static-cp-points = 6`）

使い方: python3 inputs.py（このディレクトリで）
"""

import re
import sys
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
DEPLETION = HERE / ".."

# 種類: （履歴の入力（examples/depletion から）、分岐を解く燃焼点 (MWd/kg)）。Gd の 2O は Gd が残る 10 MWd/kg
CASES = {
    "2A": ("p9-product/lattice-2A.toml", [0, 20, 60]),
    "2O": ("p2/lattice-2O.toml", [10]),
    "MOX": ("mox-assembly/lattice-mox-assembly.toml", [0, 20, 60]),
}
# 比が決めない量の候補を比べる燃焼点
CANDIDATES = {"2A": [60], "2O": [10], "MOX": [60]}

# ADR-0239 の軸の両端と公称（../p9-2A/branch-2A.toml の格子点）
MODERATOR = ([560, 580, 615], [0.752074, 0.711901, 0.608215])
FUEL = [560, 900, 1350]
NOMINAL_FUEL, NOMINAL_MODERATOR = 900, 580
WATER = "water = { density = 0.711901, boron = 650, b10 = 0.199 }"

SHIELDINGS = {
    "every": ('shielding = "every-state"', ""),
    "f1": ('shielding = "moderator-ratio"', ""),
    "f1-transfer": ('shielding = "moderator-ratio"', 'transfer = "elastic-ratio"'),
    "f2": ('shielding = "moderator-fuel-ratio"', ""),
    "f2-transfer": ('shielding = "moderator-fuel-ratio"', 'transfer = "elastic-ratio"'),
    "f1-factor": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nrings = "shielding-factor"',
    ),
    "f1-dilute": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nrings = "infinite-dilution"',
    ),
    "f1-sph": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nsph = "state"',
    ),
    "f1-average": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "state-average"',
    ),
    "f1-state": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "state"',
    ),
    "f1-state-sph": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "state"\nsph = "state"',
    ),
    "f1-exact": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "state"\nrings = "state"\nsph = "state"',
    ),
    "f1-state-sections": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "state"\nrings = "state-sections"\nsph = "state"',
    ),
    "f1-static-fuel": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "static-ratio"',
    ),
    "f1-static": (
        'shielding = "moderator-ratio"',
        'transfer = "elastic-ratio"\nfuel = "static-ratio"\nrings = "static-ratio"\nsph = "static-ratio"',
    ),
}
# 時間の内訳を測る入力: (種類, 燃焼度, 解き方)。状態は TIMING_STATES
TIMING = [("2A", 60, "f1-static"), ("2A", 60, "f1-static-cp6")]
TIMING_STATES = (MODERATOR, [NOMINAL_FUEL])
TIMING_SHIELDINGS = {
    "f1-static-cp6": (
        'shielding = "moderator-ratio"',
        (
            'transfer = "elastic-ratio"\nfuel = "static-ratio"\nrings = "static-ratio"\nsph = "static-ratio"\n'
            "static-cp-points = 6"
        ),
    ),
}
# どの燃焼点でも流すもの（ほかは CANDIDATES の燃焼点だけ）
MAIN = ["every", "f1", "f1-transfer", "f2-transfer"]


def listed(values) -> str:
    return "[" + ", ".join(f"{v:g}" for v in values) + "]"


def with_temperature(body: str, material: str, kelvin: int) -> str:
    """`[material.<material>]` の直後の `temperature = ` の行を `kelvin` にする"""
    pattern = re.compile(
        rf"(\[material\.{re.escape(material)}\]\ntemperature = )[0-9.]+\n"
    )
    if len(pattern.findall(body)) != 1:
        sys.exit(f"material {material} has no single temperature line")
    return pattern.sub(rf"\g<1>{kelvin}\n", body)


def branch(
    name: str,
    source: str,
    burnup: float,
    shielding: str,
    states=(MODERATOR, FUEL),
    prefix: str = "",
) -> str:
    path = DEPLETION / source
    text = path.read_text()
    start, end = text.find("library = "), text.find("[depletion]")
    if start < 0 or end < 0:
        sys.exit(f"{path} has no library line or [depletion]")
    data = tomllib.loads(text)
    output, depleted, library = (
        data["depletion"]["output"],
        data["depletion"]["materials"],
        data["library"],
    )
    stem = f"{prefix}{name}-{burnup:g}-{shielding}"
    moderator, fuel = states

    def relative(p: str) -> Path:
        return (path.parent / p).resolve().relative_to(HERE.resolve(), walk_up=True)

    body = text[start:end].replace(
        f'library = "{library}"',
        f'library = "{relative(library)}"\n'
        f'few-group = {{ output = "table-{stem}.msgpack", bounds = [1e7, 0.625, 1e-5], pins = true }}',
        1,
    )
    if "few-group" in data or body.count("few-group = ") != 1:
        sys.exit(f"{path} has a few-group line or no library line")
    # 減速材は公称のホウ酸水（p2 の履歴は数密度で書いてあるので書き直す）
    water = re.compile(
        r"\[material\.water\]\ntemperature = [0-9.]+\n(?:#[^\n]*\n)*(?:water|nuclides) = \{[^\n]*\}\n"
    )
    if len(water.findall(body)) != 1:
        sys.exit(f"{path} has no single water material")
    body = water.sub(
        f"[material.water]\ntemperature = {NOMINAL_MODERATOR}\n{WATER}\n", body
    )
    body = with_temperature(body, "zircaloy", NOMINAL_MODERATOR)
    for material in depleted:
        body = with_temperature(body, material, NOMINAL_FUEL)
    kind, ratio = (SHIELDINGS | TIMING_SHIELDINGS)[shielding]
    ratio = f"\n[branches.ratio]\n{ratio}\n" if ratio else ""
    fuels = ", ".join(f'"{m}"' for m in depleted)
    return (
        f"# このファイルは inputs.py が書く（手で直さない）。#0285 の段 1（ADR-0250）: {source} の {burnup:g} MWd/kg、{shielding}\n"
        f"{body.rstrip()}\n\n"
        f'[branches]\nboron = {{ ppm = [650], materials = ["water"] }}\n'
        f"moderator = {{ temperature = {listed(moderator[0])}, density = {listed(moderator[1])}, "
        f'materials = ["water", "zircaloy"] }}\n'
        f"fuel = {{ temperature = {listed(fuel)}, materials = [{fuels}] }}\n"
        f"{kind}\n{ratio}\n"
        f'[history]\nfile = "{relative(output)}"\nchain = "../chain.msgpack"\nburnup = [{burnup:g}]\n'
    )


def names() -> list[str]:
    """書く入力の名前（拡張子なし）。比べる順"""
    out = []
    for name, (_, burnups) in CASES.items():
        for burnup in burnups:
            for shielding in SHIELDINGS:
                if shielding in MAIN or burnup in CANDIDATES[name]:
                    out.append(f"{name}-{burnup:g}-{shielding}")
    return out


def main() -> None:
    for stem in names():
        name, burnup, shielding = stem.split("-", 2)
        source = CASES[name][0]
        (HERE / f"{stem}.toml").write_text(
            branch(name, source, float(burnup), shielding)
        )
    for name, burnup, shielding in TIMING:
        (HERE / f"timing-{name}-{burnup:g}-{shielding}.toml").write_text(
            branch(
                name,
                CASES[name][0],
                float(burnup),
                shielding,
                TIMING_STATES,
                "timing-",
            )
        )


if __name__ == "__main__":
    main()
