"""#0285（ADR-0250）: 比の状態のログ（target ratio_check の DEBUG。xlc の `log_transfer_check`）から、比で求めた参照ごとに、
等価なセルのリングの移行の行の、比で作った行と状態の UFG の行の差を反応率で数える。

ログの行は、ピン（代表のピン。対称操作のコピーの数で重みを付けない）、リング、核種、入射の群ごとに、群から出る移行 `out`、そのうち
上の群への移行 `up`、UFG の範囲の下の群への移行 `below`（微視的断面積）を、状態の UFG の行（state）、弾性散乱の比で縮めた行
（elastic）、静止した標的の行の比で作った行（static）で持つ。`weight`（数密度 × 中性子束 × 体積）を掛けてピンを足し、

- (リング, 核種) ごとと (リング, 核種, 群) ごとの、作った行と状態の行の反応率の差 Σ weight × (作った行 − state)
- 状態の行の反応率 Σ weight × state（相対の差の分母）

を出す。参照は、ログの比の自己遮蔽の区間（`ratio-shielding` の close）で区切った順。

使い方: python3 transfer_check.py <ログ> [<表示する (リング, 核種, 群) の数>]
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

FIELD = re.compile(r'(\w+)=("[^"]*"|\S+)')
ROWS = ["elastic", "static"]
SUMS = ["out", "up", "below"]


def blocks(path: Path) -> list[dict]:
    """参照ごとの {"rings": {(リング, 核種): 和}, "groups": {(リング, 核種, 群): 和}}。和は {"state_<量>": 状態の反応率,
    "<行>_<量>": 作った行と状態の行の反応率の差}（量は SUMS、行は ROWS）"""
    out = []

    def empty():
        return defaultdict(lambda: defaultdict(float))

    rings, groups = empty(), empty()
    for line in path.read_text().splitlines():
        if "ratio check of the ring transfer" in line:
            f = {k: v.strip('"') for k, v in FIELD.findall(line)}
            w = float(f["weight"])
            ring, nuclide, group = int(f["ring"]), f["nuclide"], int(f["group"])
            for sums in (rings[(ring, nuclide)], groups[(ring, nuclide, group)]):
                for q in SUMS:
                    state = float(f[f"state_{q}"])
                    sums[f"state_{q}"] += w * state
                    for r in ROWS:
                        sums[f"{r}_{q}"] += w * (float(f[f"{r}_{q}"]) - state)
        elif " INFO ratio-shielding: " in line and "close time" in line and rings:
            out.append({"rings": rings, "groups": groups})
            rings, groups = empty(), empty()
    return out


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    top = int(sys.argv[2]) if len(sys.argv) == 3 else 10

    def show(label: str, s: dict) -> str:
        parts = []
        for q in SUMS:
            parts.append(
                f"{q} state {s[f'state_{q}']:.4e} "
                + " ".join(f"{r} {s[f'{r}_{q}']:+.3e}" for r in ROWS)
            )
        return f"  {label}: " + "; ".join(parts)

    for b, block in enumerate(blocks(Path(sys.argv[1]))):
        print(f"reference {b}")
        for (ring, nuclide), s in sorted(block["rings"].items()):
            print(show(f"ring {ring} {nuclide:7s}", s))
        ranked = sorted(
            block["groups"].items(), key=lambda kv: -abs(kv[1]["elastic_out"])
        )
        for (ring, nuclide, group), s in ranked[:top]:
            print(show(f"ring {ring} {nuclide:7s} group {group:2d}", s))


if __name__ == "__main__":
    main()
