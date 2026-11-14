"""#0285 の段 1（ADR-0250）: 比の状態のログ（target ratio_check の DEBUG。xlc の `log_ratio_check`）から、比で求めた参照ごとに、
燃料のリングの捕獲の状態と基準の比 − 1 を、捕獲の反応率（`rate`）で重みを付けて平均したものを出す。

比べる比（どれも UFG の範囲の群）:

- equivalence: 等価原理の比（ADR-0250 の比。ピンに 1 つでリングによらない）
- state: 状態の UFG のリングの比（本当の変化）
- mean: 状態の UFG のピンの平均の比（リングの間の分布を基準のままにした形の比）
- dancoff: 基準の物質のまま状態の Dancoff 係数で解き直した UFG のリングの比（減速材の物質の変化を除いた分）

リングごと（核種ごと）と、U-238 の群ごと（リングを足す）に出す。反応率がどれも 0 の核種（個数密度が 0）は出さない。参照は、ログの比の自己遮蔽の区間（`ratio-shielding` の close）で
区切った順。

使い方: python3 ratio_check.py <ログ>
"""

import re
import sys
from collections import defaultdict
from pathlib import Path

FIELD = re.compile(r'(\w+)=("[^"]*"|\S+)')
RATIOS = ["equivalence", "state", "mean", "dancoff"]


def blocks(path: Path) -> list[dict]:
    """参照ごとの {"rings": {(リング, 核種): 平均}, "groups": {群: 平均}}。平均は {比の名前: 比 − 1}"""
    out = []
    rings: dict = defaultdict(lambda: defaultdict(float))
    groups: dict = defaultdict(lambda: defaultdict(float))

    def close() -> dict:
        def mean(sums: dict) -> dict:
            # 反応率がどの群でも 0 の核種（個数密度が 0）は、重みの平均が無いので出さない
            return {
                key: {q: s[q] / s["w"] - 1.0 for q in RATIOS}
                for key, s in sums.items()
                if s["w"] > 0.0
            }

        return {"rings": mean(rings), "groups": mean(groups)}

    for line in path.read_text().splitlines():
        if "ratio check of the ring capture" in line:
            f = {k: v.strip('"') for k, v in FIELD.findall(line)}
            w = float(f["rate"])
            keys = [(rings, (int(f["ring"]), f["nuclide"]))]
            if f["nuclide"] == "U-238":
                keys.append((groups, int(f["group"])))
            for sums, key in keys:
                for q in RATIOS:
                    sums[key][q] += w * float(f[q])
                sums[key]["w"] += w
        elif " INFO ratio-shielding: " in line and "close time" in line and rings:
            out.append(close())
            rings, groups = (
                defaultdict(lambda: defaultdict(float)),
                defaultdict(lambda: defaultdict(float)),
            )
    return out


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    for b, block in enumerate(blocks(Path(sys.argv[1]))):
        print(f"reference {b}")
        for (ring, nuclide), r in sorted(block["rings"].items()):
            print(
                f"  ring {ring} {nuclide:7s} "
                + " ".join(f"{q} {r[q]:+.3%}" for q in RATIOS)
            )
        for group, r in sorted(block["groups"].items()):
            print(
                f"  U-238 group {group:2d} "
                + " ".join(f"{q} {r[q]:+.3%}" for q in RATIOS)
            )


if __name__ == "__main__":
    main()
