"""#0293 の項目 3: 出力の規格化だけが違う OpenMC のピンセルの燃焼 2 本（openmc_depletion.py の extract の msgpack）を比べる。

`fission-q`（チェーンの核分裂の Q。ADR-0173）と `energy-deposition`（`heating-local` の tally。捕獲の γ を含む）の 2 本を、
同じ燃焼度の点で比べる。燃焼点ごとに書くもの:

- `delta`: energy-deposition の、決めた出力 / チェーンの Q の核分裂の出力 − 1（Q の外に数えたエネルギーの割合。estimates.py の
  `capture_over_fission` に当たる）。fission-q の同じ比は 1 で、外れれば止める（バランス）
- `mean_delta`、`burnup_shift`: δ の燃焼度 0 からの台形の平均 δ̄ と、燃焼度の軸のずれ B δ̄
- `predicted_pcm`: estimates.py の `add_burnup_shift` と同じ 1 次の見込み −(dk/dB) B δ̄ / k を、fission-q の履歴の k と、
  この δ で出したもの（dk/dB は 1 つ前の点との差分）
- `measured_pcm`、`measured_std_pcm`: (k_ed − k_fq) / k_fq と、2 本の k の標準誤差の二乗和の平方根 / k_fq
- `densities`: 主な核種の、リングの体積で平均した数密度の相対差（energy-deposition / fission-q − 1）

2 本が規格化のほかの前提を共有しないもの（燃焼度の点、輸送する核種、決めた出力）は止める。

使い方: python3 compare_normalization.py <fission-q の msgpack> <energy-deposition の msgpack> [--json <出力>]
"""

import argparse
import json
import math
import sys
from pathlib import Path

import msgpack

NUCLIDES = ["U235", "Pu239", "Pu240", "Pu241", "Xe135", "Sm149"]
# fission-q の決めた出力とチェーンの Q の核分裂の出力の相対差の許容（extract の和の丸めだけが残る）
CLOSURE = 1.0e-6


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def averaged(run: dict, nuclide: str, i: int) -> float:
    volume = sum(r["volume"] for r in run["rings"])
    return sum(r["volume"] * r["densities"][nuclide][i] for r in run["rings"]) / volume


def compare(fission_q: dict, deposition: dict) -> list[dict]:
    for key in ("burnup", "transportable", "power"):
        if fission_q[key] != deposition[key]:
            sys.exit(f"the two runs differ in {key}")
    power = fission_q["power"]
    for b, p in zip(fission_q["burnup"], fission_q["fission-power"]):
        if abs(p / power - 1) > CLOSURE:
            sys.exit(
                f"the fission-q run is off its set power at {b} MWd/kg: {p} / {power}"
            )
    burnup = fission_q["burnup"]
    deltas = [power / p - 1 for p in deposition["fission-power"]]
    rows = []
    area = 0.0
    for i, b in enumerate(burnup):
        k = fission_q["k"][i]
        if i > 0:
            area += (deltas[i - 1] + deltas[i]) / 2 * (b - burnup[i - 1])
        mean = area / b if b > 0 else None
        predicted = None
        if i > 0 and b > 0:
            slope = (k - fission_q["k"][i - 1]) / (b - burnup[i - 1])
            predicted = -slope * b * mean / k * 1e5
        std = math.hypot(fission_q["k-std"][i], deposition["k-std"][i])
        densities = {}
        for nuclide in NUCLIDES:
            low = averaged(fission_q, nuclide, i)
            densities[nuclide] = (
                averaged(deposition, nuclide, i) / low - 1 if low > 0 else None
            )
        rows.append(
            {
                "burnup": b,
                "k_fission_q": k,
                "k_energy_deposition": deposition["k"][i],
                "delta": deltas[i],
                "mean_delta": mean,
                "burnup_shift": b * mean if mean is not None else None,
                "predicted_pcm": predicted,
                "measured_pcm": (deposition["k"][i] - k) / k * 1e5,
                "measured_std_pcm": std / k * 1e5,
                "densities": densities,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fission_q", type=Path)
    parser.add_argument("deposition", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rows = compare(read(args.fission_q), read(args.deposition))
    for r in rows:
        shift = (
            ""
            if r["predicted_pcm"] is None
            else (
                f" mean {r['mean_delta'] * 100:.2f}% shift {r['burnup_shift']:.2f} MWd/kg"
                f" predicted {r['predicted_pcm']:+.0f} pcm"
            )
        )
        dens = " ".join(
            f"{n} {v * 100:+.2f}%" for n, v in r["densities"].items() if v is not None
        )
        print(
            f"{r['burnup']:5.1f} MWd/kg delta {r['delta'] * 100:.2f}%{shift}"
            f" measured {r['measured_pcm']:+.0f} ± {r['measured_std_pcm']:.0f} pcm | {dens}"
        )
    if args.json:
        args.json.write_text(json.dumps({"rows": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
