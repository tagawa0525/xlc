"""#0293 の項目 3: 出力の規格化だけが違う OpenMC のピンセルの燃焼 2 本（openmc_depletion.py の extract の msgpack）を比べる。

`fission-q`（チェーンの核分裂の Q。ADR-0173）と `energy-deposition`（`heating-local` の tally。捕獲の γ を含む）の 2 本を、
同じ燃焼度の点で比べる。燃焼点ごとに書くもの:

- `delta`: energy-deposition の、決めた出力 / チェーンの Q の核分裂の出力 − 1（Q の外に数えたエネルギーの割合。estimates.py の
  `capture_over_fission` に当たる）。fission-q の同じ比は 1 で、外れれば止める（バランス）
- `mean_delta`、`burnup_shift`: δ の燃焼度 0 からの台形の平均 δ̄ と、燃焼度の軸のずれ B δ̄
- `predicted_pcm`: estimates.py の `add_burnup_shift` と同じ 1 次の見込み −(dk/dB) B δ̄ / k を、fission-q の履歴の k と、
  この δ で出したもの（dk/dB は 1 つ前の点との差分）。初めの数点の dk/dB は Xe と Sm の立ち上がり（出力と時間で決まり、2 本で同じ）を
  含むので、見込みは燃焼度のずれの効きを大きく見積もる
- `measured_pcm`、`measured_std_pcm`: (k_ed − k_fq) / k_fq と、2 本の k の標準誤差の二乗和の平方根 / k_fq。2 本は同じシードで
  輸送の相関が正（燃焼度 0 は組成も同じで k が一致する）なので、差の標準誤差の上限
- `densities`: 主な核種の、リングの体積で平均した数密度の相対差（energy-deposition / fission-q − 1）

2 本がそれぞれの規格化でないもの（extract が書く `normalization`）と、規格化のほかの前提を共有しないもの（燃焼度の点、輸送する
核種、決めた出力、チェーンの核分裂の Q）は止める。

`--kappa`（normalization/casl-kappa.toml）を渡すと、VERA の燃焼のベンチマークの仕様書の
核種ごとの Kappa / Kappa0 − 1 を、fission-q の履歴の核種ごとの核分裂の出力（率 × extract が書いたチェーンの Q）で重み付けた割合（`kappa_delta`）と、
表に無い核種の核分裂の出力の割合（`kappa_uncovered`）も書く。核種ごとの定数で規格化したときに足す割合で、`delta` と並べる。

使い方: python3 compare_normalization.py <fission-q の msgpack> <energy-deposition の msgpack> [--kappa <toml>] [--json <出力>]
"""

import argparse
import json
import math
import sys
from pathlib import Path

import msgpack
import tomllib

NUCLIDES = ["U235", "Pu239", "Pu240", "Pu241", "Xe135", "Sm149"]
# fission-q の決めた出力とチェーンの Q の核分裂の出力の相対差の許容（extract の和の丸めだけが残る）
CLOSURE = 1.0e-6


def read(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def averaged(run: dict, nuclide: str, i: int) -> float:
    volume = sum(r["volume"] for r in run["rings"])
    return sum(r["volume"] * r["densities"][nuclide][i] for r in run["rings"]) / volume


def compare(fission_q: dict, deposition: dict) -> list[dict]:
    for run, mode in ((fission_q, "fission-q"), (deposition, "energy-deposition")):
        if run["normalization"] != mode:
            sys.exit(f"expected a {mode} run, got {run['normalization']}")
    for key in ("burnup", "transportable", "power", "fission-q-values"):
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


def kappa_delta(run: dict, kappa: dict) -> list[dict]:
    """燃焼点ごとの、核分裂の出力（率 × チェーンの Q）で重み付けた Kappa / Kappa0 − 1 と、`kappa` に無い核種の出力の割合。"""
    rates = run["fission-rates"]
    q = run["fission-q-values"]
    rows = []
    for i in range(len(run["burnup"])):
        total = covered = extra = 0.0
        for name, rate in rates.items():
            power = rate[i] * q[name]
            total += power
            if name in kappa:
                covered += power
                extra += power * (kappa[name]["kappa"] / kappa[name]["kappa0"] - 1)
        rows.append({"delta": extra / covered, "uncovered": 1 - covered / total})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fission_q", type=Path)
    parser.add_argument("deposition", type=Path)
    parser.add_argument("--kappa", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    fission_q = read(args.fission_q)
    rows = compare(fission_q, read(args.deposition))
    if args.kappa is not None:
        kappa = tomllib.loads(args.kappa.read_text())
        for r, k in zip(rows, kappa_delta(fission_q, kappa)):
            r["kappa_delta"] = k["delta"]
            r["kappa_uncovered"] = k["uncovered"]
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
            f" measured {r['measured_pcm']:+.0f} ± {r['measured_std_pcm']:.0f} pcm"
            + (
                f" kappa {r['kappa_delta'] * 100:.2f}% (uncovered {r['kappa_uncovered']:.1e})"
                if "kappa_delta" in r
                else ""
            )
            + f" | {dens}"
        )
    if args.json:
        args.json.write_text(json.dumps({"rows": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
