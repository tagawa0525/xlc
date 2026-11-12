"""燃料のピンの UFG の減速計算の 1 次の変化（lattice の `[sensitivity.ufg-tangent]`。ADR-0247 の段 I2-2）を、状態の間で比べる
（#0279）。暗黙の効きの係数を燃焼点ごとに公称の状態で 1 度だけ求める近似が成り立つかの材料。

向き（核種）と反応ごとに、燃料のリングの UFG の範囲の群の巨視的断面積の暗黙の分 T_s（数密度を (1 + δ) 倍したときの δ に
ついての微分）を、その状態の中性子束 φ_s で重みを付けて比べる:

- 相対の暗黙の分 A_s = Σ φ_s T_s / Σ φ_s Σ_s（リングと群の和。ピンの UFG の範囲の反応率の、暗黙の分の相対の変化）
- 公称の状態の係数を使ったときのずれ E_s = ‖φ_s (T_s − T_0)‖ / ‖φ_s T_s‖（2 ノルム）。T_0 は公称の状態の値
- その向きの暗黙の分が全体の暗黙の分に占める大きさ ‖φ_0 T_0‖ / max_d ‖φ_0 T_0,d‖（反応ごと。ずれの小さい向きと、効きの小さい向きを
  分けて読む）

使い方: python3 compare_ufg_tangent.py <公称の結果の TOML> <状態の結果の TOML>...
"""

import json
import math
import sys
from pathlib import Path

import tomllib


def pin(report: dict) -> dict:
    """結果のピン（1 本であること）"""
    pins = report["ufg_tangent"]["pins"]
    if len(pins) != 1:
        sys.exit(f"expected one pin, got {len(pins)}")
    return pins[0]


def weighted(flux: list[list[float]], values: list[list[float]]) -> list[float]:
    """[リング][群] の φ × 値を並べたもの"""
    return [
        f * v for ring_f, ring_v in zip(flux, values) for f, v in zip(ring_f, ring_v)
    ]


def norm(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values))


def compare(nominal: dict, states: dict[str, dict]) -> dict:
    """公称の状態 `nominal` と状態 `states`（名前と結果）の、向きと反応ごとの A_s、E_s と、公称の状態の大きさ"""
    raise NotImplementedError


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    paths = [Path(p) for p in sys.argv[1:]]
    reports = [tomllib.loads(p.read_text()) for p in paths]
    states = {
        p.stem.removesuffix("-sensitivity"): r for p, r in zip(paths[1:], reports[1:])
    }
    print(json.dumps(compare(reports[0], states), indent=1))


if __name__ == "__main__":
    main()
