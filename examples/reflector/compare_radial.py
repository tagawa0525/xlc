"""径方向の反射体の格子（I 型、角型）の MOC と、その定数で解いた炉心の計算（NEM）の k とノード平均の中性子束を
比べる（ADR-0132）。どちらも `RUST_LOG=debug` のログを読む（読み方は compare_axial.py と同じ）。

- 格子計算: k は1ノードの問題のログ（`one-node constants` の `k`）、ノード平均の中性子束は `node` の `average_flux`
  （格子のノードの順。全体の生成率 1 で、奥行き 1 cm あたり）
- 炉心の計算: k は最後の反復の `k`、ノード平均の中性子束は `node average flux`（直方体の位置の順。全体の生成率 1）。
  炉心の高さ（z の幅）を掛けて、格子計算と同じ奥行き 1 cm あたりにする

使い方:
  python3 compare_radial.py <格子計算のログ> <炉心の計算のログ> <炉心の高さ (cm)> <位置>...
位置は、格子のノード 0、1、… に当たる炉心の位置の番号（x が最も速く、次に y）。
"""

import re
import sys
from pathlib import Path

from compare_axial import NUMBER, lattice


def main() -> None:
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    lattice_log, core_log, height = sys.argv[1], sys.argv[2], float(sys.argv[3])
    sites = [int(s) for s in sys.argv[4:]]
    k_moc, moc = lattice(lattice_log)
    k_core, nem = core_fluxes(core_log)
    if len(sites) != len(moc):
        raise SystemExit(f"{len(moc)} lattice nodes, {len(sites)} sites given")
    print(
        f"k MOC {k_moc:.8f}, core {k_core:.8f}, dk/k {1e5 * (k_core / k_moc - 1):+.4f} pcm"
    )
    worst = [0.0, 0.0]
    for n, (b, site) in enumerate(zip(moc, sites)):
        a = [height * x for x in nem[site]]
        d = [x / y - 1 for x, y in zip(a, b)]
        print(
            f"node {n} (site {site}): MOC ({b[0]:.6e}, {b[1]:.6e}), core / MOC - 1 ({d[0]:+.2e}, {d[1]:+.2e})"
        )
        worst = [max(w, abs(x)) for w, x in zip(worst, d)]
    print(f"max |core / MOC - 1| (fast, thermal): ({worst[0]:.2e}, {worst[1]:.2e})")


def core_fluxes(log: str) -> tuple[float, dict[int, list[float]]]:
    """炉心の計算の k と、位置ごとのノード平均の中性子束（ノードのない位置は含まない）"""
    k = None
    fluxes = {}
    for line in Path(log).read_text().splitlines():
        m = re.search(rf"iteration=\d+ k=({NUMBER})", line)
        if m:
            k = float(m.group(1))
        m = re.search(
            rf"node average flux node=(\d+) flux=\[({NUMBER}), ({NUMBER})\]", line
        )
        if m:
            fluxes[int(m.group(1))] = [float(m.group(2)), float(m.group(3))]
    if k is None or not fluxes:
        raise SystemExit(f"{log}: no k or node fluxes (run with RUST_LOG=debug)")
    return k, fluxes


if __name__ == "__main__":
    main()
