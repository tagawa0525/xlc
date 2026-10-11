"""軸方向の反射体の平板（ADR-0128）の等価性: 平板の格子計算（MOC）と、その定数で解いた炉心の計算（NEM）の k と
ノード平均の中性子束を比べる。どちらも `RUST_LOG=debug` のログを読む。

- 格子計算: k は1ノードの問題のログ（`one-node constants` の `k`。MOC の k の全桁）、ノード平均の中性子束は `node` の
  `average_flux`（ノードの順。全体の生成率 1 で、奥行き 1 cm あたり）
- 炉心の計算: k は最後の反復の `k`、ノード平均の中性子束は `node average flux`（ノードの順。全体の生成率 1。ピッチを
  平板の行の高さ（1 cm）にしてあるので、格子計算と同じ規格化）

`fstar` を付けると、少数群定数のファイルの z の辺（−z、+z）の ADF と f* の印を並べる。

使い方:
  python3 compare_axial.py <格子計算のログ> <炉心の計算のログ>
  python3 compare_axial.py fstar <少数群定数のファイル>...
"""

import re
import sys
from pathlib import Path

import msgpack

NUMBER = r"[-+0-9.eE]+"


def ordered(log: str, fluxes: dict[int, list[float]]) -> list[list[float]]:
    """ノードの番号の順の中性子束。番号は 0 から抜けなく並ぶこと（1つの状態の、1回の計算のログ）"""
    if sorted(fluxes) != list(range(len(fluxes))):
        raise SystemExit(
            f"{log}: node numbers {sorted(fluxes)} are not 0 to {len(fluxes) - 1}"
        )
    return [fluxes[n] for n in range(len(fluxes))]


def lattice(log: str) -> tuple[float, list[list[float]]]:
    k = None
    fluxes = {}
    for line in Path(log).read_text().splitlines():
        if "one-node constants" in line:
            k = float(re.search(rf" k=({NUMBER})", line).group(1))
        m = re.search(rf"node node=(\d+) average_flux=\[({NUMBER}), ({NUMBER})\]", line)
        if m:
            fluxes[int(m.group(1))] = [float(m.group(2)), float(m.group(3))]
    if k is None or not fluxes:
        raise SystemExit(f"{log}: no k or node fluxes (run with RUST_LOG=debug)")
    return k, ordered(log, fluxes)


def core(log: str) -> tuple[float, list[list[float]]]:
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
    return k, ordered(log, fluxes)


def compare(lattice_log: str, core_log: str) -> None:
    k_moc, moc = lattice(lattice_log)
    k_core, nem = core(core_log)
    if len(moc) != len(nem):
        raise SystemExit(f"{len(moc)} lattice nodes, {len(nem)} core nodes")
    print(
        f"k MOC {k_moc:.8f}, core {k_core:.8f}, dk/k {1e5 * (k_core / k_moc - 1):+.4f} pcm"
    )
    worst = [0.0, 0.0]
    for n, (a, b) in enumerate(zip(nem, moc)):
        d = [x / y - 1 for x, y in zip(a, b)]
        print(
            f"node {n:2}: MOC ({b[0]:.6e}, {b[1]:.6e}), core / MOC - 1 ({d[0]:+.2e}, {d[1]:+.2e})"
        )
        worst = [max(w, abs(x)) for w, x in zip(worst, d)]
    print(f"max |core / MOC - 1| (fast, thermal): ({worst[0]:.2e}, {worst[1]:.2e})")


def fstar(files: list[str]) -> None:
    for f in files:
        table = msgpack.unpackb(Path(f).read_bytes(), raw=False)
        for state in table["states"]:
            adf = state["adf"]
            marks = state["normalized_adf"]
            sides = ", ".join(
                f"{name} ({adf[s][0]:.5f}, {adf[s][1]:.5f}){' f*' if marks[s] else ''}"
                for name, s in [("-z", 4), ("+z", 5)]
            )
            print(f"{Path(f).name}: height {state['height']}, {sides}")


def main() -> None:
    if sys.argv[1:2] == ["fstar"]:
        fstar(sys.argv[2:])
    elif len(sys.argv) == 3:
        compare(sys.argv[1], sys.argv[2])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
