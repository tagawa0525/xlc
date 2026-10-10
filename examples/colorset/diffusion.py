"""均質化した定数の炉心（外周は鏡面反射、軸方向は一様）を、細かいメッシュの差分の拡散で解き、core と同じ形の
ピン出力を書く。core の NEM とピン出力の再構成（ADR-0096）と同じ問題を、空間の離散化だけを細かくして解くので、
core との差は、NEM と再構成の方法の誤差になる。

使い方: python3 diffusion.py <cycle の台帳と並び> <出力の pins.msgpack> <セルあたりのメッシュの数>

- 並びは `cycle.toml` の loading と台帳（`inventory`）から読む（回さない並びだけ）
- ノードは集合体1体（分割なし）。定数は各ノードのテーブルの最初の状態
- ノードの間の面では、不均質な中性子束（均質な中性子束 × その辺の ADF）と正味の中性子流を連続とする
- ピン出力は Σ_g Σf,g φ_g（セルの平均）f_g（セル）。core と同じく、すべて 0 の行と列を除き、核分裂のあるセルの平均で
  1 に規格化する。行は +y（上）から
"""

import sys
from pathlib import Path

import msgpack
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import tomllib


def read_layout(cycle_path: Path):
    cycle = tomllib.loads(cycle_path.read_text())
    inventory = tomllib.loads((cycle_path.parent / cycle["inventory"]).read_text())
    rows = [line.split() for line in cycle["loading"]]
    nodes = []
    for row in rows:
        out = []
        for name in row:
            if "@" in name:
                raise SystemExit(f"turned assemblies are not handled: {name}")
            design = inventory["design"][inventory["assembly"][name]["design"]]
            if len(design["layers"]) != 1:
                raise SystemExit("only one axial layer is handled")
            table = msgpack.unpackb(
                (cycle_path.parent / design["layers"][0]).read_bytes(), raw=False
            )
            if table["split"] != 1:
                raise SystemExit("split nodes are not handled")
            out.append(table["states"][0])
        nodes.append(out)
    return nodes


def solve(nodes, per_cell: int):
    """`nodes`[行（上から）][列] の定数で解き、(k, 中性子束 [群][y（下から）][x], メッシュの情報) を返す。"""
    n_rows, n_cols = len(nodes), len(nodes[0])
    first = nodes[0][0]["pins"]
    if first is None:
        raise SystemExit("the constants have no form functions")
    groups = len(first["form"])

    # x の方向のメッシュ: ノードの列ごとに、セルの幅を per_cell 個に分ける（幅の小さい余白のセルは 1 個）
    def mesh(widths):
        h, cell = [], []
        for c, w in enumerate(widths):
            m = per_cell if w > 0.5 else 1
            h += [w / m] * m
            cell += [c] * m
        return h, cell

    xs, ys = [], []  # [(ノードの列, セルの列, 幅)]
    for col in range(n_cols):
        h, cell = mesh(nodes[0][col]["pins"]["dx"])
        xs += [(col, c, w) for w, c in zip(h, cell)]
    for row in reversed(range(n_rows)):  # y は下から
        h, cell = mesh(nodes[row][0]["pins"]["dy"])
        ys += [(row, c, w) for w, c in zip(h, cell)]
    nx, ny = len(xs), len(ys)
    size = groups * nx * ny

    def index(g, j, i):
        return (g * ny + j) * nx + i

    loss = sp.lil_matrix((size, size))
    fission = sp.lil_matrix((size, size))
    for j, (row, _, hy) in enumerate(ys):
        for i, (col, _, hx) in enumerate(xs):
            m = nodes[row][col]
            for g in range(groups):
                k = index(g, j, i)
                loss[k, k] += (m["total"][g] - m["scatter"][g][g]) * hx * hy
                for src in range(groups):
                    if src != g:
                        loss[k, index(src, j, i)] -= m["scatter"][src][g] * hx * hy
                    fission[k, index(src, j, i)] += (
                        m["chi"][g] * m["nu_fission"][src] * hx * hy
                    )
    # 面の流れ: J = (f_a φ_a − f_b φ_b) / (f_a h_a / 2D_a + f_b h_b / 2D_b)。ノードの中では f = 1
    sides = {"left": 0, "right": 1, "top": 2, "bottom": 3}

    def couple(g, a, b, length):
        (ka, ma, ha, fa), (kb, mb, hb, fb) = a, b
        c = length / (
            fa * ha / (2 * ma["diffusion"][g]) + fb * hb / (2 * mb["diffusion"][g])
        )
        loss[ka, ka] += c * fa
        loss[ka, kb] -= c * fb
        loss[kb, kb] += c * fb
        loss[kb, ka] -= c * fa

    for g in range(groups):
        for j, (row, _, hy) in enumerate(ys):
            for i in range(nx - 1):
                (ca, _, ha), (cb, _, hb) = xs[i], xs[i + 1]
                ma, mb = nodes[row][ca], nodes[row][cb]
                between = ca != cb
                fa = ma["adf"][sides["right"]][g] if between else 1.0
                fb = mb["adf"][sides["left"]][g] if between else 1.0
                couple(
                    g,
                    (index(g, j, i), ma, ha, fa),
                    (index(g, j, i + 1), mb, hb, fb),
                    hy,
                )
        for i, (col, _, hx) in enumerate(xs):
            for j in range(ny - 1):
                (ra, _, ha), (rb, _, hb) = ys[j], ys[j + 1]  # j は下、j + 1 は上
                ma, mb = nodes[ra][col], nodes[rb][col]
                between = ra != rb
                fa = ma["adf"][sides["top"]][g] if between else 1.0
                fb = mb["adf"][sides["bottom"]][g] if between else 1.0
                couple(
                    g,
                    (index(g, j, i), ma, ha, fa),
                    (index(g, j + 1, i), mb, hb, fb),
                    hx,
                )

    lu = spla.splu(loss.tocsc())
    fission = fission.tocsr()
    phi = np.ones(size)
    k = 1.0
    for _ in range(10000):
        source = fission @ phi
        new = lu.solve(source / k)
        new_k = k * (fission @ new).sum() / source.sum()
        change = np.abs(new / np.abs(new).max() - phi / np.abs(phi).max()).max()
        phi, k_old, k = new, k, new_k
        if abs(k - k_old) < 1e-10 and change < 1e-9:
            break
    else:
        raise SystemExit("power iteration did not converge")
    return k, phi.reshape(groups, ny, nx), xs, ys


def pin_powers(nodes, phi, xs, ys):
    groups = phi.shape[0]
    n_rows, n_cols = len(nodes), len(nodes[0])
    # セル（ノードの行、ノードの列、セルの行（下から）、セルの列）ごとの中性子束の面積の平均
    cols = [
        (col, c)
        for col in range(n_cols)
        for c in range(len(nodes[0][col]["pins"]["dx"]))
    ]
    rows = [
        (row, c)
        for row in reversed(range(n_rows))
        for c in range(len(nodes[row][0]["pins"]["dy"]))
    ]
    col_of = {key: n for n, key in enumerate(cols)}
    row_of = {key: n for n, key in enumerate(rows)}
    flux = np.zeros((groups, len(rows), len(cols)))
    area = np.zeros((len(rows), len(cols)))
    for j, (row, cy, hy) in enumerate(ys):
        for i, (col, cx, hx) in enumerate(xs):
            r, c = row_of[(row, cy)], col_of[(col, cx)]
            flux[:, r, c] += phi[:, j, i] * hx * hy
            area[r, c] += hx * hy
    flux /= area
    power = np.zeros((len(rows), len(cols)))
    for r, (row, cy) in enumerate(rows):
        for c, (col, cx) in enumerate(cols):
            m = nodes[row][col]
            form = m["pins"]["form"]
            power[r, c] = sum(
                m["fission"][g] * flux[g, r, c] * form[g][cy][cx] for g in range(groups)
            )
    power = power[::-1]  # 行は上から
    power = power[power.any(axis=1)][:, power.any(axis=0)]
    return power / power[power > 0].mean()


def main():
    cycle, output, per_cell = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    nodes = read_layout(cycle)
    k, phi, xs, ys = solve(nodes, per_cell)
    power = pin_powers(nodes, phi, xs, ys)
    output.write_bytes(
        msgpack.packb({"format": "xlc-pin-powers", "rows": power.tolist()})
    )
    print(f"k-eff {k:.6f} mesh {len(xs)} x {len(ys)}")


if __name__ == "__main__":
    main()
