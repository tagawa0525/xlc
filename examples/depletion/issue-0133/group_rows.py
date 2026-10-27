"""#0152 の試行: 物質の群の移行の行を、OpenMC の群の中の細かい入射の区間ごとの行（核）と、それを群にまとめる重み（群の中の
中性子束）に分け、重みを入れ替えた行を作る。xlc の行がどれだけ群の中のスペクトルの違いで OpenMC の行と違うかを測るため。

使い方:
    python3 group_rows.py rows <状態の transfer.msgpack> <核の transfer.msgpack> <xlc の balance.msgpack> <xlc の fine.msgpack>
        <物質> <最初の群> <最後の群>
    python3 group_rows.py write <状態の transfer.msgpack> <核の transfer.msgpack> <xlc の fine.msgpack> <物質> <最初の群>
        <最後の群> <重み> <出力>
    python3 group_rows.py spectrum <核の transfer.msgpack> <xlc の fine.msgpack> <物質> <下端 eV> <上端 eV> <区間の数>

<状態の transfer.msgpack> は openmc_pin_transfer.py の出力（群の行）、<核の transfer.msgpack> は細かい入射の区間を数えた同じ
スクリプトの出力（`fine`）、<xlc の balance.msgpack> は `balance` の例の出力、<xlc の fine.msgpack> は `fine_flux` の例の出力
（等価なセルの燃料と減速材のリングの格子点ごとのレサジーあたりの中性子束）。物質は xlc の名前（`water` か `fuel`）。群は 1 始まりで、
高いエネルギーから数える。

区間 i の核は r_i(h) = T_i(h) / φ_i（T_i は区間から群 h への移行、φ_i は区間の中性子束。OpenMC の値）。重み w の群の行は
Σ_i w_i r_i(h) / Σ_i w_i。重みは `openmc`（核の実行の φ_i。群の行の値に戻る）、`xlc`（xlc の中性子束の区間の積分）、`lethargy`
（区間のレサジーの幅。1/E の中性子束）。

rows は、群ごとに、群から出る移行（自分の群への移行を除く行の和）と、自分の群への移行（1/cm）を、xlc、状態の OpenMC、核の 3 つの重みで
書く。write は、状態の transfer.msgpack の物質の行のうち、群の範囲の行を、核を重みで群にまとめた行（の中性子束倍）に置き換えて書く
（`swap_transfer` の例に渡すため。ほかの物質と群はそのまま）。spectrum は、群の範囲 [下端, 上端] を同じ数の区間に分け、
区間ごとの OpenMC と xlc の中性子束（範囲の中で規格化）の比 − 1 を書く。
"""

import itertools
import sys
from pathlib import Path

import msgpack
import numpy as np

# OpenMC の物質の xlc の名前から、fine_flux の出力のリングの名前へ
FINE_RING = {"water": "moderator", "fuel": "fuel"}
WEIGHTS = ("openmc", "xlc", "lethargy")


def read(path):
    return msgpack.unpackb(Path(path).read_bytes())


class Kernel:
    """核の実行の細かい入射の区間（高いエネルギーから）"""

    def __init__(self, data, material):
        fine = data.get("fine")
        if fine is None:
            raise SystemExit("the kernel run has no fine tallies")
        self.bounds = np.array(fine["bounds"])
        entry = fine["materials"][material]
        if "transfer_mean" not in entry:
            raise SystemExit(f"the kernel run has no fine transfer of `{material}`")
        self.flux = np.array(entry["mean"])
        self.transfer = np.array(entry["transfer_mean"])
        self.group_bounds = np.array(data["bounds"])

    def bins(self, g):
        """群 g（0 始まり）の区間の番号"""
        hi, lo = self.group_bounds[g], self.group_bounds[g + 1]
        inside = np.where(
            (self.bounds[:-1] <= hi * (1 + 1e-9)) & (self.bounds[1:] >= lo * (1 - 1e-9))
        )[0]
        if len(inside) == 0 or not (
            np.isclose(self.bounds[inside[0]], hi, rtol=1e-9)
            and np.isclose(self.bounds[inside[-1] + 1], lo, rtol=1e-9)
        ):
            raise SystemExit(f"group {g + 1} is not divided into fine bins")
        return inside

    def row(self, g, weight):
        """群 g の行（1/cm）を、区間の重み `weight` でまとめたもの"""
        i = self.bins(g)
        if np.any(self.flux[i] <= 0.0):
            raise SystemExit(f"group {g + 1}: a fine bin has no flux")
        kernel = self.transfer[i] / self.flux[i][:, None]
        w = weight[i] if weight is not None else self.flux[i]
        return (w[:, None] * kernel).sum(0) / w.sum()


def xlc_weight(fine, ring, bounds):
    """xlc の格子点ごとのレサジーあたりの中性子束を、区間 [bounds[j+1], bounds[j]] で積分したもの（点の間は線形。fine_compare.py と同じ）"""
    e = np.array(fine["energies"])[::-1]
    f = np.array(fine[ring])[::-1]
    u = np.log(e)
    out = np.zeros(len(bounds) - 1)
    for j in range(len(bounds) - 1):
        ua, uc = np.log(bounds[j + 1]), np.log(bounds[j])
        if ua < u[0] - 1e-12 or uc > u[-1] + 1e-12:
            continue
        inner = (u > ua) & (u < uc)
        grid = np.concatenate([[ua], u[inner], [uc]])
        values = np.concatenate(
            [[np.interp(ua, u, f)], f[inner], [np.interp(uc, u, f)]]
        )
        out[j] = np.trapezoid(values, grid)
    return out


def weights(kernel, fine, material):
    return {
        "openmc": None,
        "xlc": xlc_weight(fine, FINE_RING[material], kernel.bounds),
        "lethargy": np.log(kernel.bounds[:-1] / kernel.bounds[1:]),
    }


def group_row(data, material, g):
    m = data["materials"][material]
    return np.array(m["transfer"][g]) / m["flux"][g]


def groups(first, last):
    first, last = int(first), int(last)
    if first < 1 or last < first:
        raise SystemExit(f"groups {first}..{last}: count from 1")
    return range(first - 1, last)


def rows(state, kernel_path, balance, fine, material, first, last):
    state, xlc = read(state), read(balance)
    kernel = Kernel(read(kernel_path), material)
    w = weights(kernel, read(fine), material)

    def split(row, g):
        return row.sum() - row[g], row[g]

    for g in groups(first, last):
        x = (
            np.array(xlc["materials"][material]["transfer"][g])
            / (xlc["materials"][material]["flux"][g])
        )
        values = {"xlc": split(x, g), "state": split(group_row(state, material, g), g)}
        for name in WEIGHTS:
            values[f"kernel-{name}"] = split(kernel.row(g, w[name]), g)
        text = " ".join(f"{k} {out:.6e} {own:.6e}" for k, (out, own) in values.items())
        print(f"group {g + 1} out/own: {text}")


def write(state_path, kernel_path, fine, material, first, last, weight, output):
    if weight not in WEIGHTS:
        raise SystemExit(f"weight must be one of {WEIGHTS}")
    state = read(state_path)
    kernel = Kernel(read(kernel_path), material)
    w = weights(kernel, read(fine), material)[weight]
    m = state["materials"][material]
    for g in groups(first, last):
        m["transfer"][g] = (kernel.row(g, w) * m["flux"][g]).tolist()
    Path(output).write_bytes(msgpack.packb(state))


def spectrum(kernel_path, fine, material, lo, hi, n):
    kernel = Kernel(read(kernel_path), material)
    lo, hi, n = float(lo), float(hi), int(n)
    b = kernel.bounds
    inside = np.where((b[:-1] <= hi * (1 + 1e-9)) & (b[1:] >= lo * (1 - 1e-9)))[0]
    po = kernel.flux[inside]
    px = xlc_weight(read(fine), FINE_RING[material], b)[inside]
    po, px = po / po.sum(), px / px.sum()
    edges = np.linspace(0, len(inside), n + 1).astype(int)
    for a, c in itertools.pairwise(edges):
        print(
            f"{b[inside[c - 1] + 1]:.4g}-{b[inside[a]]:.4g} eV: xlc {px[a:c].sum():.5f} "
            f"OpenMC {po[a:c].sum():.5f} 比 {px[a:c].sum() / po[a:c].sum() - 1:+.4f}"
        )


COMMANDS = {"rows": (rows, 7), "write": (write, 8), "spectrum": (spectrum, 6)}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        raise SystemExit(__doc__)
    command, count = COMMANDS[sys.argv[1]]
    if len(sys.argv) != count + 2:
        raise SystemExit(__doc__)
    command(*sys.argv[2:])
