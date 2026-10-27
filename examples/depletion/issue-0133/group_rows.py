"""#0152 の試行: 物質の群の移行の行を、OpenMC の群の中の細かい入射の区間ごとの行（核）と、それを群にまとめる重み（群の中の
中性子束）に分け、重みを入れ替えた行を作る。xlc の行がどれだけ群の中のスペクトルの違いで OpenMC の行と違うかを測るため。

使い方:
    python3 group_rows.py rows <状態の transfer.msgpack> <核の transfer.msgpack> <xlc の balance.msgpack> <xlc の fine.msgpack>
        <物質> <最初の群> <最後の群>
    python3 group_rows.py write <状態の transfer.msgpack> <核の transfer.msgpack> <xlc の fine.msgpack> <物質> <最初の群>
        <最後の群> <重み> <出力>
    python3 group_rows.py spectrum <核の transfer.msgpack> <xlc の fine.msgpack> <物質> <下端 eV> <上端 eV> <区間の数>
    python3 group_rows.py shift <xlc の fine.msgpack (a)> <xlc の fine.msgpack (b)> <リング> <下端 eV> <上端 eV> <区間の数>
    python3 group_rows.py summarize <設定 (TOML)> <JSON>

<状態の transfer.msgpack> は openmc_pin_transfer.py の出力（群の行）、<核の transfer.msgpack> は細かい入射の区間を数えた同じ
スクリプトの出力（`fine`）、<xlc の balance.msgpack> は `balance` の例の出力、<xlc の fine.msgpack> は `fine_flux` の例の出力
（等価なセルの燃料と減速材のリングの格子点ごとのレサジーあたりの中性子束）。物質は xlc の名前（`water` か `fuel`）。群は 1 始まりで、
高いエネルギーから数える。

区間 i の核は r_i(h) = T_i(h) / φ_i（T_i は区間から群 h への移行、φ_i は区間の中性子束。OpenMC の値）。重み w の群の行は
Σ_i w_i r_i(h) / Σ_i w_i。重みは `openmc`（核の実行の φ_i。群の行の値に戻る）、`xlc`（xlc の中性子束の区間の積分）、`lethargy`
（区間のレサジーの幅。1/E の中性子束）。

shift は xlc の 2 つの状態の、リング（`moderator` か `fuel`）の範囲の中性子束を、レサジーの等間隔の区間で積分し、範囲の中で規格化した比
a / b − 1 を書く（状態で群の中のスペクトルがどう変わるかを見るため）。

summarize は group-rows/run.sh のログ（設定のディレクトリの `logs/`）を表の JSON にまとめる。設定は `[[spec]]`（`spec` は
`<物質>:<最初の群>-<最後の群>`、`label`）、`[[kernel]]`（`key`、`label`）、`[[state]]`（`key`、`name`、`kernels` は核の名前の並び、`specs` はその状態で流した置き換え。省くとすべての `[[spec]]`。ログの置き換えが
これと違えば止める）。
寄与は swap_transfer の SPH 因子を置き換える前のままの k の変化の符号を変えたもの（pcm。swap_candidates.py と同じ）で、`openmc` は状態の
OpenMC の行、`kernels` は核ごとの `xlc` と `lethargy` の重みの行に置き換えたもの。`logs/<状態>-spectrum.txt`（spectrum の出力に
`material <物質>` を前に置いたもの）があれば、区間ごとの比 − 1 を低いエネルギーから書く。

rows は、群ごとに、群から出る移行（自分の群への移行を除く行の和）と、自分の群への移行（1/cm）を、xlc、状態の OpenMC、核の 3 つの重みで
書く。write は、状態の transfer.msgpack の物質の行のうち、群の範囲の行を、核を重みで群にまとめた行（の中性子束倍）に置き換えて書く
（`swap_transfer` の例に渡すため。ほかの物質と群はそのまま）。spectrum は、群の範囲 [下端, 上端] を同じ数の区間に分け、
区間ごとの OpenMC と xlc の中性子束（範囲の中で規格化）の比 − 1 を書く。
"""

import itertools
import json
import re
import sys
from pathlib import Path

import msgpack
import numpy as np
import tomllib
from swap_candidates import SummaryError as CandidatesError
from swap_candidates import transfer_runs

# OpenMC の物質の xlc の名前から、fine_flux の出力のリングの名前へ
FINE_RING = {"water": "moderator", "fuel": "fuel"}
WEIGHTS = ("openmc", "xlc", "lethargy")


class SummaryError(Exception):
    pass


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

    def range_bins(self, lo, hi):
        """[lo, hi] eV をちょうど埋める区間の番号（端が区間の境でなければ止める）"""
        inside = np.where(
            (self.bounds[:-1] <= hi * (1 + 1e-9)) & (self.bounds[1:] >= lo * (1 - 1e-9))
        )[0]
        if len(inside) == 0 or not (
            np.isclose(self.bounds[inside[0]], hi, rtol=1e-9)
            and np.isclose(self.bounds[inside[-1] + 1], lo, rtol=1e-9)
        ):
            raise SystemExit(f"{lo}-{hi} eV is not divided into fine bins")
        return inside

    def bins(self, g):
        """群 g（0 始まり）の区間の番号"""
        return self.range_bins(self.group_bounds[g + 1], self.group_bounds[g])

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
            raise SystemExit(
                f"{bounds[j + 1]}-{bounds[j]} eV is outside the xlc grid {e[0]}-{e[-1]} eV"
            )
        inner = (u > ua) & (u < uc)
        grid = np.concatenate([[ua], u[inner], [uc]])
        values = np.concatenate(
            [[np.interp(ua, u, f)], f[inner], [np.interp(uc, u, f)]]
        )
        out[j] = np.trapezoid(values, grid)
    return out


def weight_of(name, kernel, fine, material, bins):
    """区間 `bins`（kernel の区間の番号）の重み。ほかの区間は 0（`Kernel.row` は群の区間だけを使う）"""
    if name == "openmc":
        return None
    b = kernel.bounds
    w = np.zeros(len(b) - 1)
    if name == "xlc":
        edges = np.concatenate([b[bins], [b[bins[-1] + 1]]])
        w[bins] = xlc_weight(fine, FINE_RING[material], edges)
    elif name == "lethargy":
        w[bins] = np.log(b[bins] / b[bins + 1])
    else:
        raise SystemExit(f"weight must be one of {WEIGHTS}")
    return w


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
    fine = read(fine)

    def split(row, g):
        return row.sum() - row[g], row[g]

    for g in groups(first, last):
        values = {
            "xlc": split(group_row(xlc, material, g), g),
            "state": split(group_row(state, material, g), g),
        }
        for name in WEIGHTS:
            w = weight_of(name, kernel, fine, material, kernel.bins(g))
            values[f"kernel-{name}"] = split(kernel.row(g, w), g)
        text = " ".join(f"{k} {out:.6e} {own:.6e}" for k, (out, own) in values.items())
        print(f"group {g + 1} out/own: {text}")


def write(state_path, kernel_path, fine, material, first, last, weight, output):
    if weight not in WEIGHTS:
        raise SystemExit(f"weight must be one of {WEIGHTS}")
    state = read(state_path)
    kernel = Kernel(read(kernel_path), material)
    fine = read(fine)
    m = state["materials"][material]
    for g in groups(first, last):
        w = weight_of(weight, kernel, fine, material, kernel.bins(g))
        m["transfer"][g] = (kernel.row(g, w) * m["flux"][g]).tolist()
    Path(output).write_bytes(msgpack.packb(state))


def spectrum(kernel_path, fine, material, lo, hi, n):
    kernel = Kernel(read(kernel_path), material)
    lo, hi, n = float(lo), float(hi), int(n)
    b = kernel.bounds
    inside = kernel.range_bins(lo, hi)
    if not 1 <= n <= len(inside):
        raise SystemExit(f"{n} intervals for {len(inside)} fine bins")
    po = kernel.flux[inside]
    px = weight_of("xlc", kernel, read(fine), material, inside)[inside]
    po, px = po / po.sum(), px / px.sum()
    edges = np.linspace(0, len(inside), n + 1).astype(int)
    for a, c in itertools.pairwise(edges):
        print(
            f"{b[inside[c - 1] + 1]:.4g}-{b[inside[a]]:.4g} eV: xlc {px[a:c].sum():.5f} "
            f"OpenMC {po[a:c].sum():.5f} 比 {px[a:c].sum() / po[a:c].sum() - 1:+.4f}"
        )


def shift_lines(fine_a, fine_b, ring, lo, hi, n):
    """xlc の 2 つの状態の、リング `ring` の [lo, hi] eV の中性子束を、レサジーの等間隔の n 区間で積分し、範囲の中で規格化した
    比 a / b − 1 の行（低いエネルギーから）"""
    edges = np.geomspace(hi, lo, n + 1)
    a = xlc_weight(fine_a, ring, edges)[::-1]
    b = xlc_weight(fine_b, ring, edges)[::-1]
    a, b = a / a.sum(), b / b.sum()
    up = edges[::-1]
    return [
        f"{up[j]:.4g}-{up[j + 1]:.4g} eV: 比 {a[j] / b[j] - 1:+.4f}" for j in range(n)
    ]


def shift(fine_a, fine_b, ring, lo, hi, n):
    for line in shift_lines(
        read(fine_a), read(fine_b), ring, float(lo), float(hi), int(n)
    ):
        print(line)


def parse_rows(text):
    """rows の出力（核ごとに `kernel <名前> material <物質> <群>` を前に置いたもの）から、核ごと、`<物質>:<群>` ごとの重みごとの
    [群から出る移行, 自分の群への移行]"""
    found = {}
    kernel = material = None
    for line in text.splitlines():
        if line.startswith("kernel "):
            m = re.fullmatch(r"kernel (\S+) material (\S+) \d+-\d+", line)
            if m is None:
                raise SummaryError(f"cannot read the kernel line: {line}")
            kernel, material = m.groups()
        elif line.startswith("group "):
            if kernel is None:
                raise SummaryError(f"a group line before its kernel: {line}")
            m = re.fullmatch(r"group (\d+) out/own: (.+)", line)
            if m is None:
                raise SummaryError(f"cannot read the group line: {line}")
            words = m.group(2).split()
            if len(words) % 3:
                raise SummaryError(f"cannot read the group line: {line}")
            values = {
                words[i]: [float(words[i + 1]), float(words[i + 2])]
                for i in range(0, len(words), 3)
            }
            entries = found.setdefault(kernel, {})
            name = f"{material}:{m.group(1)}"
            if name in entries:
                raise SummaryError(f"the rows log repeats {kernel} {name}")
            entries[name] = values
        elif line.strip():
            raise SummaryError(f"cannot read the rows line: {line}")
    return found


# spectrum の区間の境（`:.4g`。指数の表記を含む）
NUMBER = r"\d+(?:\.\d+)?(?:e[+-]\d+)?"


def parse_spectrum(text):
    """spectrum の出力（物質ごとに `material <物質>` を前に置いたもの）から、物質ごとの区間の下端、上端 (eV) と比 − 1（低いエネルギーから）"""
    found = {}
    material = None
    for line in text.splitlines():
        if line.startswith("material "):
            material = line.split()[1]
            if material in found:
                raise SummaryError(f"the spectrum log repeats `{material}`")
            found[material] = []
        elif line.strip():
            m = re.fullmatch(
                rf"({NUMBER})-({NUMBER}) eV: xlc [\d.]+ OpenMC [\d.]+ 比 ([+-][\d.]+)",
                line,
            )
            if m is None or material is None:
                raise SummaryError(f"cannot read the spectrum line: {line}")
            found[material].append(
                {
                    "low": float(m.group(1)),
                    "high": float(m.group(2)),
                    "ratio": float(m.group(3)),
                }
            )
    for bins in found.values():
        bins.sort(key=lambda b: b["low"])
    return found


def contributions(path, k=None, openmc_k=None):
    """swap_transfer のログから、置き換えない k、OpenMC の k と、`<物質>:<最初の群>-<最後の群>` ごとの寄与（pcm）。`k` があれば、
    置き換えない k が同じこと（同じ入力とコードのログ）を確かめる"""
    try:
        runs = transfer_runs(Path(path).read_text())
    except CandidatesError as e:
        raise SummaryError(f"{path}: {e}") from e
    if not runs:
        raise SummaryError(f"{path}: no swap")
    bases = {r["base"] for r in runs.values()} | ({k} if k is not None else set())
    if len(bases) != 1:
        raise SummaryError(f"{path}: the unswapped k differs: {sorted(bases)}")
    references = {r["openmc"] for r in runs.values()} | (
        {openmc_k} if openmc_k is not None else set()
    )
    if len(references) != 1:
        raise SummaryError(f"{path}: the OpenMC k differs: {sorted(references)}")
    (base,) = bases
    (openmc,) = references
    found = {f"{m}:{a}-{b}": -r["change"] for (m, a, b), r in runs.items()}
    return base, openmc, found


def summarize(config):
    config = Path(config)
    c = tomllib.loads(config.read_text())
    logs = config.parent / "logs"
    specs = [s["spec"] for s in c["spec"]]
    states = []
    for state in c["state"]:
        key = state["key"]
        k, openmc_k, openmc = contributions(logs / f"{key}-transfer.txt")
        wanted = state.get("specs", specs)
        if set(openmc) != set(wanted):
            raise SummaryError(f"{key}: the OpenMC row swaps are {sorted(openmc)}")
        swaps = {spec: {"openmc": openmc[spec], "kernels": {}} for spec in wanted}
        for kernel in state["kernels"]:
            for weight in ("xlc", "lethargy"):
                _, _, found = contributions(
                    logs / f"{key}-{kernel}-{weight}-transfer.txt", k, openmc_k
                )
                if set(found) != set(wanted):
                    raise SummaryError(
                        f"{key}: the {kernel} {weight} swaps are {sorted(found)}"
                    )
                for spec, pcm in found.items():
                    swaps[spec]["kernels"].setdefault(kernel, {})[weight] = pcm
        rows = parse_rows((logs / f"{key}-rows.txt").read_text())
        if set(rows) != set(state["kernels"]):
            raise SummaryError(f"{key}: the rows log has kernels {sorted(rows)}")
        expected = {
            f"{spec.split(':')[0]}:{g + 1}"
            for spec in wanted
            for g in groups(*spec.split(":")[1].split("-"))
        }
        for kernel, entries in rows.items():
            if set(entries) != expected:
                raise SummaryError(
                    f"{key}: the {kernel} rows are {sorted(entries)}, not {sorted(expected)}"
                )
        spectrum = logs / f"{key}-spectrum.txt"
        states.append(
            {
                "key": key,
                "name": state["name"],
                "k": k,
                "openmc_k": openmc_k,
                "kernels": state["kernels"],
                "swaps": swaps,
                "rows": rows,
                "spectrum": parse_spectrum(spectrum.read_text())
                if spectrum.exists()
                else None,
            }
        )
    return {
        "format": "xlc-group-rows",
        "specs": c["spec"],
        "kernels": c["kernel"],
        "states": states,
    }


def write_summary(config, output):
    Path(output).write_text(
        json.dumps(summarize(config), ensure_ascii=False, indent=1) + "\n"
    )


COMMANDS = {
    "rows": (rows, 7),
    "write": (write, 8),
    "spectrum": (spectrum, 6),
    "summarize": (write_summary, 2),
    "shift": (shift, 6),
}

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        raise SystemExit(__doc__)
    command, count = COMMANDS[sys.argv[1]]
    if len(sys.argv) != count + 2:
        raise SystemExit(__doc__)
    command(*sys.argv[2:])
