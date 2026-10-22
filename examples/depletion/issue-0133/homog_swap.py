"""#0133 の試行: 均質な無限媒質の xlc の balance（`homogeneous` の例の `--balance`）から多群の固有値問題を組み直し、核種の
群の範囲の P0 の移行の行、核種の正味の吸収と生成（反応率 / 中性子束）、核種の χ を OpenMC のもの（`openmc_homogeneous.py` の
balance）や別のライブラリのものに置き換えて、k、群 1〜13 の中性子束（生成 1 あたり）、U-238 と U-235 の生成（生成 1 あたり）の
xlc/OpenMC − 1 を出す。

使い方: python3 homog_swap.py <xlc の balance> <OpenMC の balance> <xlc のライブラリ> <温度> [置き換え ...]

置き換え:
- `<核種>:<最初の群>-<最後の群>`: 散乱による移行の行を OpenMC の移行の率 / 中性子束に（全断面積は行の和の変化だけ動かす）
- `xs:<核種>:<最初の群>-<最後の群>`: 正味の吸収と生成を OpenMC の率 / 中性子束に
- `chi:<msgpack>`: 核種の名前から χ（群の値の列）への辞書の χ に

群は 1 始まりで高いエネルギーから数える。OpenMC の生成は群構造の外（10 MeV より上と下）を含めた和で規格化する。xlc の balance は
10 MeV より上に生まれた中性子を群 1 に含むが、OpenMC の balance の群 1 は含まない（群 1 の中性子束の差はそのぶん大きく出る）。
"""

import sys
from pathlib import Path

import msgpack
import numpy as np

FAST = 13


def read(path):
    return msgpack.unpackb(Path(path).read_bytes(), strict_map_key=False)


def main():
    xb, ob, lib = (read(p) for p in sys.argv[1:4])
    temperature = float(sys.argv[4])
    xlc, omc = xb["materials"], ob["materials"]
    names = list(xlc)
    phi = np.array(xlc[names[0]]["flux"])
    ng = len(phi)
    chis = {}
    for n in lib["nuclides"]:
        if n["name"] in xlc:
            t = next(t for t in n["temperatures"] if t["temperature"] == temperature)
            chis[n["name"]] = np.array(t["chi"])
    transfer, cross, chi_over = [], [], {}
    for spec in sys.argv[5:]:
        head, tail = spec.rsplit(":", 1)
        if head == "chi":
            chi_over = {k: np.array(v) for k, v in read(tail).items()}
            continue
        lo, hi = map(int, tail.split("-"))
        if head.startswith("xs:"):
            cross.append((head[3:], lo, hi))
        else:
            transfer.append((head, lo, hi))

    for name, _, _ in transfer + cross:
        if name not in xlc or name not in omc:
            raise SystemExit(f"{name}: not in both balances")

    def operator(swap):
        total, scatter, nu_fission = np.zeros(ng), np.zeros((ng, ng)), {}
        for name in names:
            r = xlc[name]
            t = np.array(r["transfer"]) / phi[:, None]
            a = np.array(r["absorption"]) / phi
            nf = np.array(r["production"]) / phi
            if swap:
                o = omc.get(name)
                for n, lo, hi in transfer:
                    if n == name:
                        for g in range(lo - 1, hi):
                            t[g] = np.array(o["transfer"][g]) / o["flux"][g]
                for n, lo, hi in cross:
                    if n == name:
                        for g in range(lo - 1, hi):
                            a[g] = o["absorption"][g] / o["flux"][g]
                            nf[g] = o["production"][g] / o["flux"][g]
            total += a + t.sum(axis=1)
            scatter += t
            nu_fission[name] = nf
        return np.diag(total) - scatter.T, nu_fission

    def solve(swap):
        a, nu_fission = operator(swap)
        chi = {**chis, **chi_over} if swap else chis
        f, k = phi.copy(), 1.0
        for _ in range(2000):
            rates = {n: (nu_fission[n] * f).sum() for n in nu_fission}
            source = sum(chi[n] * rates[n] for n in nu_fission if rates[n] > 0.0)
            f_new = np.linalg.solve(a, source / k)
            produced = sum((nu_fission[n] * f_new).sum() for n in nu_fission)
            k_new = k * produced / sum(rates.values())
            f = f_new / produced
            done = abs(k_new - k) < 1e-12
            k = k_new
            if done:
                break
        return k, f, {n: nu_fission[n] * f for n in nu_fission}

    k0, _, _ = solve(False)
    k, f, production = solve(True)
    outside = ob["outside"]
    produced = sum(np.sum(m["production"]) for m in omc.values())
    produced += outside["above"]["production"] + outside["below"]["production"]
    flux_o = np.array(omc[names[0]]["flux"]) / produced
    print(
        f"k: xlc {k0:.6f} -> {k:.6f} ({(k - k0) * 1e5:+.1f} pcm), OpenMC {ob['k']:.6f}"
    )
    print(
        f"flux per production xlc/OpenMC - 1 (%), groups 1-{FAST}: "
        + " ".join(f"{(f[g] / flux_o[g] - 1) * 100:+.2f}" for g in range(FAST))
    )
    for n in ("U-238", "U-235"):
        p_o = np.array(omc[n]["production"]) / produced
        print(
            f"{n} production per production: groups 1-{FAST} "
            f"{(production[n][:FAST].sum() / p_o[:FAST].sum() - 1) * 100:+.3f}%, "
            f"all {(production[n].sum() / p_o.sum() - 1) * 100:+.3f}%"
        )


main()
