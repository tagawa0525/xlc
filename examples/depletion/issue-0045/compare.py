"""#0045 の試行（make_inputs.py）の、間引いた格子（`t1e-3`、`t1e-4`）と間引かない格子（`full`）の差を出す（ADR-0185）。

使い方: python3 compare.py（このディレクトリで。make_inputs.py の入力を実行した出力を読む）

- 燃焼（2A-step、1C-rings1、1C-rings4）: 燃焼点ごとの k∞ の差（pcm、Δk/k）、数密度の差（絶対値の最大を初期の数密度の最大で
  割ったもの。examples/depletion/verification.md と同じ）、主要な核種の数密度の相対差の最大（領域の中の最大）、2 群の少数群定数の
  相対差の最大、時間の内訳（履歴の `seconds`）
- 1 状態（2A-60）: k∞ の差と 2 群の少数群定数の相対差の最大、自己遮蔽と全体の時間（ログ）
- 2A の燃焼（40 点）: 間引いた格子の履歴と examples/depletion/p2/lattice-2A.msgpack（間引かない格子）の差
"""

import pathlib
import re

import msgpack

HERE = pathlib.Path(__file__).resolve().parent
VARIANTS = ("t1e-3", "t1e-4")
MAJOR = ("U-235", "U-238", "Pu-239", "Pu-240", "Pu-241", "Xe-135", "Sm-149")
FIELDS = ("total", "absorption", "nu_fission", "diffusion", "scatter")


def load(path):
    return msgpack.unpackb(path.read_bytes())


def relative(a, b):
    worst = 0.0
    if isinstance(a, dict):
        for k in a:
            if isinstance(a[k], (list, dict, float)):
                worst = max(worst, relative(a[k], b[k]))
        return worst
    if isinstance(a, list):
        for x, y in zip(a, b, strict=True):
            worst = max(worst, relative(x, y))
        return worst
    if a is None or b is None:
        return 0.0
    return abs(a - b) if b == 0 else abs(a / b - 1)


def few_group(a, b):
    return " ".join(f"{k} {relative(a[k], b[k]):.1e}" for k in FIELDS)


def log_times(stem):
    text = (HERE / f"{stem}.log").read_text()
    shielding = [
        float(m.group(1)) * {"s": 1.0, "ms": 1e-3}[m.group(2)]
        for m in re.finditer(r"shielding: close time\.busy=([\d.]+)(ms|s) ", text)
    ]
    points = re.findall(r"points_sum=(\d+)", text)
    return shielding, points


def depletion(stem, base_path=None, variants=VARIANTS):
    base = load(base_path or HERE / f"{stem}-full.msgpack")
    names = base["nuclides"]
    major = [names.index(n) for n in MAJOR if n in names]
    initial = max(max(d) for d in base["points"][0]["densities"])
    print(f"== {stem}: full の時間 {fmt(base['seconds'])}")
    for lst in variants:
        other = load(HERE / f"{stem}-{lst}.msgpack")
        print(f"-- {lst}: 時間 {fmt(other['seconds'])}")
        for po, pb in zip(other["points"], base["points"], strict=True):
            dk = (po["k_inf"] - pb["k_inf"]) / pb["k_inf"] * 1e5
            dens = max(
                abs(x - y)
                for ro, rb in zip(po["densities"], pb["densities"], strict=True)
                for x, y in zip(ro, rb, strict=True)
            )
            rel = {
                names[i]: max(
                    abs(ro[i] / rb[i] - 1)
                    for ro, rb in zip(po["densities"], pb["densities"], strict=True)
                    if rb[i] > 0
                )
                for i in major
                if any(rb[i] > 0 for rb in pb["densities"])
            }
            fg = ""
            if po.get("few_group"):
                fg = few_group(
                    po["few_group"]["states"][0], pb["few_group"]["states"][0]
                )
            worst = max(rel, key=rel.get) if rel else "-"
            print(
                f"burnup {po['burnup']:5.2f} k {pb['k_inf']:.6f} dk {dk:+.2f} pcm "
                f"density {dens / initial:.1e} major {worst} {rel.get(worst, 0):.1e} {fg}"
            )


def fmt(seconds):
    return {k: round(v, 1) for k, v in seconds.items()}


def state(path):
    """少数群定数のファイル（分岐のテーブルの形。1 状態）の状態"""
    return load(path)["states"][0]


def static(stem):
    base = state(HERE / f"{stem}-full.fg.msgpack")
    shielding, points = log_times(f"{stem}-full")
    print(
        f"== {stem}: full k {base['k_inf']:.6f} 自己遮蔽 {shielding} s、UFG の点 {points}"
    )
    for lst in VARIANTS:
        other = state(HERE / f"{stem}-{lst}.fg.msgpack")
        dk = (other["k_inf"] - base["k_inf"]) / base["k_inf"] * 1e5
        shielding, points = log_times(f"{stem}-{lst}")
        print(
            f"-- {lst}: dk {dk:+.2f} pcm {few_group(other, base)} 自己遮蔽 {shielding} s、UFG の点 {points}"
        )


if __name__ == "__main__":
    for stem in ("2A-step", "1C-rings1", "1C-rings4"):
        depletion(stem)
    static("2A-60")
    depletion(
        "2A-depletion",
        HERE.parent / "p2" / "lattice-2A.msgpack",
        [v for v in VARIANTS if (HERE / f"2A-depletion-{v}.msgpack").exists()],
    )
