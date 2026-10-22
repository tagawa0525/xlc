"""#0045 の試行（make_inputs.py）の、間引いた格子（`t1e-3`、`t1e-4`）と間引かない格子（`full`）の差を出す（ADR-0185）。

使い方: python3 compare.py（このディレクトリで。make_inputs.py の入力を実行した出力を読む）

- 燃焼（2A-step、1C-rings1、1C-rings4）: 燃焼点ごとの k∞ の差（pcm、Δk/k）、数密度の差（絶対値の最大を初期の数密度の最大で
  割ったもの。examples/depletion/verification.md と同じ）、主要な核種の数密度の相対差の最大（領域の中の最大）、2 群の少数群定数の
  相対差の最大、時間の内訳（履歴の `seconds`）
- 1 状態（2A-60）: k∞ の差と 2 群の少数群定数の相対差の最大、自己遮蔽と全体の時間（ログ）
- 2A の燃焼（40 点）: 間引いた格子の履歴と examples/depletion/p2/lattice-2A.msgpack（間引かない格子）の差
- ピンごとの実効断面積（`pin_xs` の出力 `<入力>.pinxs.msgpack`。2A-60、2A-step）: 核種ごとの UFG の範囲（4 eV〜20 keV）の捕獲
  （全 − 弾性散乱 − 核分裂）と核分裂の、群のレサジー幅の重み（1/E）の和の相対差の、ピンの中の最大
"""

import itertools
import math
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


def resonance_weights(library):
    """ライブラリの群のうち UFG の範囲（4 eV〜20 keV）の群のレサジー幅（ほかの群は 0）"""
    bounds = load(library)["group_bounds"]
    return [
        math.log(hi / lo) if lo >= 4.0 and hi <= 2.0e4 else 0.0
        for hi, lo in itertools.pairwise(bounds)
    ]


def pin_rates(path, weights):
    """{核種: [(捕獲, 核分裂)]}（燃料のあるピンごとの、UFG の範囲の 1/E の重みの和）"""
    out = {}
    for name, x in load(path)["nuclides"].items():
        pins = []
        for t_row, e_row, f_row in zip(
            x["total"], x["elastic"], x["fission"], strict=True
        ):
            for t, e, f in zip(t_row, e_row, f_row, strict=True):
                if any(t):
                    capture = sum(
                        w * (a - b - c)
                        for w, a, b, c in zip(weights, t, e, f, strict=True)
                    )
                    fission = sum(w * c for w, c in zip(weights, f, strict=True))
                    pins.append((capture, fission))
        out[name] = pins
    return out


def captures(directory, base, others, library, top=12):
    weights = resonance_weights(library)
    b = pin_rates(directory / f"{base}.pinxs.msgpack", weights)
    for other in others:
        o = pin_rates(directory / f"{other}.pinxs.msgpack", weights)
        rows = []
        for name, pins in b.items():
            worst = [0.0, 0.0]
            for (cb, fb), (co, fo) in zip(pins, o[name], strict=True):
                for q, (x, y) in enumerate(((co, cb), (fo, fb))):
                    if y != 0:
                        worst[q] = max(worst[q], abs(x / y - 1))
            rows.append((worst[0], worst[1], name))
        rows.sort(reverse=True)
        print(
            f"-- {other} / {base}: UFG の範囲の捕獲と核分裂の相対差の、ピンの中の最大（大きい順）"
        )
        print("   " + "、".join(f"{n} {c:.1e}/{f:.1e}" for c, f, n in rows[:top]))


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
    library = HERE.parent / "p2" / "depletion.msgpack"
    for stem in ("2A-60", "2A-step"):
        if (HERE / f"{stem}-full.pinxs.msgpack").exists():
            captures(HERE, f"{stem}-full", [f"{stem}-{v}" for v in VARIANTS], library)
    depletion(
        "2A-depletion",
        HERE.parent / "p2" / "lattice-2A.msgpack",
        [v for v in VARIANTS if (HERE / f"2A-depletion-{v}.msgpack").exists()],
    )
