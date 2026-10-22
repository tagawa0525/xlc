"""#0129 の試行の、表（table）と格子点ごと（exact）の差を出す。

使い方: python3 compare.py（このディレクトリで。make_inputs.py の入力を実行した出力を読む）

- 燃焼: 燃焼点ごとの k∞ の差（pcm、Δk/k）、2 群の少数群定数（全、吸収、νΣf、拡散係数、散乱）の相対差の最大、
  数密度の差（絶対値の最大を初期の数密度の最大で割ったもの。examples/depletion/verification.md と同じ）、時間の内訳
- 燃焼なし: 標準出力の k（6 桁）の差と、媒質の断面積の相対差の最大
"""

import pathlib

import msgpack

HERE = pathlib.Path(__file__).resolve().parent


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


def depletion(stem):
    t, e = (load(HERE / f"{stem}-{k}.msgpack") for k in ("table", "exact"))
    print(f"== {stem}")
    initial = max(max(d) for d in e["points"][0]["densities"])
    for pt, pe in zip(t["points"], e["points"], strict=True):
        dk = (pt["k_inf"] - pe["k_inf"]) / pe["k_inf"] * 1e5
        dens = max(
            abs(x - y)
            for rt, re_ in zip(pt["densities"], pe["densities"], strict=True)
            for x, y in zip(rt, re_, strict=True)
        )
        fg = ""
        if pt.get("few_group"):
            st, se = pt["few_group"]["states"][0], pe["few_group"]["states"][0]
            parts = {
                key: relative(st[key], se[key])
                for key in ("total", "absorption", "nu_fission", "diffusion", "scatter")
            }
            fg = " ".join(f"{k} {v:.1e}" for k, v in parts.items())
        print(
            f"burnup {pt['burnup']:5.2f} k {pe['k_inf']:.6f} dk {dk:+.3f} pcm "
            f"density {dens / initial:.1e} {fg}"
        )
    print(
        "seconds table",
        {k: round(v, 2) for k, v in t["seconds"].items()},
        "solves",
        t["solves"],
    )
    print(
        "seconds exact",
        {k: round(v, 2) for k, v in e["seconds"].items()},
        "solves",
        e["solves"],
    )


def static(stem):
    def k(kind):
        for line in (HERE / f"{stem}-{kind}.out").read_text().splitlines():
            if line.startswith("k "):
                return float(line.split()[1])
        return None

    kt, ke = k("table"), k("exact")
    mt, me = (load(HERE / f"{stem}-{x}.medium.msgpack") for x in ("table", "exact"))
    st, se = mt["states"][0], me["states"][0]
    worst = {
        key: relative(st[key], se[key])
        for key in st
        if isinstance(st[key], (list, float))
    }
    dk = (kt - ke) / ke * 1e5
    print(
        f"== {stem}: k {ke} dk {dk:+.2f} pcm (6 桁)",
        {k: f"{v:.1e}" for k, v in worst.items()},
    )


for stem in ("1C-5", "gd-fresh", "gd-mid"):
    static(stem)
for stem in ("1C-rings4", "2O-step"):
    depletion(stem)
