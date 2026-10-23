"""ADR-0195 の確認: 前の解から始めた出力を、平らな初期値からの出力と比べる。

使い方: python3 compare_warm_start.py history <平らな初期値からの履歴> <前の解からの履歴>
        python3 compare_warm_start.py table <平らな初期値からのテーブル> <前の解からのテーブル>
"""

import sys

import msgpack


def load(path):
    with open(path, "rb") as f:
        return msgpack.unpackb(f.read())


def leaves(value, prefix=""):
    """入れ子の値の数値の葉（名前と値）。"""
    if isinstance(value, dict):
        for k, v in value.items():
            yield from leaves(v, f"{prefix}.{k}")
    elif isinstance(value, list):
        for i, v in enumerate(value):
            yield from leaves(v, f"{prefix}[{i}]")
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield prefix, float(value)


def largest_relative(a, b, floor):
    """|a| か |b| が floor を超える葉の、相対差の最大と、その名前。"""
    worst = (0.0, None)
    for (name, x), (other, y) in zip(leaves(a), leaves(b), strict=True):
        if name != other:
            raise SystemExit(f"{name} vs {other}")
        scale = max(abs(x), abs(y))
        if scale > floor and abs(x - y) / scale > worst[0]:
            worst = (abs(x - y) / scale, name)
    return worst


def history(main, trial):
    main, trial = load(main), load(trial)
    if main["nuclides"] != trial["nuclides"]:
        raise SystemExit("the nuclides differ")
    for p, q in zip(main["points"], trial["points"]):
        if p["burnup"] != q["burnup"]:
            raise SystemExit(f"burnup {p['burnup']} vs {q['burnup']}")
        largest = max(max(abs(x) for x in n) for n in p["densities"])
        density = largest_relative(p["densities"], q["densities"], 1e-10 * largest)
        pins = largest_relative(p["pin_powers"], q["pin_powers"], 0.0)
        print(
            f"burnup {q['burnup']:6} k {p['k_inf']:.8f} {q['k_inf']:.8f} "
            f"dk {q['k_inf'] - p['k_inf']:+.2e} densities {density[0]:.2e} {density[1]} "
            f"pins {pins[0]:.2e}"
        )


def table(main, trial):
    main, trial = load(main), load(trial)
    names = [a["name"] for a in main["axes"]]
    if names != [a["name"] for a in trial["axes"]]:
        raise SystemExit(f"axes {names}")
    counts = [len(a["points"]) for a in trial["axes"]]
    for s, state in enumerate(trial["states"]):
        rest, at = s, []
        for n in reversed(counts):
            at.append(rest % n)
            rest //= n
        at.reverse()
        values = [a["points"][i] for a, i in zip(trial["axes"], at)]
        index = 0
        for axis, value in zip(main["axes"], values):
            index = index * len(axis["points"]) + axis["points"].index(value)
        other = main["states"][index]
        worst = largest_relative(other, state, 0.0)
        print(
            f"{dict(zip(names, values))} k {other['k_inf']:.8f} {state['k_inf']:.8f} "
            f"dk {state['k_inf'] - other['k_inf']:+.2e} largest {worst[0]:.2e} {worst[1]}"
        )


{"history": history, "table": table}[sys.argv[1]](*sys.argv[2:])
