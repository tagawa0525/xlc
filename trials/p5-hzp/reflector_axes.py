"""試行: 反射体の格子の分岐のテーブルで、状態の軸ごとに定数がどれだけ動くか（反射体のノードと燃料のノード）。

状態の並びは boron、moderator、fuel の順で fuel が最も速い（state = 4 i_boron + 2 i_moderator + i_fuel）。
軸ごとに、その軸だけが違う 4 組の相対差の最大を出す。

`--single <単一集合体のテーブル>` を前に置くと、続くテーブル（外周の燃料のノード）の燃料の温度の比（570 K / 560 K）と、
単一集合体の同じ比との差の最大も出す（分岐を単一集合体の比で近似したときの誤差）。

使い方: python3 reflector_axes.py [--single <テーブル>] <テーブル>...
"""

import sys

import msgpack

AXES = {"fuel": 1, "moderator": 2, "boron": 4}
KEYS = ["total", "absorption", "nu_fission", "diffusion", "scatter", "adf"]


def flat(x):
    if isinstance(x, list):
        return [v for item in x for v in flat(item)]
    return [] if x is None else [x]


def load(path: str) -> dict:
    return msgpack.unpackb(open(path, "rb").read(), raw=False, strict_map_key=False)


def ratio_difference(states: list, single: list) -> None:
    for key in KEYS:
        worst = 0.0
        for a in (0, 2, 4, 6):
            for x0, x1, y0, y1 in zip(
                flat(states[a][key]),
                flat(states[a + 1][key]),
                flat(single[a][key]),
                flat(single[a + 1][key]),
            ):
                if x0 and y0:
                    worst = max(worst, abs(x1 / x0 - y1 / y0))
        print(f"  fuel ratio vs single  {key} {100 * worst:.5f}%")


def main() -> None:
    args = sys.argv[1:]
    single = None
    if args[:1] == ["--single"]:
        single, args = load(args[1])["states"], args[2:]
    for path in args:
        t = load(path)
        names = [a["name"] for a in t["axes"]]
        if names != ["boron", "moderator", "fuel"] or len(t["states"]) != 8:
            raise SystemExit(f"{path}: axes {names}")
        states = t["states"]
        print(path)
        for axis, step in AXES.items():
            pairs = [(s, s + step) for s in range(8) if not s & step]
            worst = {}
            for key in KEYS:
                m = 0.0
                for a, b in pairs:
                    for x, y in zip(flat(states[a][key]), flat(states[b][key])):
                        if x != 0.0:
                            m = max(m, abs(y / x - 1))
                worst[key] = m
            print(
                f"  {axis:9s} "
                + "  ".join(f"{k} {100 * v:.4f}%" for k, v in worst.items())
            )
        if single is not None:
            ratio_difference(states, single)


if __name__ == "__main__":
    main()
