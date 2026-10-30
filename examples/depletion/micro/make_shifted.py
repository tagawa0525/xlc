"""燃焼点 k の組成を、点 k - 1 の組成のうち追う核種だけを点 k の値にしたものにした燃焼の履歴を書く（ADR-0216 の線形の補正の検査）。

    python3 make_shifted.py [--kind shifted|spread|previous] <履歴> <追う核種を micro に挙げた lattice の入力> <出力>

履歴は lattice の燃焼の履歴（`xlc-depletion-history`。ADR-0173）。追う核種は入力の `[history]` の `micro`（ADR-0216）で、どれも履歴が
輸送した核種であること（違えば止まる）。燃焼度 0 の点は、lattice が入力の組成と比べるのでそのまま。置き換えるのは各燃焼点の各領域の
数密度だけで、k∞ などのほかの項目はそのまま（`[history]` の分岐は数密度だけを読む。ADR-0187）。点 k の輸送の k∞ は、点 k - 1 の
少数群定数に追う核種の線形の補正 Σᵢ σᵢ(k - 1)(N̄ᵢ(k) − N̄ᵢ(k - 1)) を足した予測と比べる（compare_linearization.py）。

`--kind` は書く履歴（#0191 の候補を分ける輸送。compare_split.py）:

- shifted（既定）: 上の組成
- spread: shifted のうち追う核種の点 k の数密度を、燃焼の領域の体積で重みを付けた平均（N̄ᵢ(k)）を保ち、領域の間の比を点 k - 1 の
  ものにした組成（Nᵢ,r(k - 1) N̄ᵢ(k) / N̄ᵢ(k - 1)）。点 k - 1 で平均が 0 の核種は燃焼の領域で一様（微視的断面積の重み 1。ADR-0216 の
  案 Z1）。線形の補正の σᵢ は数密度で重みを付けた平均なので、この組成では数密度の差の空間の分布が補正の重みと同じ。負の数密度は
  lattice と同じく 0 として扱う
- previous: 点 k の数密度を点 k - 1 のもの（すべての核種）にした履歴。燃焼度は元のまま。lattice の `[history]` の `shielding-file`
  に渡し、自己遮蔽を前の点の組成で解く
"""

import sys

import msgpack
import tomllib


def shifted(history, tracked):
    if history.get("format") != "xlc-depletion-history":
        raise SystemExit("the history must be xlc-depletion-history")
    missing = [n for n in tracked if n not in history["transportable"]]
    if missing:
        raise SystemExit(f"not transported in the history: {missing}")
    indices = [history["nuclides"].index(n) for n in tracked]
    points = history["points"]
    # 後ろの点から書き換えると、点 k - 1 の組成はまだ元のまま
    for k in range(len(points) - 1, 0, -1):
        before = points[k - 1]["densities"]
        after = points[k]["densities"]
        for r, (old, new) in enumerate(zip(before, after)):
            row = list(old)
            for i in indices:
                row[i] = new[i]
            after[r] = row
    return history


def spread(history, tracked):
    history = shifted(history, tracked)
    indices = [history["nuclides"].index(n) for n in tracked]
    volumes = [d["volume"] for d in history["regions"]]
    total = sum(volumes)

    def mean(densities, i):
        return sum(v * max(n[i], 0.0) for v, n in zip(volumes, densities)) / total

    points = history["points"]
    # 点 k - 1 は shifted で書き換えたもの（追う核種は点 k - 1 の値）なので、後ろの点から書き換えて k - 1 の元の分布を使う
    for k in range(len(points) - 1, 0, -1):
        before = points[k - 1]["densities"]
        after = points[k]["densities"]
        for i in indices:
            target, base = mean(after, i), mean(before, i)
            for r, row in enumerate(after):
                row[i] = (
                    target * max(before[r][i], 0.0) / base if base > 0.0 else target
                )
    return history


def previous(history):
    if history.get("format") != "xlc-depletion-history":
        raise SystemExit("the history must be xlc-depletion-history")
    points = history["points"]
    for k in range(len(points) - 1, 0, -1):
        points[k]["densities"] = [list(row) for row in points[k - 1]["densities"]]
    return history


def tracked(text):
    micro = tomllib.loads(text).get("history", {}).get("micro")
    if not micro:
        raise SystemExit("the input lists no micro nuclides in [history]")
    return micro


def main():
    args = sys.argv[1:]
    kind = "shifted"
    if args[:1] == ["--kind"] and len(args) > 1:
        kind, args = args[1], args[2:]
    if len(args) != 3 or kind not in ("shifted", "spread", "previous"):
        raise SystemExit(__doc__)
    source, toml, output = args
    with open(source, "rb") as f:
        history = msgpack.unpack(f, raw=False, strict_map_key=False)
    with open(toml, encoding="utf-8") as f:
        nuclides = tracked(f.read())
    if kind == "previous":
        written = previous(history)
    elif kind == "spread":
        written = spread(history, nuclides)
    else:
        written = shifted(history, nuclides)
    with open(output, "wb") as f:
        msgpack.pack(written, f)


if __name__ == "__main__":
    main()
