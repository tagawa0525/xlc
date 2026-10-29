"""燃焼点 k の組成を、点 k - 1 の組成のうち追う核種だけを点 k の値にしたものにした燃焼の履歴を書く（ADR-0216 の線形の補正の検査）。

    python3 make_shifted.py <履歴> <追う核種を micro に挙げた lattice の入力> <出力>

履歴は lattice の燃焼の履歴（`xlc-depletion-history`。ADR-0173）。追う核種は入力の `[history]` の `micro`（ADR-0216）で、どれも履歴が
輸送した核種であること（違えば止まる）。燃焼度 0 の点は、lattice が入力の組成と比べるのでそのまま。置き換えるのは各燃焼点の各領域の
数密度だけで、k∞ などのほかの項目はそのまま（`[history]` の分岐は数密度だけを読む。ADR-0187）。点 k の輸送の k∞ は、点 k - 1 の
少数群定数に追う核種の線形の補正 Σᵢ σᵢ(k - 1)(N̄ᵢ(k) − N̄ᵢ(k - 1)) を足した予測と比べる（compare_linearization.py）。
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


def tracked(text):
    micro = tomllib.loads(text).get("history", {}).get("micro")
    if not micro:
        raise SystemExit("the input lists no micro nuclides in [history]")
    return micro


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    source, toml, output = sys.argv[1:]
    with open(source, "rb") as f:
        history = msgpack.unpack(f, raw=False, strict_map_key=False)
    with open(toml, encoding="utf-8") as f:
        nuclides = tracked(f.read())
    with open(output, "wb") as f:
        msgpack.pack(shifted(history, nuclides), f)


if __name__ == "__main__":
    main()
