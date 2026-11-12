"""自己遮蔽を省いた 2A の分岐のテーブルを、すべての点で解いたテーブルと比べ、差を暗黙の効きの 1 次の見込みと並べる（#0180）。

- 燃焼の軸（`skip`）: 2A-every.toml と 2A-skip-sparse.toml のテーブル（公称の 1 状態、27 燃焼点）。燃焼点ごとに、反応度の差
  Δρ = 1/k_every − 1/k_skip（pcm）、少数群定数の項目（`ITEMS`）とピン出力の形状関数（`pins.form`）の相対差の最大
  （compare_branch_shielding.py の `relative`）。見込みは Σᵢ Sᵢ δᵢ、打ち消しを数えない大きさは Σᵢ |Sᵢ δᵢ|。Sᵢ は核種 i の暗黙の効き（自己遮蔽を解き直した k の直接の摂動と
  固定した摂動の差を、数密度の相対の変化あたりの反応度にしたもの。compare_few_group.py の `k_pcm` / 0.1）で、燃焼度の近い
  ピンセルの結果（SENSITIVITY）、δᵢ = (N̄ᵢ(S) − N̄ᵢ(B)) / N̄ᵢ(B) は自己遮蔽した点 S と燃焼点 B の、履歴のすべての燃焼の領域を
  体積で均質化した数密度の相対の差（自己遮蔽の組成が状態の組成からずれる向き）
- 微分の摂動の軸（`derivative`）: 2A-derivative.toml と 2A-derivative-fixed.toml のテーブルの Pu-240 の数密度の微分の差を、
  Pu-240 の数密度を 1 割変えたときの 2 群の k∞ の反応度の差にする（compare_few_group.py の `k_gradient`。微分 × 0.1 N̄ が定数の変化）

使い方: python3 compare_skip.py --json <要約の出力先>（このディレクトリで）
"""

import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import msgpack
import tomllib

HERE = Path(__file__).resolve().parent
HISTORY = HERE / "../../p9-product/lattice-2A.msgpack"
ITEMS = ["absorption", "nu_fission", "diffusion", "scatter", "adf"]
# 燃焼度 (MWd/kg) ごとの暗黙の効きの結果。燃焼点にはいちばん近い燃焼度のものを使い、2 つから同じ距離なら低い方（12 は 4、
# 40 は 20。`nearest`）
SENSITIVITY = {
    4.0: HERE / "pincell-few-group-4-sensitivity.toml",
    20.0: HERE / "pincell-few-group-20-sensitivity.toml",
    60.0: HERE / "../../issue-0279/pincell-few-group-sensitivity.toml",
}
FRACTION = 0.1


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def homogenized(history: dict, burnup: float) -> dict[str, float]:
    """燃焼点 `burnup` の、すべての燃焼の領域を体積で均質化した数密度（負の値は 0。分岐が解く物質と同じ）"""
    point = next((p for p in history["points"] if p["burnup"] == burnup), None)
    if point is None:
        sys.exit(f"the history has no point at {burnup} MWd/kg")
    volume = sum(r["volume"] for r in history["regions"])
    return {
        name: sum(
            max(d[i], 0.0) * r["volume"]
            for r, d in zip(history["regions"], point["densities"])
        )
        / volume
        for i, name in enumerate(history["nuclides"])
    }


def estimate(
    implicit: dict[str, float], shield: dict[str, float], state: dict[str, float]
) -> tuple[float, list[tuple[str, float]]]:
    """暗黙の効き `implicit`（核種ごと、数密度の相対の変化あたりの pcm）による、組成 `shield` で自己遮蔽した組成 `state` の
    反応度の差の 1 次の見込み (pcm) と、核種ごとの項（絶対値の大きい順）"""
    terms = []
    for name, s in implicit.items():
        if state[name] <= 0.0:
            if shield[name] != state[name]:
                sys.exit(
                    f"{name} has no density in the state but some in the shielding"
                )
            continue
        terms.append((name, s * (shield[name] - state[name]) / state[name]))
    terms.sort(key=lambda t: -abs(t[1]))
    return sum(t for _, t in terms), terms


def bound(terms: list[tuple[str, float]]) -> float:
    """項 `terms` の打ち消しを数えない大きさ Σ |項| (pcm)。項ごとの効きの誤差が打ち消しを崩しても超えない目安"""
    return sum(abs(t) for _, t in terms)


def shielding_burnup(path: Path) -> list[float]:
    """テーブルを作った入力 `path` の `[history] shielding-burnup`（自己遮蔽を解いた燃焼点）"""
    points = tomllib.loads(path.read_text()).get("history", {}).get("shielding-burnup")
    if points is None:
        sys.exit(f"{path} has no [history] shielding-burnup")
    return points


def implicit_effects(few_group, path: Path) -> dict[str, float]:
    """`[sensitivity.few-group]` の結果 `path` の、自己遮蔽を解き直した核種の暗黙の効き（相対の変化あたりの pcm）"""
    report = tomllib.loads(path.read_text())
    return {
        n["name"]: n["k_pcm"] / FRACTION for n in few_group.compare(report)["implicit"]
    }


def nearest(burnup: float) -> float:
    return min(SENSITIVITY, key=lambda b: abs(b - burnup))


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def table(path: Path) -> dict:
    return msgpack.unpackb(path.read_bytes(), strict_map_key=False)


def burnups(t: dict) -> list[float]:
    axis = t["axes"][0]
    if axis["name"] != "burnup":
        sys.exit("the first axis of the table is not the burnup")
    if len(t["states"]) != len(axis["points"]):
        sys.exit("the comparison takes one state per burnup point")
    return axis["points"]


def skip(shielding: dict, inputs: list) -> list[dict]:
    """燃焼の軸の比較（燃焼点ごと）"""
    branch = load("branch_shielding", HERE / "../../compare_branch_shielding.py")
    few_group = load("few_group", HERE / "../../issue-0279/compare_few_group.py")
    every_path, skip_path = (
        HERE / "table-2A-every.msgpack",
        HERE / "table-2A-skip-sparse.msgpack",
    )
    every, skipped = table(every_path), table(skip_path)
    points = burnups(every)
    if burnups(skipped) != points:
        sys.exit("the tables have different burnup points")
    effects = {b: implicit_effects(few_group, p) for b, p in SENSITIVITY.items()}
    with HISTORY.open("rb") as f:
        history = msgpack.unpackb(f.read(), strict_map_key=False)
    inputs += [
        {"path": os.path.relpath(p, HERE), "md5": md5(p)}
        for p in [every_path, skip_path, HISTORY, *SENSITIVITY.values()]
    ]
    rows = []
    for b, e, s in zip(points, every["states"], skipped["states"]):
        source = max(p for p in shielding if p <= b)
        state = homogenized(history, b)
        total, terms = estimate(
            effects[nearest(b)], homogenized(history, source), state
        )
        rows.append(
            {
                "burnup": b,
                "shielded_at": source,
                "drho_pcm": (1.0 / e["k_inf"] - 1.0 / s["k_inf"]) * 1e5,
                "relative": {item: branch.relative(e[item], s[item]) for item in ITEMS}
                | {"pins": branch.relative(e["pins"]["form"], s["pins"]["form"])},
                "estimate_pcm": total,
                "bound_pcm": bound(terms),
                "sensitivity_burnup": nearest(b),
                "terms": terms[:3],
            }
        )
    return rows


def two_group(
    state: dict, change: dict | None = None, scale: float = 0.0
) -> dict[str, float]:
    """状態 `state` の 2 群の定数を compare_few_group.py の応答の名前で（`change` があれば `scale` 倍を足す）"""

    def at(item: str, *index: int) -> float:
        value = state[item]
        delta = change[item] if change is not None else None
        for i in index:
            value = value[i]
            delta = delta[i] if delta is not None else None
        return value + (scale * delta if delta is not None else 0.0)

    values = {}
    for g in range(2):
        values[f"total-{g + 1}"] = at("total", g)
        values[f"nu-fission-{g + 1}"] = at("nu_fission", g)
        for h in range(2):
            values[f"scatter-{g + 1}-{h + 1}"] = at("scatter", g, h)
    return values


def derivative(inputs: list) -> list[dict]:
    """微分の摂動の軸の比較（燃焼点ごと）: Pu-240 の数密度を 1 割変えたときの 2 群の k∞ の反応度の変化を、摂動の自己遮蔽を
    解き直した微分と解き直さない微分で求め、その差（解き直さない − 解き直す）を書く"""
    few_group = load("few_group", HERE / "../../issue-0279/compare_few_group.py")
    paths = (
        HERE / "table-2A-derivative.msgpack",
        HERE / "table-2A-derivative-fixed.msgpack",
    )
    solved, fixed = (table(p) for p in paths)
    inputs += [{"path": os.path.relpath(p, HERE), "md5": md5(p)} for p in paths]
    points = burnups(solved)
    if burnups(fixed) != points:
        sys.exit("the derivative tables have different burnup points")
    # 固定した微分は、すべての燃焼点で自己遮蔽を解いた（その点の組成に固定した）入力のもの
    if shielding_burnup(HERE / "2A-derivative-fixed.toml") != points:
        sys.exit("2A-derivative-fixed.toml does not shield at every burnup point")
    rows = []
    for k, b in enumerate(points):
        row = {"burnup": b}
        for name, t in [("solved", solved), ("fixed", fixed)]:
            micro = t["micro"]
            i = micro["nuclides"].index("Pu-240")
            density = micro["densities"][k][i]
            state = t["states"][k]
            chi = state["chi"]
            k0 = few_group.two_group_k(two_group(state), chi)
            changed = few_group.two_group_k(
                two_group(
                    state, micro["states"][k][i]["derivative"], FRACTION * density
                ),
                chi,
            )
            row[name] = {
                "k_inf": state["k_inf"],
                "drho_pcm": (1.0 / k0 - 1.0 / changed) * 1e5,
            }
        row["difference_pcm"] = row["fixed"]["drho_pcm"] - row["solved"]["drho_pcm"]
        rows.append(row)
    return rows


def main() -> None:
    args = sys.argv[1:]
    if len(args) != 2 or args[0] != "--json":
        raise SystemExit(__doc__)
    inputs = []
    sparse = shielding_burnup(HERE / "2A-skip-sparse.toml")
    summary = {
        "format": "xlc-compare-shielding-skip",
        "shielding_burnup": sparse,
        "skip": skip(sparse, inputs),
        "derivative": derivative(inputs),
        "inputs": inputs,
    }
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
    print(
        f"{'burnup':>6s} {'from':>5s} {'drho':>9s} {'estimate':>9s} {'ratio':>6s} {'bound':>7s} {'pins':>8s} {'nu_fiss':>8s} {'absorb':>8s}  main terms"
    )
    for r in summary["skip"]:
        ratio = r["drho_pcm"] / r["estimate_pcm"] if r["estimate_pcm"] else float("nan")
        terms = ", ".join(f"{n} {t:+.2f}" for n, t in r["terms"])
        print(
            f"{r['burnup']:6g} {r['shielded_at']:5g} {r['drho_pcm']:+9.3f} {r['estimate_pcm']:+9.3f} {ratio:6.2f} {r['bound_pcm']:7.2f} "
            f"{r['relative']['pins']:8.1e} {r['relative']['nu_fission']:8.1e} {r['relative']['absorption']:8.1e}  {terms}"
        )
    print(
        "Pu-240 density +10 % reactivity from the derivatives: solved, fixed, fixed - solved (pcm)"
    )
    for r in summary["derivative"]:
        print(
            f"  {r['burnup']:g} MWd/kg: {r['solved']['drho_pcm']:+.3f} {r['fixed']['drho_pcm']:+.3f} "
            f"{r['difference_pcm']:+.3f} (k {r['solved']['k_inf']:.6f} {r['fixed']['k_inf']:.6f})"
        )
    for i in inputs:
        print(f"input {i['path']} md5 {i['md5']}")


if __name__ == "__main__":
    main()
