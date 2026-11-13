"""燃料の種類ごとに、微分の摂動の自己遮蔽を解き直したテーブルと固定したテーブル（inputs.py）の Pu-240 の数密度の微分を比べる
（#0180、ADR-0249）。

燃焼点ごとに、Pu-240 の数密度を 1 割変えたときの 2 群の k∞ の反応度（微分 × 0.1 N̄ が定数の変化。../reuse/compare_skip.py の
`derivative` と同じ量）を両方のテーブルで求め、差（固定 − 解き直し）と、解き直した値に対する割合を書く。

使い方: python3 compare.py --json <要約の出力先>（このディレクトリで）
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FRACTION = 0.1


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SKIP = load("compare_skip", HERE / "../reuse/compare_skip.py")
FEW_GROUP = load("few_group", HERE / "../../issue-0279/compare_few_group.py")
CASES = load("inputs", HERE / "inputs.py").CASES


def reactivity(t: dict, k: int, nuclide: str = "Pu-240") -> dict:
    """テーブル `t` の燃焼点 `k` の、核種 `nuclide` の数密度を FRACTION 倍だけ変えたときの 2 群の k∞ の反応度の変化（pcm）"""
    micro = t["micro"]
    i = micro["nuclides"].index(nuclide)
    density = micro["densities"][k][i]
    state = t["states"][k]
    chi = state["chi"]
    k0 = FEW_GROUP.two_group_k(SKIP.two_group(state), chi)
    changed = FEW_GROUP.two_group_k(
        SKIP.two_group(state, micro["states"][k][i]["derivative"], FRACTION * density),
        chi,
    )
    return {"k_inf": state["k_inf"], "drho_pcm": (1.0 / k0 - 1.0 / changed) * 1e5}


def share(difference: float, solved: float) -> float | None:
    """差の、解き直した値に対する割合（解き直した値が 0 なら割れないので None）"""
    return difference / solved if solved != 0.0 else None


def same_nominal(solved: dict, fixed: dict) -> bool:
    """2 つのテーブルの状態（k∞ と少数群定数のすべての項目）と、追う核種と均質化した数密度が同じか（摂動の自己遮蔽だけが違い、
    状態は同じに解くはず）"""
    return solved["states"] == fixed["states"] and all(
        solved["micro"][key] == fixed["micro"][key] for key in ("nuclides", "densities")
    )


def compare(name: str, inputs: list) -> list[dict]:
    paths = [HERE / f"table-{name}-{v}.msgpack" for v in ("solved", "fixed")]
    solved, fixed = (SKIP.table(p) for p in paths)
    inputs += [{"path": p.name, "md5": SKIP.md5(p)} for p in paths]
    points = SKIP.burnups(solved)
    if SKIP.burnups(fixed) != points:
        sys.exit(f"{name}: the tables have different burnup points")
    if not same_nominal(solved, fixed):
        sys.exit(f"{name}: the tables solve different states")
    rows = []
    for k, b in enumerate(points):
        row = {
            "burnup": b,
            "solved": reactivity(solved, k),
            "fixed": reactivity(fixed, k),
        }
        row["difference_pcm"] = row["fixed"]["drho_pcm"] - row["solved"]["drho_pcm"]
        row["fraction"] = share(row["difference_pcm"], row["solved"]["drho_pcm"])
        rows.append(row)
    return rows


def main() -> None:
    args = sys.argv[1:]
    if len(args) != 2 or args[0] != "--json":
        raise SystemExit(__doc__)
    inputs = []
    summary = {
        "format": "xlc-compare-derivative-shielding-fuels",
        "cases": {name: compare(name, inputs) for name in CASES},
        "inputs": inputs,
    }
    Path(args[1]).write_text(json.dumps(summary, indent=1) + "\n")
    print(
        "Pu-240 density +10 % reactivity: solved, fixed, fixed - solved (pcm), fraction"
    )
    for name, rows in summary["cases"].items():
        for r in rows:
            print(
                f"  {name} {r['burnup']:g} MWd/kg: {r['solved']['drho_pcm']:+.3f} {r['fixed']['drho_pcm']:+.3f} "
                f"{r['difference_pcm']:+.3f} {'-' if r['fraction'] is None else format(r['fraction'], '+.2%')} (k {r['solved']['k_inf']:.6f} {r['fixed']['k_inf']:.6f})"
            )
    for i in inputs:
        print(f"input {i['path']} md5 {i['md5']}")


if __name__ == "__main__":
    main()
