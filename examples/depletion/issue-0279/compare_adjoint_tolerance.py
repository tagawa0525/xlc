"""同じ状態の少数群定数の数密度の微分（lattice の `[sensitivity.few-group]`。ADR-0247 の段 4）を、一般化随伴の判定の違う
2 つの結果の間で比べ、判定が微分に入れる差の大きさを書く（#0279）。

- 定数ごと: 一般化随伴の値 G（`relative`）の差を定数の値 R で割った |ΔG| / R の最大と、その核種と定数
- 核種ごと: 数密度を 1 割変えたときの反応度の差 Σ_R ∂k∞/∂R ΔG / k² (pcm) と、項の絶対値の和（k∞ と FRACTION は
  compare_few_group.py と同じ）

使い方: python3 compare_adjoint_tolerance.py <結果の TOML> <結果の TOML>
"""

import sys
from pathlib import Path

import tomllib
from compare_few_group import FRACTION, k_gradient


def differences(coarse: dict, fine: dict) -> dict:
    """`coarse` と `fine` の G の差の、定数に対する最大と、核種ごとの反応度の差（反応度の差の絶対値の大きい順）"""
    a, b = coarse["few_group"], fine["few_group"]
    # 判定のほかが同じ状態の結果でなければ、差が判定の分にならない
    if (
        coarse["k"] != fine["k"]
        or a["responses"] != b["responses"]
        or a["values"] != b["values"]
    ):
        sys.exit("the two results are not of the same state")
    names = [n["name"] for n in a["nuclides"]]
    if names != [n["name"] for n in b["nuclides"]]:
        sys.exit("the two results have different nuclides")
    responses = a["responses"]
    gradient = k_gradient(dict(zip(responses, a["values"])), a["chi"])
    scale = 1e5 * FRACTION / coarse["k"] ** 2
    largest = None
    nuclides = []
    for na, nb in zip(a["nuclides"], b["nuclides"]):
        terms = []
        for response, value, ga, gb in zip(
            responses, a["values"], na["relative"], nb["relative"]
        ):
            if value == 0.0:
                if ga != 0.0 or gb != 0.0:
                    sys.exit(
                        f"{na['name']} {response}: a derivative of a zero constant"
                    )
                continue
            relative = abs(ga - gb) / abs(value)
            if largest is None or relative > largest["relative"]:
                largest = {
                    "name": na["name"],
                    "response": response,
                    "relative": relative,
                }
            terms.append(scale * gradient[response] * (ga - gb))
        nuclides.append(
            {
                "name": na["name"],
                "pcm": sum(terms),
                "absolute": sum(abs(t) for t in terms),
            }
        )
    nuclides.sort(key=lambda n: -abs(n["pcm"]))
    return {"largest": largest, "nuclides": nuclides}


def main(arguments: list[str]) -> int:
    if len(arguments) != 2:
        sys.exit(__doc__)
    coarse, fine = (tomllib.loads(Path(p).read_text()) for p in arguments)
    for path, report in zip(arguments, (coarse, fine)):
        few = report["few_group"]
        print(
            f"{path}: GMRES tolerance {few['tolerance']:.0e}, sweeps {sum(few['sweeps'])}, "
            f"generalized {few['seconds']['generalized']:.1f} s"
        )
    result = differences(coarse, fine)
    largest = result["largest"]
    print(
        f"largest |dG| / R: {largest['relative']:.2e} ({largest['name']} {largest['response']})"
    )
    top = result["nuclides"][0]
    widest = max(result["nuclides"], key=lambda n: n["absolute"])
    print(
        f"reactivity of a 10 % density change: largest {top['pcm']:+.2e} pcm ({top['name']}), "
        f"largest sum of absolute terms {widest['absolute']:.2e} pcm ({widest['name']})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
