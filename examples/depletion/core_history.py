"""core の燃焼の履歴（`xlc-core-history`。ADR-0190、ADR-0218）を読む比較のツールの共通の部分。"""

# 2 つの燃焼度（MWd/kg）を同じ燃焼点とみなす相対の許容（1 MWd/kg 未満は絶対の許容）。core の燃焼度は出力密度 × EFPD で、
# 格子計算の履歴は燃焼度の点から時間を求めて燃焼させ、炉心の燃焼度はノードの燃焼度の重金属の質量の重みの平均なので、
# 同じ点でも丸めが違う
BURNUP_TOLERANCE = 1e-9


def same_burnup(a: float, b: float) -> bool:
    """燃焼度 `a` が `b` と同じ点か（`BURNUP_TOLERANCE`）"""
    return abs(a - b) <= BURNUP_TOLERANCE * max(1.0, abs(b))


def rho_pcm(base: float, k: float) -> float:
    """`base` から `k` への反応度の差 (pcm)"""
    return (1.0 / base - 1.0 / k) * 1e5


def default_midpoints(core: dict) -> None:
    """中点を初めの点のホウ素で解いた燃焼（ADR-0191 の案 S1）と、Xe を平衡にした燃焼（ADR-0219）なら止める（中点でも探し、Xe を時間で
    追う燃焼を比べるツールで比べない）"""
    if core.get("midpoint_boron") == "start":
        raise SystemExit(
            "the core history solves the midpoints at the start boron (ADR-0191 S1); this tool compares the default search"
        )
    if core.get("xenon") == "equilibrium":
        raise SystemExit(
            "the core history balances its xenon at each solve (ADR-0219); this tool compares xenon tracked in time"
        )


def single_state(core: dict) -> dict:
    """状態が 1 つで定格の出力のままの燃焼の、その状態（履歴の `states` の値）。状態が 1 つでないか、出力の割合が 1 でないステップ
    （崩壊だけのステップは 0）があれば止める（ステップごとに状態や出力の変わる燃焼は、すべての点を定格の燃焼の履歴と 1 つの状態の
    テーブルの値と比べるツールで比べない。ADR-0218）"""
    default_midpoints(core)
    states = core["states"]
    if len(states) != 1:
        raise SystemExit(
            f"the core history has the states {sorted(states)}; this tool compares a depletion of one state"
        )
    powers = sorted({p["power"] for p in core["points"] if p.get("power", 1.0) != 1.0})
    if powers:
        raise SystemExit(
            f"the core history has steps at the power fractions {powers}; this tool compares a depletion at the rated power"
        )
    (state,) = states.values()
    return state
