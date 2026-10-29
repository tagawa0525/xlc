"""core の燃焼の履歴（`xlc-core-history`。ADR-0190、ADR-0218）を読む比較のツールの共通の部分。"""


def single_state(core: dict) -> dict:
    """状態が 1 つで定格の出力のままの燃焼の、その状態（履歴の `states` の値）。状態が 1 つでないか、出力の割合が 1 でないステップ
    （崩壊だけのステップは 0）があれば止める（ステップごとに状態や出力の変わる燃焼は、すべての点を定格の燃焼の履歴と 1 つの状態の
    テーブルの値と比べるツールで比べない。ADR-0218）"""
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
