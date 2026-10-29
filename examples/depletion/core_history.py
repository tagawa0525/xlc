"""core の燃焼の履歴（`xlc-core-history`。ADR-0190、ADR-0218）を読む比較のツールの共通の部分。"""


def single_state(core: dict) -> dict:
    """状態が 1 つの燃焼の、その状態（履歴の `states` の値）。状態が 1 つでなければ止める（ステップごとに状態の変わる燃焼は、
    すべての点を 1 つの状態のテーブルの値と比べるツールで比べない）"""
    states = core["states"]
    if len(states) != 1:
        raise SystemExit(
            f"the core history has the states {sorted(states)}; this tool compares a depletion of one state"
        )
    (state,) = states.values()
    return state
