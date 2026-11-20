"""版ごとの `xs check` の報告（ADR-0266）を、燃焼チェーンの核種とその外に分けて数える。"""


def nuclide_key(path: str) -> tuple[int, int, int]:
    return (0, 0, 0)


def summarize(versions: list) -> dict:
    return {}
