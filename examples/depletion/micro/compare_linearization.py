"""ミクロ燃焼の線形の補正（ADR-0216）が、追う核種の組成の差の輸送の k∞ をどれだけ再現するかを、燃焼度の隣り合う点の組ごとに出す。"""


def k_inf(state):
    raise NotImplementedError


def corrected(state, sigmas, deltas):
    raise NotImplementedError


def compare(nominal, shifted):
    raise NotImplementedError
