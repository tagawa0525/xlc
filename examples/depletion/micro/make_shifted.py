"""燃焼点 k の組成を、点 k - 1 の組成のうち追う核種だけを点 k の値にしたものにした燃焼の履歴を書く（ADR-0216 の線形の補正の検査）。

python3 make_shifted.py <履歴> <追う核種を micro に挙げた lattice の入力> <出力>
"""


def shifted(history, tracked):
    raise NotImplementedError


def tracked(text):
    raise NotImplementedError


def main():
    raise NotImplementedError


if __name__ == "__main__":
    main()
