"""燃焼チェーンの入力と多群ライブラリ（xs）の入力が、天然の元素の評価（A = 0）を当てる核種で揃っていることを確かめる
（ADR-0262 の帰結。チェーンも xs も相手の入力を読まないので、どちらのツールも食い違いに気づかない）。

使い方: python3 check_natural_isotopes.py <チェーンの入力> <チェーン（msgpack）> <xs の入力>...

xs の入力ごとに、各核種の評価の MF1/MT451 の ZA（HEAD）と LISO（2 つ目の記録）を読み、次を確かめる。食い違いがあればすべてを標準エラーに出して 1 で終わる。

- チェーンの `neutron-directory` と xs の `endf-directory` が同じディレクトリ（同じ版の中性子のファイル）を指す
- 核種がチェーンにある（チェーンの核種の名前、Z、A、状態で照らす）
- 評価が天然の元素（A = 0）なら、同じ Z の核種で、チェーンの `natural-isotopes` が挙げる。核種ごとの評価なら、ZA と LISO がその核種の
  Z、A、状態のもので、`natural-isotopes` が挙げない

標準出力には、xs の入力ごとに天然の元素の評価を当てた核種と、`natural-isotopes` のうちライブラリに無い核種を書く。
"""

import itertools
import sys
from pathlib import Path

import msgpack
import tomllib


def endf_number(field: str) -> float:
    """ENDF の 11 桁の数（`9.223500+4` のように指数の e を省く形も読む）"""
    s = field.strip()
    for i in range(1, len(s)):
        if s[i] in "+-" and s[i - 1] not in "eE":
            return float(s[:i] + "e" + s[i:])
    return float(s)


def endf_identity(path: Path) -> tuple[int, int, int]:
    """評価の MF1/MT451 の最初の記録（HEAD）の ZA と 2 つ目の記録の LISO を (Z, A, LISO) にする（ファイルの先頭の 2 つの記録だけを読む）"""
    with path.open() as f:
        records = list(
            itertools.islice(
                (
                    line
                    for line in f
                    if line[70:72].strip() == "1" and line[72:75].strip() == "451"
                ),
                2,
            )
        )
    if len(records) < 2:
        raise SystemExit(f"{path}: no MF1/MT451 HEAD and CONT records")
    za = round(endf_number(records[0][0:11]))
    try:
        liso = int(records[1][33:44])
    except ValueError:
        raise SystemExit(f"{path}: no LISO in the second MF1/MT451 record") from None
    return za // 1000, za % 1000, liso


def check(
    chain_input: dict,
    chain_dir: Path,
    chain: dict,
    xs_input: dict,
    xs_dir: Path,
) -> dict:
    """1 つの xs の入力を、チェーンの入力とチェーンに照らす。`errors` が空なら揃っている"""
    errors = []
    neutrons = (chain_dir / chain_input["neutron-directory"]).resolve()
    endf = (xs_dir / xs_input["endf-directory"]).resolve()
    if neutrons != endf:
        errors.append(f"the chain reads {neutrons}, the library {endf}")
    listed = list(chain_input.get("natural-isotopes", []))
    nuclides = {n["name"]: (n["z"], n["a"], n["state"]) for n in chain["nuclides"]}
    natural = []
    for name, entry in xs_input["nuclide"].items():
        if name not in nuclides:
            errors.append(f"{name} is not in the chain")
            continue
        z, a, liso = endf_identity(endf / entry["endf"])
        if a == 0:
            if z != nuclides[name][0]:
                errors.append(
                    f"{name}: {entry['endf']} is the natural evaluation of Z = {z}"
                )
            elif name not in listed:
                errors.append(
                    f"{name}: the library takes the natural evaluation {entry['endf']}, "
                    "but natural-isotopes of the chain does not list it"
                )
            else:
                natural.append([name, entry["endf"]])
        elif (z, a, liso) != nuclides[name]:
            errors.append(
                f"{name}: {entry['endf']} is the evaluation of ZA {z}{a:03d} LISO {liso}"
            )
        elif name in listed:
            errors.append(
                f"{name}: natural-isotopes of the chain lists it, but the library takes "
                f"its own evaluation {entry['endf']}"
            )
    in_library = set(xs_input["nuclide"])
    return {
        "natural": natural,
        "chain_only": [n for n in listed if n not in in_library],
        "errors": errors,
    }


def main() -> None:
    args = sys.argv[1:]
    if len(args) < 3:
        raise SystemExit(__doc__)
    chain_path, chain_msgpack, *xs_paths = map(Path, args)
    chain_input = tomllib.loads(chain_path.read_text())
    chain = msgpack.unpackb(chain_msgpack.read_bytes())
    failed = False
    for xs_path in xs_paths:
        result = check(
            chain_input,
            chain_path.parent,
            chain,
            tomllib.loads(xs_path.read_text()),
            xs_path.parent,
        )
        natural = ", ".join(f"{n} ← {f}" for n, f in result["natural"]) or "なし"
        print(f"{xs_path}: 天然の元素の評価を当てた核種: {natural}")
        print(
            f"{xs_path}: natural-isotopes のうちライブラリに無い核種: "
            + (", ".join(result["chain_only"]) or "なし")
        )
        for e in result["errors"]:
            print(f"{xs_path}: {e}", file=sys.stderr)
            failed = True
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
