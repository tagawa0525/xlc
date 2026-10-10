"""格子の入力の設定を変えた写しを作る（トラック間隔やリングの数の依存を測るため）。

使い方: python3 lattice_variant.py <格子の入力> <出力ディレクトリ> <名前> [変更 ...]

<出力ディレクトリ>/<名前>.toml に写しを書く。ライブラリは元の入力からの絶対パスに、ピン出力は
<出力ディレクトリ>/pins-<名前>.msgpack にする。変更は次の形で、どれも元の入力にちょうど 1 つ（rings:<物質> は 1 つ以上）
当てはまらなければ止まる:

- `azimuthal=128`、`spacing=0.01` など: 最上位の `キー = 値` の行を置き換える
- `rings:<物質>=<数>`: 半径のある `<物質>` の帯のリングの数を変える（その物質の帯すべて）
- `rings@<半径>=<数>`: 外側の半径が <半径> の帯のリングの数を変える

解くのは `lattice <出力ディレクトリ>/<名前>.toml`（k は標準出力）。
"""

import re
import sys
from pathlib import Path


def change(text: str, edit: str) -> str:
    key, value = edit.split("=")
    if key.startswith("rings:"):
        material = re.escape(key[len("rings:") :])
        pattern = rf'(\{{ material = "{material}", radius = [0-9.]+, rings = )\d+'
        expected = None
    elif key.startswith("rings@"):
        radius = re.escape(key[len("rings@") :])
        pattern = rf'(\{{ material = "[^"]+", radius = {radius}, rings = )\d+'
        expected = 1
    else:
        text, n = re.subn(
            rf"^{re.escape(key)} = .*$", f"{key} = {value}", text, flags=re.MULTILINE
        )
        if n != 1:
            raise SystemExit(f"{edit}: expected one top-level `{key}`, found {n}")
        return text
    text, n = re.subn(pattern, rf"\g<1>{value}", text)
    if n == 0 or (expected is not None and n != expected):
        raise SystemExit(f"{edit}: matched {n} zones")
    return text


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    source = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2]).resolve()
    name = sys.argv[3]
    text = source.read_text()

    library = re.search(r'^library = "([^"]*)"$', text, flags=re.MULTILINE)
    if library is None:
        raise SystemExit(f"{source}: no `library`")
    absolute = (source.parent / library.group(1)).resolve()
    text = text.replace(library.group(0), f'library = "{absolute}"')
    text, n = re.subn(
        r'^pin-powers = "[^"]*"$',
        f'pin-powers = "{out / f"pins-{name}.msgpack"}"',
        text,
        flags=re.MULTILINE,
    )
    if n > 1:
        raise SystemExit(f"{source}: more than one `pin-powers`")
    for edit in sys.argv[4:]:
        text = change(text, edit)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{name}.toml").write_text(text)


if __name__ == "__main__":
    main()
