"""DVC の build のステージが bin/ に置く実行ファイルから、ソースの行と列の番号を消す（ADR-0255）。

使い方: python3 scripts/normalize_binary.py <ビルドした実行ファイル> <出力>

計算を変えない編集（注釈の行を足す、関数を並べ替える）でも、ソースの行の番号は実行ファイルの次の 3 か所に入り、
ステージの deps（ADR-0159）を変える。どれも計算のコードからは読まれず、パニックの文とログの記述にだけ使われる:

- パニックの位置 `core::panic::Location`（配列の境界の検査など）の行と列
- `tracing` のイベントとスパンの `Metadata` の `line: Option<u32>`
- `tracing::enabled!` の名前 `"enabled <ファイル>:<行>"` の行の数字（イベントの名前は xlc_log のマクロが行を含めない）

この 3 か所の番号を 0 にしたコピーを書く。ほかのバイトは変えない。構造体のレイアウトは `repr(Rust)` で保証されないので、
値どうしの整合で見つける: Location は「NUL で終わる `.rs` の文字列を指す &str の後に 0 でない行と列」、Metadata は
「名前 `"event <ファイル>"` の &str の直前に Some(行)、後ろの決まった距離に同じファイルの &str」。距離はイベントから
決め、スパンにも使う。整合が崩れたとき（Location か Metadata が見つからない、距離がイベントごとに違う、1 つの &str が
Location と Metadata のファイル名の両方に当たる）は、
コンパイラか `tracing` のレイアウトが変わったとみなしてエラーで止める。黙って消し残すと、注釈の行でステージが
古くなるのに気付けず、黙って別のバイトを消すと計算が変わりうる。
crates/ のイベントの名前が行を含む（`tracing` のイベントのマクロを直接呼んだ）ときもエラーで止める。
"""

import re
import struct
import sys
from pathlib import Path

# tracing の Metadata の名前。xlc_log のイベントは "event <ファイル>"、tracing::enabled! は "enabled <ファイル>:<行>"。
# tracing のイベントのマクロを名前なしで直接呼ぶと "event <ファイル>:<行>" になる
NAME = re.compile(rb"(event|enabled) ([^\x00]+?\.rs)(?::([0-9]+))?")
# &str とみなす長さの上限。ソースのパスと名前はこれより短く、ほかの整数を長さと取り違える範囲を狭める
MAX_TEXT = 4096
# Metadata の名前の slot からファイルの slot までの距離の上限（今のレイアウトでは 80 バイト）
MAX_DISTANCE = 256
SHT_NOBITS = 8
SHF_ALLOC = 2


class Elf:
    """ELF64 リトルエンディアンの節の表だけを読む。"""

    def __init__(self, data: bytes):
        if data[:4] != b"\x7fELF" or data[4] != 2 or data[5] != 1:
            raise ValueError("ELF64 のリトルエンディアンでない")
        (shoff,) = struct.unpack_from("<Q", data, 0x28)
        shentsize, shnum, shstrndx = struct.unpack_from("<HHH", data, 0x3A)
        headers = [
            struct.unpack_from("<IIQQQQ", data, shoff + i * shentsize)
            for i in range(shnum)
        ]
        names = headers[shstrndx][4]
        self.sections = {}
        for name, kind, flags, addr, offset, size in headers:
            label = data[names + name : data.index(b"\x00", names + name)].decode()
            # NOBITS（.bss など）はファイルに中身が無い
            if kind != SHT_NOBITS:
                self.sections[label] = (flags, addr, offset, size)
        self.loaded = [
            (addr, offset, size)
            for flags, addr, offset, size in self.sections.values()
            if flags & SHF_ALLOC and size
        ]

    def file_offset(self, addr: int, size: int) -> int | None:
        """アドレス addr から size バイトが 1 つの載る節に収まるなら、そのファイルの位置。"""
        for base, offset, length in self.loaded:
            if base <= addr and addr + size <= base + length:
                return offset + addr - base
        return None


def normalize(data: bytes) -> tuple[bytes, dict[str, int]]:
    elf = Elf(data)
    if ".data.rel.ro" not in elf.sections:
        raise ValueError(".data.rel.ro が無い（位置の情報の置き場が変わった）")
    _, _, relro, relro_size = elf.sections[".data.rel.ro"]
    out = bytearray(data)

    def text(slot: int) -> tuple[bytes, int] | None:
        """slot の 16 バイトを &str（ptr、長さ）とみたときの中身と、中身のファイルの位置。"""
        ptr, length = struct.unpack_from("<QQ", data, slot)
        if not 0 < length <= MAX_TEXT:
            return None
        at = elf.file_offset(ptr, length)
        if at is None:
            return None
        return data[at : at + length], at

    def location(slot: int) -> bool:
        """slot が Location か: ファイル名が NUL で終わり（file_as_c_str のため）、行と列が 1 から数えた値。"""
        if slot + 24 > end:
            return False
        ptr, length = struct.unpack_from("<QQ", data, slot)
        at = elf.file_offset(ptr, length + 1)
        line, column = struct.unpack_from("<II", data, slot + 16)
        return at is not None and data[at + length] == 0 and line > 0 and column > 0

    counts = {"location": 0, "metadata": 0, "enabled": 0}
    end = relro + relro_size
    first = relro + (-relro) % 8
    # .rs で終わる &str（Location か Metadata のファイル名）と、Metadata の名前
    files = {}
    names = []
    for slot in range(first, end - 15, 8):
        found = text(slot)
        if found is None:
            continue
        body, at = found
        # 名前 "event <ファイル>" も .rs で終わるので、ファイル名より先に見る
        match = NAME.fullmatch(body)
        if match is not None:
            if slot - 8 >= first:
                names.append((slot, match, at))
        elif body.endswith(b".rs"):
            files[slot] = body

    distances = set()
    for slot, match, at in names:
        kind, path, digits = match.groups()
        tag, line = struct.unpack_from("<II", data, slot - 8)
        if tag != 1:
            raise ValueError(
                f"{match.group(0)!r} の直前が Some(行) でない（Metadata のレイアウトが変わった）"
            )
        if digits is not None:
            if kind == b"event" and path.startswith(b"crates/"):
                raise ValueError(
                    f"{match.group(0)!r}: tracing のイベントのマクロを名前なしで呼んでいる。xlc_log のマクロを使う"
                )
            if line != int(digits):
                raise ValueError(
                    f"{match.group(0)!r} の行と Metadata の line {line} が違う"
                )
            start, stop = match.span(3)
            out[at + start : at + stop] = b"0" * (stop - start)
            if kind == b"enabled":
                counts["enabled"] += 1
        # 同じファイルの次の Metadata のファイルも後ろにあるので、最も近いものをとる
        near = [
            other - slot
            for other in range(slot + 8, slot + MAX_DISTANCE + 8, 8)
            if files.get(other) == path
        ]
        if not near:
            raise ValueError(
                f"{match.group(0)!r} の Metadata にファイルの &str が見つからない"
            )
        distances.add(near[0])
    if len(distances) != 1:
        raise ValueError(
            f"Metadata の名前からファイルまでの距離が 1 つに決まらない: {sorted(distances)}（イベントが無いか、レイアウトが変わった）"
        )
    (distance,) = distances

    # Metadata のファイル名を先に決める。file!() の文字列の後ろが偶然 NUL だと、ファイル名の後ろのフィールドが
    # Location の行と列に見えるので、Location はそれ以外から探す
    for slot in files:
        name = slot - distance
        if name - 8 < first or text(name) is None:
            continue
        tag, line = struct.unpack_from("<II", data, name - 8)
        if tag != 1 or line == 0:
            continue
        if location(slot):
            raise ValueError(
                f"0x{slot:x} の &str が Location にも Metadata のファイル名にも当たる（レイアウトが変わったか、偶然の一致）"
            )
        struct.pack_into("<I", out, name - 4, 0)
        counts["metadata"] += 1
        files[slot] = None
    for slot, body in files.items():
        if body is not None and location(slot):
            struct.pack_into("<II", out, slot + 16, 0, 0)
            counts["location"] += 1

    if counts["location"] == 0:
        raise ValueError("Location が見つからない（Location のレイアウトが変わった）")
    return bytes(out), counts


def main(argv: list[str]) -> None:
    if len(argv) != 2:
        raise SystemExit("usage: python3 scripts/normalize_binary.py <built> <output>")
    src, dst = Path(argv[0]), Path(argv[1])
    data, counts = normalize(src.read_bytes())
    dst.write_bytes(data)
    dst.chmod(src.stat().st_mode)
    print(f"{dst}: " + " ".join(f"{k}={v}" for k, v in counts.items()))


if __name__ == "__main__":
    main(sys.argv[1:])
