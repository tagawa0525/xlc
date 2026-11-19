"""normalize_binary.py の検査（bin/ の実行ファイルから行と列の番号を除くコピー）。

実行: python3 -m unittest discover -s scripts -p 'test_normalize_binary.py'

実行ファイルの代わりに、節の表と .rodata、.data.rel.ro だけの小さな ELF を組み立てる。構造体の並びは rustc 1.99.0 と
tracing 0.1.44 のもの（ADR-0255 の「背景」）: Location は（ファイル名の ptr、長さ、行、列）、Metadata は先頭に line の
Option<u32>、8 バイト目に名前の &str、名前から 0x50 バイト後ろにファイルの &str。
"""

import os
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import normalize_binary as norm

RODATA_ADDR, RODATA_OFFSET = 0x1000, 0x1000
# .data.rel.ro はアドレスとファイルの位置を違えて、アドレスからファイルの位置への読み替えも確かめる
RELRO_ADDR, RELRO_OFFSET = 0x5000, 0x2000
SIZE = 0x1000
METADATA = 0x70
FILE_DISTANCE = 0x50


class Image:
    """.rodata に文字列、.data.rel.ro に構造体を置いた ELF64 を組み立てる。"""

    def __init__(self):
        self.rodata = bytearray()
        self.relro = bytearray()

    def string(self, text: bytes, nul: bool = False) -> tuple[int, int]:
        """文字列を置き、(アドレス、長さ)。nul でなければ後ろに 'x' を置き、偶然の NUL を避ける。"""
        addr = RODATA_ADDR + len(self.rodata)
        self.rodata += text + (b"\x00" if nul else b"x")
        return addr, len(text)

    def last_string(self, text: bytes) -> tuple[int, int]:
        """.rodata の最後のバイトで終わるように文字列を置く。"""
        self.rodata += b"y" * (SIZE - len(text) - len(self.rodata))
        addr = RODATA_ADDR + len(self.rodata)
        self.rodata += text
        return addr, len(text)

    def place(self, data: bytes) -> int:
        """8 バイト境界に構造体を置き、.data.rel.ro の中の位置。"""
        self.relro += b"\x00" * (-len(self.relro) % 8)
        at = len(self.relro)
        self.relro += data
        return at

    def location(self, path: bytes, line: int, column: int) -> int:
        ptr, length = self.string(path, nul=True)
        return self.place(struct.pack("<QQII", ptr, length, line, column))

    def metadata(self, name: bytes, path: bytes, line: int) -> int:
        """tracing の Metadata。名前の後ろのフィールド（target、fields、callsite、module_path）は 0 でない値で埋める。"""
        name_ptr, name_len = self.string(name)
        file_ptr, file_len = self.string(path)
        module_ptr, module_len = self.string(b"xlc_test::module")
        data = bytearray(METADATA)
        struct.pack_into("<II", data, 0, 1, line)
        struct.pack_into("<QQ", data, 8, name_ptr, name_len)
        struct.pack_into("<QQ", data, 0x18, module_ptr, module_len)
        struct.pack_into("<QQQQ", data, 0x28, 0x7777, 2, 0x8888, 0x9999)
        struct.pack_into("<QQ", data, 0x48, module_ptr, module_len)
        struct.pack_into("<QQ", data, 8 + FILE_DISTANCE, file_ptr, file_len)
        struct.pack_into("<Q", data, 0x68, 0x0102)
        return self.place(bytes(data))

    def build(self) -> bytes:
        names = b"\x00.rodata\x00.data.rel.ro\x00.shstrtab\x00"
        shstrtab_offset = RELRO_OFFSET + SIZE
        shoff = shstrtab_offset + len(names) + (-(shstrtab_offset + len(names)) % 8)
        out = bytearray(shoff + 4 * 64)
        out[:4] = b"\x7fELF"
        out[4], out[5], out[6] = 2, 1, 1
        struct.pack_into("<Q", out, 0x28, shoff)
        struct.pack_into("<HHH", out, 0x3A, 64, 4, 3)
        out[RODATA_OFFSET : RODATA_OFFSET + len(self.rodata)] = self.rodata
        out[RELRO_OFFSET : RELRO_OFFSET + len(self.relro)] = self.relro
        out[shstrtab_offset : shstrtab_offset + len(names)] = names
        sections = [
            (0, 0, 0, 0, 0, 0),
            (1, 1, 2, RODATA_ADDR, RODATA_OFFSET, SIZE),
            (9, 1, 3, RELRO_ADDR, RELRO_OFFSET, SIZE),
            (22, 3, 0, 0, shstrtab_offset, len(names)),
        ]
        for i, header in enumerate(sections):
            struct.pack_into("<IIQQQQ", out, shoff + i * 64, *header)
        return bytes(out)


def u32(data: bytes, relro_at: int) -> int:
    return struct.unpack_from("<I", data, RELRO_OFFSET + relro_at)[0]


def changed(before: bytes, after: bytes) -> set[int]:
    return {i for i, (a, b) in enumerate(zip(before, after, strict=True)) if a != b}


def stage_image() -> tuple[Image, dict[str, int]]:
    """ステージの実行ファイルの形: Location、xlc_log のイベント、スパン、enabled!"""
    image = Image()
    at = {
        "location": image.location(b"crates/a/src/x.rs", 120, 17),
        "event": image.metadata(b"event crates/a/src/x.rs", b"crates/a/src/x.rs", 40),
        "span": image.metadata(b"sweep", b"crates/a/src/x.rs", 77),
        "enabled": image.metadata(
            b"enabled crates/a/src/y.rs:853", b"crates/a/src/y.rs", 853
        ),
    }
    return image, at


class NormalizeTest(unittest.TestCase):
    def test_zeroes_the_line_and_column_of_locations(self):
        image, at = stage_image()
        data, counts = norm.normalize(image.build())
        self.assertEqual(
            (u32(data, at["location"] + 16), u32(data, at["location"] + 20)), (0, 0)
        )
        self.assertEqual(counts["location"], 1)

    def test_zeroes_the_line_of_event_span_and_enabled_metadata(self):
        image, at = stage_image()
        data, counts = norm.normalize(image.build())
        for kind in ("event", "span", "enabled"):
            self.assertEqual(u32(data, at[kind]), 1, kind)
            self.assertEqual(u32(data, at[kind] + 4), 0, kind)
        self.assertEqual(counts["metadata"], 3)

    def test_zeroes_the_line_digits_of_enabled_names(self):
        image, _ = stage_image()
        data, counts = norm.normalize(image.build())
        self.assertIn(b"enabled crates/a/src/y.rs:000x", data)
        self.assertNotIn(b":853", data)
        self.assertEqual(counts["enabled"], 1)

    def test_changes_no_other_byte(self):
        image, at = stage_image()
        before = image.build()
        data, _ = norm.normalize(before)
        expected = {RELRO_OFFSET + at["location"] + i for i in range(16, 24)}
        for kind in ("event", "span", "enabled"):
            expected |= {RELRO_OFFSET + at[kind] + i for i in range(4, 8)}
        name = before.index(b"enabled crates/a/src/y.rs:853")
        expected |= set(range(name + 26, name + 29))
        self.assertTrue(changed(before, data) <= expected)
        self.assertEqual(len(data), len(before))

    def test_same_code_on_shifted_lines_gives_the_same_copy(self):
        first, _ = stage_image()
        shifted = Image()
        shifted.location(b"crates/a/src/x.rs", 121, 17)
        shifted.metadata(b"event crates/a/src/x.rs", b"crates/a/src/x.rs", 41)
        shifted.metadata(b"sweep", b"crates/a/src/x.rs", 78)
        shifted.metadata(b"enabled crates/a/src/y.rs:854", b"crates/a/src/y.rs", 854)
        self.assertEqual(
            norm.normalize(first.build())[0], norm.normalize(shifted.build())[0]
        )

    def test_text_without_nul_or_with_zero_column_is_not_a_location(self):
        image, _ = stage_image()
        ptr, length = image.string(b"crates/a/src/z.rs")
        plain = image.place(struct.pack("<QQII", ptr, length, 5, 6))
        ptr, length = image.string(b"crates/a/src/w.rs", nul=True)
        zero = image.place(struct.pack("<QQII", ptr, length, 5, 0))
        data, counts = norm.normalize(image.build())
        self.assertEqual((u32(data, plain + 16), u32(data, plain + 20)), (5, 6))
        self.assertEqual((u32(data, zero + 16), u32(data, zero + 20)), (5, 0))
        self.assertEqual(counts["location"], 1)

    def test_location_in_the_last_slot_of_the_section_is_zeroed(self):
        image, _ = stage_image()
        ptr, length = image.string(b"crates/a/src/z.rs", nul=True)
        image.relro += b"\x00" * (SIZE - 24 - len(image.relro))
        last = image.place(struct.pack("<QQII", ptr, length, 9, 3))
        self.assertEqual(last + 24, SIZE)
        data, counts = norm.normalize(image.build())
        self.assertEqual((u32(data, last + 16), u32(data, last + 20)), (0, 0))
        self.assertEqual(counts["location"], 2)

    def test_metadata_whose_strings_end_the_section_is_zeroed(self):
        image, _ = stage_image()
        # スパンの名前を .rodata の最後に置く（後ろに読めるバイトが無い）
        span = image.metadata(b"sweep", b"crates/a/src/x.rs", 77)
        ptr, length = image.last_string(b"coarse")
        image.relro[span + 8 : span + 24] = struct.pack("<QQ", ptr, length)
        data, counts = norm.normalize(image.build())
        self.assertEqual(u32(data, span + 4), 0)
        self.assertEqual(counts["metadata"], 4)

    def test_metadata_file_followed_by_nul_keeps_the_fields_after_it(self):
        image, _ = stage_image()
        # スパンだけのファイル。file!() の文字列の後ろが偶然 NUL で、後ろのフィールドが 0 でない
        span = image.metadata(b"homogenize", b"crates/a/src/only_spans.rs", 12)
        ptr, length = image.string(b"crates/a/src/only_spans.rs", nul=True)
        file = span + 8 + FILE_DISTANCE
        image.relro[file : file + 16] = struct.pack("<QQ", ptr, length)
        image.relro[file + 16 : file + 24] = struct.pack("<II", 3, 4)
        before = image.build()
        data, counts = norm.normalize(before)
        self.assertEqual(u32(data, span + 4), 0)
        self.assertEqual((u32(data, file + 16), u32(data, file + 20)), (3, 4))
        self.assertEqual(counts["location"], 1)

    def test_stage_event_named_by_line_is_an_error(self):
        image, _ = stage_image()
        image.metadata(b"event crates/a/src/x.rs:12", b"crates/a/src/x.rs", 12)
        with self.assertRaisesRegex(ValueError, "xlc_log"):
            norm.normalize(image.build())

    def test_dependency_event_named_by_line_is_kept_apart_from_the_rule(self):
        image, _ = stage_image()
        dep = image.metadata(
            b"event /registry/dep/src/lib.rs:9", b"/registry/dep/src/lib.rs", 9
        )
        data, _ = norm.normalize(image.build())
        self.assertEqual(u32(data, dep + 4), 0)
        self.assertIn(b"/registry/dep/src/lib.rs:0x", data)

    def test_name_not_after_some_line_is_an_error(self):
        image, at = stage_image()
        image.relro[at["event"] : at["event"] + 4] = struct.pack("<I", 0)
        with self.assertRaisesRegex(ValueError, "Some"):
            norm.normalize(image.build())

    def test_enabled_name_and_line_that_disagree_are_an_error(self):
        image, _ = stage_image()
        image.metadata(b"enabled crates/a/src/y.rs:10", b"crates/a/src/y.rs", 11)
        with self.assertRaisesRegex(ValueError, "line 11"):
            norm.normalize(image.build())

    def test_metadata_whose_file_is_at_another_distance_is_an_error(self):
        image, _ = stage_image()
        other = image.metadata(b"event crates/a/src/v.rs", b"crates/a/src/v.rs", 3)
        image.relro[other + 8 + FILE_DISTANCE : other + 8 + FILE_DISTANCE + 16] = (
            b"\x00" * 16
        )
        ptr, length = image.string(b"crates/a/src/v.rs")
        image.relro[other + 0x68 : other + 0x78] = struct.pack("<QQ", ptr, length)
        with self.assertRaisesRegex(ValueError, "距離"):
            norm.normalize(image.build())

    def test_binary_without_events_is_an_error(self):
        image = Image()
        image.location(b"crates/a/src/x.rs", 1, 1)
        with self.assertRaisesRegex(ValueError, "距離"):
            norm.normalize(image.build())

    def test_binary_without_locations_is_an_error(self):
        image = Image()
        image.metadata(b"event crates/a/src/x.rs", b"crates/a/src/x.rs", 40)
        with self.assertRaisesRegex(ValueError, "Location"):
            norm.normalize(image.build())

    def test_not_elf64_little_endian_is_an_error(self):
        data = bytearray(stage_image()[0].build())
        data[5] = 2
        with self.assertRaisesRegex(ValueError, "ELF64"):
            norm.normalize(bytes(data))


class MainTest(unittest.TestCase):
    def test_writes_the_copy_with_the_mode_of_the_built_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp, "built"), Path(tmp, "bin")
            src.write_bytes(stage_image()[0].build())
            src.chmod(0o755)
            result = subprocess.run(
                [sys.executable, str(Path(norm.__file__)), str(src), str(dst)],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertEqual(dst.read_bytes(), norm.normalize(src.read_bytes())[0])
            self.assertEqual(os.stat(dst).st_mode & 0o777, 0o755)
            self.assertIn("location=1 metadata=3 enabled=1", result.stdout)


if __name__ == "__main__":
    unittest.main()
