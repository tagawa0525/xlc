"""rod_table.py の結合（制御棒の値ごとのテーブルを、制御棒の軸を持つ 1 つのテーブルに）。

実行: python3 -m unittest discover -s examples/core -p 'test_rod_table.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rod_table


def table(k: float, **extra) -> dict:
    """燃焼度 2 点の軸だけのテーブル（状態は k∞ だけを持つ）"""
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "axes": [{"name": "burnup", "points": [0.0, 10.0]}],
        "split": 1,
        "states": [{"k_inf": k}, {"k_inf": k / 2}],
        **extra,
    }


def merge(*tables: dict) -> dict:
    """値 out、aic、… のテーブルを書いて rod_table.py で結合した結果"""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        args = []
        for value, t in zip(["out", "aic", "b4c"], tables):
            path = root / f"{value}.msgpack"
            path.write_bytes(msgpack.packb(t, use_bin_type=True))
            args.append(f"{value}={path}")
        output = root / "merged.msgpack"
        argv = sys.argv
        sys.argv = ["rod_table.py", str(output), *args]
        try:
            rod_table.main()
        finally:
            sys.argv = argv
        return msgpack.unpackb(output.read_bytes())


class Merge(unittest.TestCase):
    def test_states_interleave_with_the_rod_axis_last(self):
        merged = merge(table(1.2), table(1.0))
        self.assertEqual(merged["axes"][-1], {"name": "rod", "values": ["out", "aic"]})
        self.assertEqual([s["k_inf"] for s in merged["states"]], [1.2, 1.0, 0.6, 0.5])

    def test_the_heavy_metal_mass_is_kept_when_every_table_has_the_same(self):
        merged = merge(table(1.2, heavy_metal=1258.0), table(1.0, heavy_metal=1258.0))
        self.assertEqual(merged["heavy_metal"], 1258.0)

    def test_tables_with_other_heavy_metal_masses_are_errors(self):
        with self.assertRaises(SystemExit):
            merge(table(1.2, heavy_metal=1258.0), table(1.0, heavy_metal=1250.0))
        with self.assertRaises(SystemExit):
            merge(table(1.2), table(1.0, heavy_metal=1258.0))


if __name__ == "__main__":
    unittest.main()
