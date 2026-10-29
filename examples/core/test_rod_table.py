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
        "points": [[0], [1]],
        "states": [{"k_inf": k}, {"k_inf": k / 2}],
        **extra,
    }


def boron(k: float, nominal: float = 1300.0) -> dict:
    """燃焼度 2 点 × ホウ素 2 点（公称は `nominal`）のテーブル"""
    return {
        "format": "xlc-table",
        "bounds": [1e7, 0.625, 1e-5],
        "axes": [
            {"name": "burnup", "points": [0.0, 10.0]},
            {
                "name": "boron",
                "points": [1000.0, 1300.0],
                "b10": 0.199,
                "nominal": nominal,
                "transform": "linear",
            },
        ],
        "split": 1,
        "points": [[0, 0], [0, 1], [1, 0], [1, 1]],
        "states": [{"k_inf": k + i} for i in range(4)],
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
        self.assertEqual(merged["points"], [[0, 0], [0, 1], [1, 0], [1, 1]])

    def test_rodded_tables_may_solve_other_points_and_some_burnups(self):
        # 制御棒を入れた組（ADR-0200）は、燃焼度 10 の点の公称のホウ素の状態だけ。燃焼度の番号は棒のないテーブルの軸で数え直す
        rodded = boron(5.0)
        rodded["axes"][0]["points"] = [10.0]
        rodded["points"] = [[0, 1]]
        rodded["states"] = [{"k_inf": 0.7}]
        merged = merge(boron(1.0), rodded)
        self.assertEqual(merged["axes"][0]["points"], [0.0, 10.0])
        self.assertEqual(
            merged["points"],
            [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
        )
        self.assertEqual(
            [s["k_inf"] for s in merged["states"]], [1.0, 2.0, 3.0, 4.0, 0.7]
        )

    def test_burnups_outside_the_first_table_are_errors(self):
        rodded = table(1.0)
        rodded["axes"][0]["points"] = [0.0, 20.0]
        with self.assertRaises(SystemExit):
            merge(table(1.2), rodded)

    def test_axes_with_other_nominal_values_or_points_are_errors(self):
        with self.assertRaises(SystemExit):
            merge(boron(1.0), boron(1.0, nominal=1000.0))
        other = boron(1.0)
        other["axes"][1]["points"] = [900.0, 1300.0]
        with self.assertRaises(SystemExit):
            merge(boron(1.0), other)

    def test_tables_without_points_are_errors(self):
        bare = table(1.0)
        del bare["points"]
        with self.assertRaises(SystemExit):
            merge(table(1.2), bare)

    def test_the_heavy_metal_mass_is_kept_when_every_table_has_the_same(self):
        merged = merge(table(1.2, heavy_metal=1258.0), table(1.0, heavy_metal=1258.0))
        self.assertEqual(merged["heavy_metal"], 1258.0)

    def test_tables_with_other_heavy_metal_masses_are_errors(self):
        with self.assertRaises(SystemExit):
            merge(table(1.2, heavy_metal=1258.0), table(1.0, heavy_metal=1250.0))
        with self.assertRaises(SystemExit):
            merge(table(1.2), table(1.0, heavy_metal=1258.0))

    def test_micro_states_follow_the_points_and_keep_the_first_densities(self):
        # 微視的断面積（ADR-0216）は点ごとに 1 つで、状態と同じく並べ替える。数密度は燃焼度の点ごとで、制御棒の値によらない
        out = table(1.2, micro=micro(["o0", "o1"], [[7e-4], [6e-4]]))
        rodded = table(1.0, micro=micro(["r0", "r1"], [[7e-4], [6e-4]]))
        rodded["axes"][0]["points"] = [10.0]
        rodded["points"] = [[0]]
        rodded["states"] = [{"k_inf": 0.5}]
        rodded["micro"]["states"] = [["r1"]]
        rodded["micro"]["densities"] = [[6e-4]]
        merged = merge(out, rodded)
        self.assertEqual(merged["micro"]["nuclides"], ["U-235"])
        self.assertEqual(merged["micro"]["densities"], [[7e-4], [6e-4]])
        self.assertEqual(merged["micro"]["states"], [["o0"], ["o1"], ["r1"]])

    def test_micro_data_must_agree_across_the_tables(self):
        def with_micro(k, densities, nuclides=("U-235",)):
            return table(k, micro=micro(["s0", "s1"], densities, list(nuclides)))

        # 片方だけが微視的断面積を持つ
        with self.assertRaises(SystemExit):
            merge(with_micro(1.2, [[7e-4], [6e-4]]), table(1.0))
        with self.assertRaises(SystemExit):
            merge(table(1.2), with_micro(1.0, [[7e-4], [6e-4]]))
        # 追う核種が違う
        with self.assertRaises(SystemExit):
            merge(
                with_micro(1.2, [[7e-4], [6e-4]]),
                with_micro(1.0, [[7e-4], [6e-4]], ("U-238",)),
            )
        # 数密度の NaN は、どちらの側にあっても合わないものとして止める
        with self.assertRaises(SystemExit):
            merge(
                with_micro(1.2, [[7e-4], [6e-4]]),
                with_micro(1.0, [[7e-4], [float("nan")]]),
            )
        with self.assertRaises(SystemExit):
            merge(
                with_micro(1.2, [[7e-4], [float("nan")]]),
                with_micro(1.0, [[7e-4], [6e-4]]),
            )
        # 同じ燃焼度の数密度が違う（同じ履歴なら丸めの外で違わない）
        with self.assertRaises(SystemExit):
            merge(
                with_micro(1.2, [[7e-4], [6e-4]]), with_micro(1.0, [[7e-4], [6.1e-4]])
            )


def micro(states: list, densities: list, nuclides: list | None = None) -> dict:
    """1 核種の微視的断面積（状態は名前の文字列で代える）と、燃焼度の点ごとの数密度"""
    return {
        "nuclides": nuclides or ["U-235"],
        "densities": densities,
        "states": [[s] for s in states],
    }


if __name__ == "__main__":
    unittest.main()
