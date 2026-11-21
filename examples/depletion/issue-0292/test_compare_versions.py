"""compare_versions.py（#0292、#0303: 核データの版ごとの全核種の燃焼の比較）。

実行: python3 -m unittest discover -s examples/depletion/issue-0292 -p 'test_compare_versions.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_versions


def capture(sigma):
    return {"fission": 0.0, "reactions": [{"mt": 102, "channels": [sigma]}]}


# 2 群の状態: 散乱が無く、核分裂の中性子はすべて速い群に出るので、中性子束は [1, 0]、P = νΣf,1
def state(k, nu_fission=2.0):
    return {
        "k_inf": k,
        "total": [1.0, 1.0],
        "scatter": [[0.0, 0.0], [0.0, 0.0]],
        "chi": [1.0, 0.0],
        "absorption": [1.0, 0.0],
        "nu_fission": [nu_fission, 0.0],
    }


def history(names, initial, later, ks, transportable=None):
    """燃焼度 0、10、30、60 の点を持つ、体積 0.5 の燃焼の領域 1 つの履歴（燃焼度 0 の後は同じ組成 `later`）"""
    return {
        "format": "xlc-depletion-history",
        "nuclides": names,
        "transportable": transportable or names,
        "regions": [{"material": "fuel", "volume": 0.5}],
        "points": [
            {"burnup": b, "k_inf": k, "densities": [initial if b == 0 else later]}
            for b, k in zip([0.0, 10.0, 30.0, 60.0], ks)
        ],
    }


def table(names, sigmas, ks, nu_fission=2.0):
    """燃焼度 0、10、30、60 の基準条件の 1 状態ずつのテーブル（冷却計算の捕獲の断面積はどの点も `sigmas`）"""
    return {
        "format": "xlc-table",
        "axes": [{"name": "burnup", "points": [0.0, 10.0, 30.0, 60.0]}],
        "points": [[0], [1], [2], [3]],
        "states": [state(k, nu_fission) for k in ks],
        "cooling": {
            "nuclides": names,
            "points": [[capture(s) for s in sigmas] for _ in range(4)],
        },
    }


BASE_NAMES = ["U-235", "Co-58", "Sm-149"]


def versions():
    """基準の版 B と、Co-58 の捕獲が大きく、チェーンに X-1 を持つ版 A"""
    b = history(BASE_NAMES, [1e-3, 0.0, 0.0], [1e-3, 1e-8, 1e-6], [1.3, 1.2, 1.1, 1.0])
    a = history(
        ["X-1", *BASE_NAMES],
        [0.0, 1e-3, 0.0, 0.0],
        [5e-9, 1e-3, 2e-8, 1e-6],
        [1.3, 1.19, 1.08, 0.97],
    )
    tb = table(BASE_NAMES, [10.0, 100.0, 1000.0], [1.3, 1.2, 1.1, 1.0])
    ta = table(
        ["X-1", *BASE_NAMES], [50.0, 20.0, 1000.0, 1000.0], [1.3, 1.19, 1.08, 0.97]
    )
    return {"B": (b, tb, 1.0), "A": (a, ta, 1.0)}


class Reactivity(unittest.TestCase):
    def test_the_difference_to_the_base_is_given_at_every_point(self):
        s = compare_versions.summarize(versions(), "B")
        a = s["versions"]["A"]
        self.assertEqual(s["burnup"], [0.0, 10.0, 30.0, 60.0])
        self.assertAlmostEqual(a["reactivity_pcm"][3], 1e5 * (1 / 1.0 - 1 / 0.97))
        self.assertNotIn("reactivity_pcm", s["versions"]["B"])

    def test_the_marks_are_the_record_burnups(self):
        s = compare_versions.summarize(versions(), "B")
        self.assertEqual([m["burnup"] for m in s["marks"]], [0.0, 10.0, 30.0, 60.0])

    def test_histories_of_other_burnups_or_compositions_stop(self):
        v = versions()
        v["A"][0]["points"][1]["burnup"] = 11.0
        with self.assertRaises(SystemExit):
            compare_versions.summarize(v, "B")
        v = versions()
        v["A"][0]["points"][0]["densities"][0][1] = 2e-3
        with self.assertRaises(SystemExit):
            compare_versions.summarize(v, "B")


class Densities(unittest.TestCase):
    def test_every_nuclide_of_a_chain_has_its_density_at_the_marks(self):
        s = compare_versions.summarize(versions(), "B")
        co = s["densities"]["Co-58"]
        # 燃料のピン全体の平均（領域が 1 つなので領域の値）
        self.assertAlmostEqual(co["B"][3], 1e-8)
        self.assertAlmostEqual(co["A"][3], 2e-8)
        self.assertIsNone(s["densities"]["X-1"]["B"][3])
        self.assertAlmostEqual(s["densities"]["X-1"]["A"][3], 5e-9)

    def test_nuclides_without_a_density_in_any_version_are_left_out(self):
        s = compare_versions.summarize(versions(), "B")
        self.assertNotIn("Sm-150", s["densities"])


class Split(unittest.TestCase):
    def test_a_term_is_the_change_of_the_absorption_per_production(self):
        # a_j = N̄ⱼ σₐⱼ φ̄ / P。N̄ は体積の比 0.5 / 1.0 を掛けたもの、φ̄ = 1、P = 2
        s = compare_versions.summarize(versions(), "B")
        p = s["versions"]["A"]["split"][3]
        self.assertEqual(p["burnup"], 60.0)
        co = p["nuclides"]["Co-58"]
        self.assertAlmostEqual(
            co["pcm"], 1e5 * (1e-8 * 0.5 * 100.0 - 2e-8 * 0.5 * 1000.0) / 2.0
        )
        # 基準の版のテーブルに無い核種は、基準の吸収を 0 とする
        x = p["nuclides"]["X-1"]
        self.assertAlmostEqual(x["pcm"], -1e5 * 5e-9 * 0.5 * 50.0 / 2.0)

    def test_the_remainder_closes_the_direct_difference_of_the_tables(self):
        s = compare_versions.summarize(versions(), "B")
        for p in s["versions"]["A"]["split"]:
            terms = sum(v["pcm"] for v in p["nuclides"].values())
            self.assertAlmostEqual(p["remainder_pcm"], p["direct_pcm"] - terms)
            self.assertAlmostEqual(
                p["direct_pcm"],
                1e5
                * (
                    1 / s["versions"]["B"]["table_k_inf"][p["index"]]
                    - 1 / s["versions"]["A"]["table_k_inf"][p["index"]]
                ),
            )

    def test_the_remainder_is_the_absorption_outside_the_terms_and_the_balance(self):
        # 残り = 格子の吸収 / 生成のうち燃料の核種の項の外の差 + 1/k と格子の吸収 / 生成の差。テーブルの状態は
        # 吸収 / 生成 = 1 / 2 で、k はそれと別の値なので、2 つ目が 0 にならない
        s = compare_versions.summarize(versions(), "B")
        for p in s["versions"]["A"]["split"]:
            terms = sum(v["pcm"] for v in p["nuclides"].values())
            # 2 つの版で吸収 / 生成が同じなので、項の外の差は項の和の符号を変えたもの
            self.assertAlmostEqual(p["outside_pcm"], -terms)
            self.assertAlmostEqual(p["balance_pcm"], p["direct_pcm"])
            self.assertAlmostEqual(
                p["remainder_pcm"], p["outside_pcm"] + p["balance_pcm"]
            )

    def test_the_report_gives_the_parts_of_the_remainder(self):
        s = compare_versions.summarize(versions(), "B")
        s["inputs"] = []
        line = next(
            x for x in compare_versions.report(s) if x.startswith("A − B の 60 MWd/kg")
        )
        p = s["versions"]["A"]["split"][3]
        self.assertIn(f"項の外の吸収 {p['outside_pcm']:+.1f}", line)
        self.assertIn(f"1/k と吸収 / 生成の差 {p['balance_pcm']:+.1f}", line)

    def test_the_terms_of_nuclides_in_one_table_only_are_summed(self):
        # 評価の有無の分の 1 次の目安: 片方の版のテーブルにしか無い核種（X-1 は版 A だけ）の項の和
        s = compare_versions.summarize(versions(), "B")
        p = s["versions"]["A"]["split"][3]
        self.assertAlmostEqual(p["only_version_pcm"], p["nuclides"]["X-1"]["pcm"])
        self.assertEqual(p["only_base_pcm"], 0.0)

    def test_the_report_gives_the_terms_of_one_version_only(self):
        s = compare_versions.summarize(versions(), "B")
        s["inputs"] = []
        line = next(
            x for x in compare_versions.report(s) if x.startswith("A − B の 60 MWd/kg")
        )
        p = s["versions"]["A"]["split"][3]
        self.assertIn(f"A だけの核種 {p['only_version_pcm']:+.1f}", line)
        self.assertIn(f"B だけの核種 {p['only_base_pcm']:+.1f}", line)

    def test_the_burnup_zero_point_is_split_too(self):
        # 燃焼度 0 の組成は同じで、版の差は断面積の差だけ
        s = compare_versions.summarize(versions(), "B")
        p = s["versions"]["A"]["split"][0]
        self.assertEqual(p["burnup"], 0.0)
        u = p["nuclides"]["U-235"]
        self.assertAlmostEqual(u["pcm"], 1e5 * 1e-3 * 0.5 * (10.0 - 20.0) / 2.0)


class Guards(unittest.TestCase):
    def test_a_lattice_volume_that_makes_the_fuel_absorb_more_than_the_lattice_stops(
        self,
    ):
        h, t, _ = versions()["A"]
        with self.assertRaises(SystemExit):
            compare_versions.absorptions(h, t, 1e-6)

    def test_a_history_with_two_points_at_a_table_burnup_stops(self):
        h, t, volume = versions()["A"]
        h["points"].append({**h["points"][-1]})
        with self.assertRaises(SystemExit):
            compare_versions.absorptions(h, t, volume)

    def places(
        self, root, library="../v/library/depletion.msgpack", chain="../c.msgpack"
    ):
        """版の置き場 v（library と燃焼の履歴 u）と、テーブルの入力 ufg/all.toml"""
        (root / "v" / "u").mkdir(parents=True)
        (root / "ufg" / "out").mkdir(parents=True)
        (root / "v" / "u" / "nominal.toml").write_text(
            '[depletion]\nchain = "../../c.msgpack"\n'
        )
        (root / "ufg" / "all.toml").write_text(
            'few-group = { output = "out/all.msgpack" }\n'
            f'library = "{library}"\n'
            f'[history]\nfile = "../v/u/nominal.msgpack"\nchain = "{chain}"\n'
        )
        return (
            root / "v" / "u" / "nominal.msgpack",
            root / "ufg" / "out" / "all.msgpack",
            root / "ufg" / "all.toml",
        )

    def test_the_table_input_must_use_the_library_and_chain_of_the_history(self):
        for library, chain in [
            ("../w/library/depletion.msgpack", "../c.msgpack"),
            ("../v/library/depletion.msgpack", "../d.msgpack"),
        ]:
            with (
                self.subTest(library=library, chain=chain),
                tempfile.TemporaryDirectory() as d,
            ):
                files = self.places(Path(d), library, chain)
                with self.assertRaises(SystemExit):
                    compare_versions.check_inputs(*files)

    def test_the_shares_check_the_inputs_too(self):
        with tempfile.TemporaryDirectory() as d:
            _, table, table_input = self.places(Path(d))
            names = Path(d) / "names.txt"
            names.write_text("Co-58\n")
            argv = [
                "compare_versions.py",
                "shares",
                "--json",
                str(Path(d) / "o.json"),
                str(Path(d) / "other.msgpack"),
                str(table),
                str(table_input),
                str(names),
            ]
            with mock.patch.object(sys, "argv", argv), self.assertRaises(SystemExit):
                compare_versions.main()

    def test_the_table_input_must_solve_the_given_history_into_the_given_table(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            history, table, table_input = self.places(root)
            compare_versions.check_inputs(history, table, table_input)
            with self.assertRaises(SystemExit):
                compare_versions.check_inputs(
                    root / "other.msgpack", table, table_input
                )
            with self.assertRaises(SystemExit):
                compare_versions.check_inputs(history, root / "x.msgpack", table_input)


class Shares(unittest.TestCase):
    def test_the_shares_take_the_listed_nuclides_of_the_version(self):
        h, t, volume = versions()["A"]
        out = compare_versions.version_shares(h, t, volume, ["Co-58"])
        p = out["points"][-1]
        self.assertEqual(p["burnup"], 60.0)
        # 格子の吸収 A = Σa,1 φ1 = 1 に占める割合
        self.assertAlmostEqual(
            p["untransported_all"]["Co-58"]["fraction"], 2e-8 * 0.5 * 1000.0
        )
        self.assertEqual(p["untransported_all"]["Co-58"]["captured_ratio"], 1.0)


class Sets(unittest.TestCase):
    def dilution(self, rows):
        return {"nuclides": [{"nuclide": n, "max_pcm": v} for n, v in rows]}

    def test_the_sets_of_a_version_are_compared_with_the_base_version(self):
        base = [self.dilution([("A", 0.5), ("B", 0.2)])]
        other = [self.dilution([("A", 0.3), ("C", 0.15), ("B", 0.05)])]
        d = compare_versions.set_difference(base, other, 0.1)
        self.assertEqual(d["base"], ["A", "B"])
        self.assertEqual(d["version"], ["A", "C"])
        self.assertEqual(d["only_base"], ["B"])
        self.assertEqual(d["only_version"], ["C"])


if __name__ == "__main__":
    unittest.main()
