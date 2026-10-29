"""compare_core_steps.py の要約（ステップの幅を変えた core の燃焼の差とその比、中点のホウ素の探索と案 S1 の差、固有値計算の回数）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_steps.py'
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_core_steps

# 基準のサイクル: (終わりの EFPD, 出力の割合)。0 は崩壊だけのステップ（1 日）
BASE = [(10.0, 1.0), (20.0, 0.5), (20.0, 0.0), (30.0, 1.0)]


def core(split, method="search", error=1.0, solves=3):
    """集合体 2 体・軸方向のノード 2 つ（上下で同じ値）の燃焼の履歴。出力のあるステップを `split` 個に分け、点の値は真の値に誤差 `error` / split² を足したもの
    （中点の形が 2 次の燃焼）"""
    points = []
    efpd, days = 0.0, 0.0

    def point(power):
        e = error / split**2
        b = efpd * 0.04
        return {
            "efpd": efpd,
            "days": days,
            **({} if power is None else {"power": power}),
            "state": "hfp",
            "burnup": b,
            "k": 1.0,
            "boron": 1000.0 - 10.0 * efpd + 2.0 * e,
            "solves": solves,
            "assembly_burnups": [b + e, b - e],
            "burnups": [[b + e] * 2, [b - e] * 2],
            "powers": [[0.25 * (1.0 + 0.01 * e)] * 2, [0.25 * (1.0 - 0.01 * e)] * 2],
        }

    points.append(point(None))
    for end, power in BASE:
        if power == 0.0:
            days += 1.0
            points.append(point(0.0))
            continue
        start = efpd
        for i in range(1, split + 1):
            efpd = start + (end - start) * i / split
            days += (end - start) / split / power
            points.append(point(power))
    return {
        "format": "xlc-core-history",
        "specific_power": 40.0,
        "states": {"hfp": {"search": {"boron": [0.0, 1300.0]}, "b10": 0.199}},
        "axial_mesh": [0.0, 10.0, 20.0],
        "assemblies": ["A1", "A2"],
        "heavy_metal": [[2000.0, 2000.0], [2000.0, 2000.0]],
        "midpoint_boron": method,
        "points": points,
    }


def log(core, solves=5):
    """燃焼のログ: 出力のあるステップごとの中点の行"""
    steps = sum(
        1 for a, b in zip(core["points"], core["points"][1:]) if b["efpd"] != a["efpd"]
    )
    return "".join(
        f"2026-10-29T00:00:00Z  INFO depletion step{{efpd={i}}}: xlc_core::depletion: depletion midpoint "
        f"k=1.0 boron=900.0 solves={solves}\n"
        for i in range(steps)
    )


def summarize(runs):
    return compare_core_steps.summarize(
        {k: c for k, c in runs.items()},
        {k: log(c, 1 if k[0] == "start" else 5) for k, c in runs.items()},
    )


def search_runs(splits=(1, 2, 4), **kwargs):
    return {("search", s): core(s, **kwargs) for s in splits}


class Steps(unittest.TestCase):
    def test_second_order_differences_shrink_by_four(self):
        s = summarize(search_runs())
        for p in s["steps"]["search"][1:]:
            (r,) = p["ratios"]
            self.assertEqual(r["splits"], [1, 4])
            for q in ("boron", "burnup"):
                self.assertAlmostEqual(r[q], 4.0, places=9)
            # 出力の割合は細かい方の値で割った相対差なので、比は 4 から誤差の大きさの分ずれる
            self.assertAlmostEqual(r["power"], 4.0, delta=0.01)
            first = p["differences"][0]
            self.assertEqual(first["splits"], [1, 2])
            # 粗い − 細かい: 2 (1 − 1/4)
            self.assertAlmostEqual(first["boron"], 1.5)
            self.assertAlmostEqual(first["burnup"], 0.75)

    def test_points_are_matched_by_efpd_and_days(self):
        s = summarize(search_runs())
        points = s["steps"]["search"]
        self.assertEqual(
            [(p["efpd"], p["days"]) for p in points],
            [(0.0, 0.0), (10.0, 10.0), (20.0, 30.0), (20.0, 31.0), (30.0, 41.0)],
        )
        # 崩壊だけのステップの点も比べる
        self.assertEqual(points[3]["power"], 0.0)

    def test_boron_ratios_keep_the_sign_of_the_differences(self):
        # 隣り合う差の符号が変われば比は負（漸近の範囲の外）
        runs = search_runs()
        runs[("search", 4)] = core(4, error=8.0)
        s = summarize(runs)
        self.assertLess(s["steps"]["search"][1]["ratios"][0]["boron"], 0.0)

    def test_equal_points_give_no_ratio(self):
        s = summarize(search_runs(error=0.0))
        self.assertIsNone(s["steps"]["search"][1]["ratios"][0]["boron"])

    def test_start_minus_search_at_the_splits_of_both(self):
        runs = search_runs(splits=(1, 2, 4, 8))
        runs.update({("start", sp): core(sp, "start") for sp in (1, 2)})
        s = summarize(runs)
        self.assertEqual(s["splits"], {"search": [1, 2, 4, 8], "start": [1, 2]})
        self.assertEqual(sorted(s["midpoint"]), ["1", "2"])
        self.assertEqual(len(s["steps"]["search"][1]["ratios"]), 2)
        self.assertEqual(s["steps"]["start"][1]["ratios"], [])

    def test_start_minus_search_at_each_split(self):
        runs = search_runs()
        runs.update({("start", sp): core(sp, "start", error=3.0) for sp in (1, 2, 4)})
        s = summarize(runs)
        self.assertEqual(s["methods"], ["search", "start"])
        p = s["midpoint"]["1"][1]
        self.assertAlmostEqual(p["boron"], 2.0 * (3.0 - 1.0))
        self.assertAlmostEqual(p["burnup"], 2.0)
        self.assertEqual(sorted(s["midpoint"]), ["1", "2", "4"])

    def test_solves_add_the_points_and_the_logged_midpoints(self):
        runs = search_runs()
        runs[("start", 1)] = core(1, "start")
        runs[("start", 2)] = core(2, "start")
        runs[("start", 4)] = core(4, "start")
        s = summarize(runs)
        by = {(x["method"], x["split"]): x for x in s["solves"]}
        # 基準: 点 5 つ × 3、出力のあるステップ 3 つの中点 × 5（start は 1）
        self.assertEqual(by[("search", 1)]["point_solves"], 15)
        self.assertEqual(by[("search", 1)]["midpoint_solves"], 15)
        self.assertEqual(by[("start", 1)]["midpoint_solves"], 3)
        self.assertEqual(by[("search", 4)]["points"], 14)
        self.assertEqual(by[("search", 4)]["total"], 14 * 3 + 12 * 5)


def axial_core(split, offsets, method="search"):
    """`core` の点 i の、上（10〜20 cm）と下（0〜10 cm）のノードの出力の割合の差を `offsets(i)` にした履歴"""
    c = core(split, method)
    for i, p in enumerate(c["points"]):
        a = offsets(i)
        p["powers"] = [[0.25 * (1.0 - a), 0.25 * (1.0 + a)]] * 2
    return c


class Axial(unittest.TestCase):
    def test_offsets_of_every_point_are_the_upper_minus_the_lower_half(self):
        runs = {("search", s): axial_core(s, lambda i: 0.1 * (-1) ** i) for s in (1, 2)}
        s = summarize(runs)
        by = {(x["method"], x["split"]): x for x in s["axial"]}
        self.assertEqual(sorted(by), [("search", 1), ("search", 2)])
        fine = by[("search", 2)]
        points = runs[("search", 2)]["points"]
        # 点ごとの日数（1e-6 日に丸める）と偏り（記録の表の桁の 1e-4 に丸める。要約を数十 KB に収める）
        self.assertEqual(fine["days"], [round(p["days"], 6) for p in points])
        self.assertEqual(
            fine["offset"], [round(0.1 * (-1) ** i, 4) for i in range(len(points))]
        )
        self.assertEqual(set(fine), {"method", "split", "days", "offset"})

    def test_the_mesh_needs_a_boundary_at_half_height(self):
        runs = {("search", s): axial_core(s, lambda i: 0.0) for s in (1, 2)}
        for c in runs.values():
            c["axial_mesh"] = [0.0, 12.0, 20.0]
        with self.assertRaises(SystemExit):
            summarize(runs)


class Methods(unittest.TestCase):
    def test_equilibrium_xenon_histories_are_their_own_method(self):
        runs = search_runs()
        runs.update({("start", sp): core(sp, "start") for sp in (1, 2)})
        runs.update({("equilibrium", sp): core(sp, error=0.5) for sp in (1, 2, 4)})
        s = summarize(runs)
        self.assertEqual(s["methods"], ["equilibrium", "search", "start"])
        self.assertEqual(s["splits"]["equilibrium"], [1, 2, 4])
        # 探索と案 S1 の差は、平衡 Xe の燃焼があっても出す
        self.assertEqual(sorted(s["midpoint"]), ["1", "2"])
        r = s["steps"]["equilibrium"][1]["ratios"][0]
        self.assertAlmostEqual(r["burnup"], 4.0, places=9)

    def test_methods_are_named_by_the_midpoint_boron_and_the_xenon(self):
        name = compare_core_steps.method_name
        micro = {"nuclides": ["Xe-135"]}
        # ミクロ燃焼しない履歴は Xe の扱いを持たない
        self.assertEqual(name({"midpoint_boron": "search"}), "search")
        self.assertEqual(
            name({"midpoint_boron": "search", "xenon": "tracked", **micro}), "search"
        )
        self.assertEqual(
            name({"midpoint_boron": "start", "xenon": "tracked", **micro}), "start"
        )
        self.assertEqual(
            name({"midpoint_boron": "search", "xenon": "equilibrium", **micro}),
            "equilibrium",
        )
        for bad in (
            {"midpoint_boron": "start", "xenon": "equilibrium", **micro},
            {"midpoint_boron": "search", "xenon": "other", **micro},
            {"midpoint_boron": "other"},
            # ミクロ燃焼の履歴は Xe の扱いを持ち、ミクロ燃焼しない履歴は持たない（ADR-0219）
            {"midpoint_boron": "search", **micro},
            {"midpoint_boron": "search", "xenon": "tracked"},
        ):
            with self.assertRaises(SystemExit):
                name(bad)


class Errors(unittest.TestCase):
    def test_runs_need_the_base_split(self):
        with self.assertRaises(SystemExit):
            summarize(search_runs(splits=(2, 4)))

    def test_runs_must_share_the_core(self):
        runs = search_runs()
        runs[("search", 2)]["heavy_metal"] = [[2000.0, 2000.0], [2000.0, 2001.0]]
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_runs_must_track_the_same_nuclides(self):
        runs = search_runs()
        runs[("search", 2)]["nuclides"] = ["Xe-135"]
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_finer_runs_hold_the_base_points(self):
        runs = search_runs()
        del runs[("search", 2)]["points"][2]
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_logs_hold_one_midpoint_per_efpd_step(self):
        runs = search_runs()
        logs = {k: log(c) for k, c in runs.items()}
        logs[("search", 2)] = logs[("search", 1)]
        with self.assertRaises(SystemExit):
            compare_core_steps.summarize(runs, logs)

    def test_histories_must_hold_the_points_of_their_split(self):
        # 分割 4 の履歴を分割 2 として渡すと、基準の点はそろうが、点の数が分割 2 のサイクルと合わない
        runs = search_runs()
        runs[("search", 2)] = core(4)
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_power_differences_stop_at_nodes_without_fine_power(self):
        runs = search_runs()
        runs[("search", 2)]["points"][2]["powers"] = [[0.5, 0.5], [0.0, 0.0]]
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_matched_points_must_share_the_state_and_power(self):
        runs = search_runs()
        runs[("search", 2)]["points"][2]["state"] = "rodded"
        with self.assertRaises(SystemExit):
            summarize(runs)

    def test_splits_must_double(self):
        with self.assertRaises(SystemExit):
            summarize(search_runs(splits=(1, 2, 8)))

    def test_start_needs_the_base_split(self):
        runs = search_runs()
        runs[("start", 2)] = core(2, "start")
        with self.assertRaises(SystemExit):
            summarize(runs)


if __name__ == "__main__":
    unittest.main()
