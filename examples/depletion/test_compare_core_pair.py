"""compare_core_pair.py の要約（入力の 1 項目だけを変えた core の燃焼の差と、固有値計算の数、反復の回数、時間）。

実行: python3 -m unittest discover -s examples/depletion -p 'test_compare_core_pair.py'
"""

import sys
import tempfile
import unittest
from pathlib import Path

import msgpack

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_core_pair import ITEMS, compare, pair_runs, parse, report

GUESS = ITEMS["initial-guess"]
XENON_BORON = ITEMS["xenon-boron"]
ACCELERATION = ITEMS["acceleration"]


def core(k=1.0, boron=1000.0, burnup=1.0, power=0.5, solves=(1, 2)):
    """集合体 2 体・軸方向のノード 1 つの燃焼の履歴。点は 0 と 10 EFPD。2 体目のノードは核分裂がない"""
    return {
        "specific_power": 40.0,
        "states": {"hfp": {}},
        "axial_mesh": [0.0, 10.0],
        "assemblies": ["A", "R"],
        "heavy_metal": [[1.0], [0.0]],
        "points": [
            {
                "efpd": efpd,
                "days": efpd,
                "state": "hfp",
                "k": k,
                "boron": boron,
                "solves": n,
                "burnups": [[burnup * efpd], [0.0]],
                "powers": [[power], [0.0]],
            }
            for efpd, n in zip((0.0, 10.0), solves)
        ],
    }


def log(iterations, midpoints=((1.0, 900.0, 1),), total=5.0):
    """固有値計算ごとの反復の回数 `iterations`（1 回 0.25 s）、中点の行 (k, ホウ素, 回数)、全体の時間の core のログ"""
    lines = [
        f"2026-11-06T00:00:00Z  INFO nem{{iterations={i} seconds=0.25}}: xlc_core::nem: close time.busy=250ms time.idle=1µs"
        for i in iterations
    ]
    lines += [
        f"2026-11-06T00:00:00Z  INFO depletion step{{efpd=10.0}}: xlc_core::depletion: depletion midpoint k={k} boron={b} solves={n}"
        for k, b, n in midpoints
    ]
    lines.append(f"2026-11-06T00:00:00Z  INFO core: core seconds={total}")
    return "\n".join(lines) + "\n"


def write_run(
    directory: Path,
    name: str,
    top: str = "",
    depletion: str = "",
    history=None,
    text=None,
) -> Path:
    """core の入力 `<name>.toml` と、それが名前を書く履歴 `<name>.msgpack`、そのログ `<name>.log`。`top` と `depletion` は
    ケースと `[depletion]` に足す行、`text` はログ（省けば `log([1] * 4)`）"""
    source = directory / f"{name}.toml"
    source.write_text(
        f'cycle = "cycle.toml"\n{top}\n[depletion]\nspecific-power = 40.0\n'
        f'history = "{name}.msgpack"\n{depletion}\n'
    )
    (directory / f"{name}.msgpack").write_bytes(msgpack.packb(history or core()))
    (directory / f"{name}.log").write_text(log([1] * 4) if text is None else text)
    return source


class Inputs(unittest.TestCase):
    """2 本の燃焼の役は、引数の順でなく、入力の項目の値で決まる"""

    def test_reads_the_history_and_the_log_that_each_input_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            nested = write_run(
                tmp, "nested", depletion='xenon-boron = "nested"', history=core(k=0.9)
            )
            joint = write_run(tmp, "joint")
            base_core, base_log, other, _, files = pair_runs(XENON_BORON, nested, joint)
        self.assertEqual(base_core["points"][0]["k"], 0.9)
        self.assertEqual(other["points"][0]["k"], 1.0)
        self.assertEqual(base_log, log([1] * 4))
        self.assertEqual(
            [f.name for f in files],
            [
                "nested.toml",
                "nested.msgpack",
                "nested.log",
                "joint.toml",
                "joint.msgpack",
                "joint.log",
            ],
        )

    def test_inputs_given_in_the_other_order_are_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            nested = write_run(tmp, "nested", depletion='xenon-boron = "nested"')
            joint = write_run(tmp, "joint")
            with self.assertRaises(SystemExit) as stop:
                pair_runs(XENON_BORON, joint, nested)
        self.assertIn("xenon-boron", str(stop.exception))

    def test_inputs_with_the_same_value_of_the_item_are_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            a = write_run(tmp, "a")
            b = write_run(tmp, "b")
            with self.assertRaises(SystemExit):
                pair_runs(GUESS, a, b)

    def test_inputs_that_differ_in_another_item_are_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            flat = write_run(
                tmp,
                "flat",
                top='acceleration = "none"',
                depletion='initial-guess = "flat"',
            )
            default = write_run(tmp, "default")
            with self.assertRaises(SystemExit) as stop:
                pair_runs(GUESS, flat, default)
        self.assertIn("acceleration", str(stop.exception))

    def test_the_log_is_named_after_the_input(self):
        """記録のコマンドは core の標準エラーを `<入力の名前>.log` に書く"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            nested = write_run(tmp, "nested", depletion='xenon-boron = "nested"')
            joint = write_run(tmp, "joint")
            (tmp / "joint.toml").write_text(
                (tmp / "joint.toml")
                .read_text()
                .replace("joint.msgpack", "other.msgpack")
            )
            (tmp / "joint.msgpack").rename(tmp / "other.msgpack")
            *_, files = pair_runs(XENON_BORON, nested, joint)
        self.assertEqual(
            [f.name for f in files[3:]], ["joint.toml", "other.msgpack", "joint.log"]
        )

    def test_a_case_item_is_read_from_the_top_of_the_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            none = write_run(tmp, "none", top='acceleration = "none"')
            cmfd = write_run(tmp, "cmfd")
            pair_runs(ACCELERATION, none, cmfd)
            with self.assertRaises(SystemExit):
                pair_runs(ACCELERATION, cmfd, none)
            # `[depletion]` に書いた `acceleration` は core の入力の項目でない
            misplaced = write_run(tmp, "misplaced", depletion='acceleration = "none"')
            with self.assertRaises(SystemExit):
                pair_runs(ACCELERATION, misplaced, cmfd)


class InitialGuess(unittest.TestCase):
    def test_differences_are_the_default_minus_the_flat_start(self):
        flat = core()
        default = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4, burnup=1.0 + 1e-9)
        c = compare(
            GUESS,
            flat,
            log([100, 100, 100, 100], midpoints=((1.0, 900.0, 1),)),
            default,
            log(
                [100, 40, 30, 20], midpoints=((1.0 - 1e-9, 900.0 + 2e-4, 1),), total=4.0
            ),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(later["burnup"], 1e-9, delta=1e-15)
        # 核分裂のないノードは除く
        self.assertEqual(later["power"], 0.0)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(
            c["runs"]["default"],
            {
                "solves": 4,
                "iterations": 190,
                "eigenvalue_seconds": 1.0,
                "seconds": 4.0,
            },
        )
        self.assertEqual(c["runs"]["flat"]["iterations"], 400)

    def test_solve_counts_must_match_the_nem_spans(self):
        # 点 1 + 2 と中点 1 の 4 回に、span が 3 つしかない
        with self.assertRaises(SystemExit):
            compare(GUESS, core(), log([1, 1, 1]), core(), log([1, 1, 1]))

    def test_histories_of_other_points_are_errors(self):
        other = core()
        other["points"][1]["days"] = 11.0
        with self.assertRaises(SystemExit):
            compare(GUESS, core(), log([1] * 4), other, log([1] * 4))

    def test_power_in_one_history_only_is_an_error(self):
        other = core()
        other["points"][1]["powers"][1] = [0.1]
        with self.assertRaises(SystemExit):
            compare(GUESS, core(), log([1] * 4), other, log([1] * 4))


class XenonBoron(unittest.TestCase):
    def test_differences_are_the_joint_minus_the_nested_iteration(self):
        nested = core(solves=(1, 6))
        joint = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4, solves=(1, 3))
        c = compare(
            XENON_BORON,
            nested,
            log([100] * 8, midpoints=((1.0, 900.0, 1),), total=9.0),
            joint,
            log([100, 40, 30, 20, 10], midpoints=((1.0, 900.0 + 2e-4, 1),), total=4.0),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(sorted(c["runs"]), ["joint", "nested"])
        self.assertEqual(c["runs"]["nested"]["solves"], 8)
        self.assertEqual(
            c["runs"]["joint"],
            {
                "solves": 5,
                "iterations": 200,
                "eigenvalue_seconds": 1.25,
                "seconds": 4.0,
            },
        )

    def test_report_gives_the_solves_and_their_ratio(self):
        c = compare(
            XENON_BORON,
            core(solves=(1, 6)),
            log([100] * 8, total=9.0),
            core(solves=(1, 3)),
            log([50] * 5, total=4.0),
        )
        summary = {"format": XENON_BORON.format, "cases": {"small": c}}
        self.assertIn("固有値計算 8 → 5 回（0.625 倍）", "\n".join(report(summary)))


class Acceleration(unittest.TestCase):
    def test_differences_are_the_cmfd_minus_the_unaccelerated_solves(self):
        plain = core()
        cmfd = core(k=1.0 + 2e-9, boron=1000.0 - 1e-4)
        c = compare(
            ACCELERATION,
            plain,
            log([4000, 1200, 900, 3000], midpoints=((1.0, 900.0, 1),), total=9.0),
            cmfd,
            log([60, 20, 15, 40], midpoints=((1.0, 900.0 + 2e-4, 1),), total=1.5),
        )
        later = c["points"][1]
        self.assertAlmostEqual(later["k"], 2e-9, delta=1e-15)
        self.assertAlmostEqual(later["boron"], -1e-4, delta=1e-10)
        self.assertAlmostEqual(c["midpoints"][0]["boron"], 2e-4, delta=1e-10)
        self.assertEqual(sorted(c["runs"]), ["cmfd", "none"])
        self.assertEqual(
            c["runs"]["cmfd"],
            {
                "solves": 4,
                "iterations": 135,
                "eigenvalue_seconds": 1.0,
                "seconds": 1.5,
                "per_solve": {"iterations": [15, 60], "seconds": [0.25, 0.25]},
            },
        )
        self.assertEqual(c["runs"]["none"]["per_solve"]["iterations"], [900, 4000])

    def test_report_gives_the_ratios_and_the_range_per_solve(self):
        c = compare(
            ACCELERATION,
            core(),
            log([4000, 1200, 900, 3000], total=9.0),
            core(),
            log([60, 20, 15, 40], total=1.5),
        )
        summary = {"format": ACCELERATION.format, "cases": {"small": c}}
        text = "\n".join(report(summary))
        self.assertIn("外側の反復 9,100 → 135 回（0.015 倍）", text)
        self.assertIn("1 回の外側の反復 900〜4,000 → 15〜60 回", text)


class Parse(unittest.TestCase):
    def test_an_empty_path_in_the_list_is_rejected(self):
        with self.assertRaises(SystemExit):
            parse("case=a,")

    def test_a_case_names_two_inputs(self):
        self.assertEqual(
            parse("case=a.toml,b.toml"), ("case", [Path("a.toml"), Path("b.toml")])
        )
        with self.assertRaises(SystemExit):
            parse("case=a.msgpack,a.log,b.msgpack,b.log")


if __name__ == "__main__":
    unittest.main()
