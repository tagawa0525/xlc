"""rewrite-lock-paths.py の検査（ADR-0260 の「確認」）。

使い捨ての git と DVC のリポジトリでステージを流してコミットし、ファイルを動かしてツールを実行し、`dvc status` を動かす前と
比べる。PATH の dvc（default の devShell の flake.nix のラッパー）を本物で使う。

実行: python3 -m unittest discover -s scripts -p 'test_rewrite_lock_paths.py'
"""

import importlib.util
import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path

import yaml

SCRIPTS = Path(__file__).resolve().parent
TOOL = SCRIPTS / "rewrite-lock-paths.py"
# 使い捨てのリポジトリに、利用者の git の設定（hook、署名）を持ち込まない
GIT_ENV = {
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_SYSTEM": os.devnull,
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
}


class Repo:
    def __init__(self, root: Path):
        self.root = root
        self.env = dict(os.environ, **GIT_ENV)
        self.env.pop("PYTHONPATH", None)
        self.run("git", "init", "-q", "-b", "t")
        self.run("dvc", "init", "-q")
        self.run("dvc", "config", "core.autostage", "false")

    def run(self, *args, check=True):
        result = subprocess.run(
            args,
            cwd=self.root,
            env=self.env,
            capture_output=True,
            text=True,
            check=False,
        )
        if check and result.returncode != 0:
            raise AssertionError(
                f"{' '.join(args)}: {result.returncode}\n{result.stdout}\n{result.stderr}"
            )
        return result

    def write(self, path, text):
        p = self.root / path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(textwrap.dedent(text))

    def read(self, path):
        return (self.root / path).read_text()

    def lock(self, path="dvc.lock"):
        return yaml.safe_load(self.read(path))

    def commit(self, message="chore: snapshot"):
        self.run("git", "add", "--all", ".")
        self.run("git", "commit", "-q", "--allow-empty", "-m", message)

    def repro(self, *targets):
        self.run("dvc", "repro", "-q", *targets)

    def status(self, *targets):
        return json.loads(self.run("dvc", "status", "--json", *targets).stdout or "{}")

    def mv(self, src, dst):
        (self.root / dst).parent.mkdir(parents=True, exist_ok=True)
        if (
            self.run("git", "ls-files", "--error-unmatch", src, check=False).returncode
            == 0
        ):
            self.run("git", "mv", src, dst)
        else:
            shutil.move(self.root / src, self.root / dst)

    def tool(self, *args):
        return self.run("python3", str(TOOL), *args, check=False)

    def remove_outputs(self, *paths):
        for p in paths:
            target = self.root / p
            if target.is_dir():
                shutil.rmtree(target)
            else:
                target.unlink()


def lines(result, prefix):
    return [
        line[len(prefix) :].strip()
        for line in result.stdout.splitlines()
        if line.startswith(prefix)
    ]


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Repo(Path(self.tmp.name))
        self.r = self.repo

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class MovingOnlyRewritesPaths(Base):
    def setUp(self):
        super().setUp()
        r = self.r
        r.write("a/in.txt", "same\n")
        r.write("a/twin.txt", "same\n")
        r.write("a/data/x.txt", "x\n")
        r.write("a/p.toml", "threads = 2\n[t]\nv = 1\n")
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              make:
                cmd: cat in.txt twin.txt > o.txt && mkdir -p od && cp data/x.txt od/ && cat p.toml > ptxt
                deps:
                  - in.txt
                  - twin.txt
                  - data
                params:
                  - p.toml:
                      - threads
                      - t
                outs:
                  - o.txt
                  - od
                  - ptxt
              each:
                foreach:
                  - in.txt
                do:
                  cmd: cat ${item} o.txt > ${item}.out
                  deps: ["${item}", o.txt]
                  outs: ["${item}.out"]
            """,
        )
        r.repro("a/dvc.yaml")
        r.commit()
        self.before = r.lock("a/dvc.lock")

    def test_moving_inputs_dirs_params_and_foreach_items_rewrites_only_paths(self):
        r = self.r
        r.mv("a/in.txt", "a/inputs/in.txt")
        r.mv("a/twin.txt", "a/inputs/twin.txt")
        r.mv("a/data", "a/inputs/data")
        r.mv("a/p.toml", "a/conf/p.toml")
        r.mv("a/od", "a/out/od")
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              make:
                cmd: cat inputs/in.txt inputs/twin.txt > o.txt && mkdir -p out/od && cp inputs/data/x.txt out/od/ && cat conf/p.toml > ptxt
                deps:
                  - inputs/in.txt
                  - inputs/twin.txt
                  - inputs/data
                params:
                  - conf/p.toml:
                      - threads
                      - t
                outs:
                  - o.txt
                  - out/od
                  - ptxt
              each:
                foreach:
                  - inputs/in.txt
                do:
                  cmd: cat ${item} o.txt > ${item}.out
                  deps: ["${item}", o.txt]
                  outs: ["${item}.out"]
            """,
        )
        r.mv("a/in.txt.out", "a/inputs/in.txt.out")
        self.assertNotEqual(r.status("a/dvc.yaml"), {})

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(r.status("a/dvc.yaml"), {})
        after = r.lock("a/dvc.lock")
        self.assertEqual(set(after["stages"]), {"make", "each@inputs/in.txt"})
        make = after["stages"]["make"]
        self.assertEqual(
            [d["path"] for d in make["deps"]],
            ["inputs/data", "inputs/in.txt", "inputs/twin.txt"],
        )
        self.assertEqual(list(make["params"]), ["conf/p.toml"])
        self.assertEqual(
            make["params"]["conf/p.toml"],
            self.before["stages"]["make"]["params"]["p.toml"],
        )
        old_md5 = sorted(d["md5"] for d in self.before["stages"]["make"]["deps"])
        self.assertEqual(sorted(d["md5"] for d in make["deps"]), old_md5)
        self.assertEqual(
            {o["path"]: o["md5"] for o in make["outs"]}["out/od"],
            {o["path"]: o["md5"] for o in self.before["stages"]["make"]["outs"]}["od"],
        )
        self.assertEqual(len(lines(result, "REWRITTEN:")), 2)

    def test_frozen_stage_is_rewritten_the_same_way(self):
        r = self.r
        text = r.read("a/dvc.yaml").replace("  make:\n", "  make:\n    frozen: true\n")
        r.write("a/dvc.yaml", text)
        r.commit()
        r.mv("a/twin.txt", "a/twin2.txt")
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("twin.txt", "twin2.txt"))
        before = r.status("a/dvc.yaml")

        result = r.tool()

        self.assert_ok(result)
        self.assertNotEqual(before, {})
        self.assertEqual(r.status("a/dvc.yaml"), {})
        deps = [d["path"] for d in r.lock("a/dvc.lock")["stages"]["make"]["deps"]]
        self.assertIn("twin2.txt", deps)

    def test_renaming_a_file_inside_a_dir_dep_is_not_rewritten(self):
        r = self.r
        r.mv("a/data", "a/data2")
        r.mv("a/data2/x.txt", "a/data2/y.txt")
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace("data/x.txt", "data2/y.txt")
            .replace("- data\n", "- data2\n"),
        )

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertEqual(
            [l.split(":")[1] for l in lines(result, "UNVERIFIED:")], ["make"]
        )
        self.assertEqual(
            r.lock("a/dvc.lock")["stages"]["make"], self.before["stages"]["make"]
        )

    def test_cmd_change_other_than_paths_is_not_rewritten(self):
        r = self.r
        r.mv("a/twin.txt", "a/twin2.txt")
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace("twin.txt", "twin2.txt")
            .replace("cat in.txt", "cat -n in.txt"),
        )

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertIn("a/dvc.yaml:make", result.stdout)
        self.assertEqual(
            r.lock("a/dvc.lock")["stages"]["make"], self.before["stages"]["make"]
        )

    def test_output_in_the_work_tree_with_another_md5_is_not_rewritten(self):
        r = self.r
        r.mv("a/o.txt", "a/o2.txt")
        r.write("a/o2.txt", "different\n")
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("o.txt", "o2.txt"))

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertIn("UNVERIFIED: a/dvc.yaml:make", result.stdout)

    def test_duplicating_a_stage_does_not_make_a_lock_entry_for_the_copy(self):
        r = self.r
        text = r.read("a/dvc.yaml")
        copy = (
            text.split("  each:")[0]
            .split("  make:", 1)[1]
            .replace("o.txt", "o3.txt")
            .replace(" od", " od3")
            .replace("ptxt", "ptxt3")
        )
        r.write(
            "a/dvc.yaml", text.replace("  each:", "  make-copy:" + copy + "  each:")
        )

        result = r.tool()

        self.assert_ok(result)
        self.assertNotIn("make-copy", r.lock("a/dvc.lock")["stages"])
        self.assertEqual(lines(result, "NO-ENTRY:"), ["a/dvc.yaml:make-copy"])


class UpstreamDirectory(Base):
    def test_cmd_naming_an_upstream_output_directory_follows_its_outs(self):
        r = self.r
        r.write("a/x.txt", "x\n")
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              up:
                cmd: mkdir -p od && cp x.txt od/x.txt
                deps: [x.txt]
                outs: [od]
              down:
                cmd: ls od > down.out && cat od/x.txt >> down.out
                deps: [od/x.txt]
                outs: [down.out]
            """,
        )
        r.repro("a/dvc.yaml")
        r.commit()
        r.remove_outputs("a/od", "a/down.out")
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("od", "res/od"))

        result = r.tool()

        self.assert_ok(result)
        r.run("dvc", "checkout", "-q")
        self.assertEqual(r.status("a/dvc.yaml"), {})


class ChangedBeforeMoving(Base):
    def setUp(self):
        super().setUp()
        r = self.r
        r.write("a/x.txt", "x\n")
        r.write("a/y.txt", "y\n")
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              up:
                cmd: cat x.txt y.txt > up.out
                deps: [x.txt, y.txt]
                outs: [up.out]
              down:
                cmd: cat up.out > down.out
                deps: [up.out]
                outs: [down.out]
            """,
        )
        r.repro("a/dvc.yaml")
        r.commit()
        # 動かす前から changed: x.txt を lock と違う中身にしてコミットする
        r.write("a/x.txt", "x2\n")
        r.commit()

    def test_stage_keeps_its_reason_and_outputs_return_to_the_new_paths(self):
        r = self.r
        before = r.status("a/dvc.yaml")
        self.assertEqual(list(before), ["a/dvc.yaml:up"])
        r.remove_outputs("a/up.out", "a/down.out")
        # lock と違う x.txt も、lock と同じ y.txt も、出力も動かす
        r.mv("a/x.txt", "a/in/x.txt")
        r.mv("a/y.txt", "a/in/y.txt")
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace("x.txt", "in/x.txt")
            .replace("y.txt", "in/y.txt")
            .replace("up.out", "res/up.out"),
        )

        result = r.tool()

        self.assert_ok(result)
        r.run("dvc", "checkout", "-q")
        self.assertTrue((r.root / "a/res/up.out").exists())
        self.assertEqual(
            r.status("a/dvc.yaml"),
            {"a/dvc.yaml:up": [{"changed deps": {"a/in/x.txt": "modified"}}]},
        )
        self.assertEqual(
            r.lock("a/dvc.lock")["stages"]["up"]["deps"][0]["md5"],
            yaml.safe_load(r.run("git", "show", "HEAD:a/dvc.lock").stdout)["stages"][
                "up"
            ]["deps"][0]["md5"],
        )

    def test_without_outputs_upstream_outs_and_downstream_deps_are_rewritten(self):
        r = self.r
        r.remove_outputs("a/up.out", "a/down.out")
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace("up.out", "res/up.out")
            .replace("down.out", "res/down.out"),
        )

        result = r.tool()

        self.assert_ok(result)
        stages = r.lock("a/dvc.lock")["stages"]
        self.assertEqual(stages["down"]["deps"][0]["path"], "res/up.out")
        self.assertEqual(stages["down"]["outs"][0]["path"], "res/down.out")
        r.run("dvc", "checkout", "-q")
        self.assertEqual(list(r.status("a/dvc.yaml")), ["a/dvc.yaml:up"])


class MovingDvcYaml(Base):
    def setUp(self):
        super().setUp()
        r = self.r
        r.write("shared.txt", "s\n")
        r.write("a/in.txt", "i\n")
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              one:
                cmd: cat ../shared.txt in.txt > one.out && echo ..#attr --out=one.out
                deps: [../shared.txt, in.txt]
                outs: [one.out]
              two:
                cmd: cat one.out ../shared.txt > two.out
                deps: [one.out, ../shared.txt]
                outs: [two.out]
            """,
        )
        r.repro("a/dvc.yaml")
        r.commit()

    def test_moving_the_dvc_yaml_deeper_rewrites_the_relative_paths(self):
        r = self.r
        for f in ("dvc.yaml", "dvc.lock", "in.txt", "one.out", "two.out"):
            r.mv(f"a/{f}", f"a/deep/{f}")
        r.write(
            "a/deep/dvc.yaml",
            r.read("a/deep/dvc.yaml")
            .replace("../shared.txt", "../../shared.txt")
            .replace("..#attr", "../..#attr"),
        )
        self.assertNotEqual(r.status("a/deep/dvc.yaml"), {})

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(r.status("a/deep/dvc.yaml"), {})
        deps = [d["path"] for d in r.lock("a/deep/dvc.lock")["stages"]["one"]["deps"]]
        self.assertEqual(deps, ["../../shared.txt", "in.txt"])

    def test_a_glob_in_the_cmd_of_a_moved_dvc_yaml_is_not_rewritten(self):
        r = self.r
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("cat one.out", "cat *.out"))
        r.repro("a/dvc.yaml")
        r.commit()
        for f in ("dvc.yaml", "dvc.lock", "in.txt", "one.out", "two.out"):
            r.mv(f"a/{f}", f"a/deep/{f}")
        r.write(
            "a/deep/dvc.yaml",
            r.read("a/deep/dvc.yaml")
            .replace("../shared.txt", "../../shared.txt")
            .replace("..#attr", "../..#attr"),
        )

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertEqual(
            [l.split(":")[1] for l in lines(result, "UNVERIFIED:")], ["two"]
        )
        self.assertIn("one", r.lock("a/deep/dvc.lock")["stages"])

    def test_a_stage_whose_definition_did_not_change_is_left_as_it_was(self):
        r = self.r
        # 動かす前から changed command（dvc.yaml の cmd を変えて流していない）
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml").replace(
                "cat one.out ../shared.txt", "cat ../shared.txt one.out"
            ),
        )
        r.commit()
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("in.txt", "in2.txt"))
        r.mv("a/in.txt", "a/in2.txt")

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(
            lines(result, "REWRITTEN:"), ["a/dvc.yaml:one（← a/dvc.lock:one）"]
        )
        self.assertEqual(list(r.status("a/dvc.yaml")), ["a/dvc.yaml:two"])

    def test_a_lock_changed_since_the_base_commit_stops_the_tool(self):
        r = self.r
        r.write("a/dvc.lock", r.read("a/dvc.lock").replace("one.out", "one2.out"))

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertIn("a/dvc.lock が動かす前のコミットと違う", result.stderr)

    def test_moving_a_stage_to_another_dvc_yaml(self):
        r = self.r
        text = r.read("a/dvc.yaml")
        one, two = text.split("  two:")
        r.write("a/dvc.yaml", one)
        r.write(
            "b/dvc.yaml",
            "stages:\n  two:"
            + two.replace("one.out", "../a/one.out")
            .replace("../shared.txt", "../shared.txt")
            .replace("> two.out", "> two.out"),
        )
        r.mv("a/two.out", "b/two.out")

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(r.status("a/dvc.yaml", "b/dvc.yaml"), {})
        self.assertEqual(list(r.lock("a/dvc.lock")["stages"]), ["one"])
        self.assertEqual(
            sorted(d["path"] for d in r.lock("b/dvc.lock")["stages"]["two"]["deps"]),
            ["../a/one.out", "../shared.txt"],
        )


NORMALIZE = SCRIPTS / "normalize_toml.py"


class TomlPaths(Base):
    """入力の TOML の中のパスを書き換えた移動（正規形の toml-values と、params）"""

    def setUp(self):
        super().setUp()
        r = self.r
        (r.root / "scripts").mkdir()
        shutil.copy(NORMALIZE, r.root / "scripts/normalize_toml.py")
        r.write("lib/lib.txt", "library\n")
        r.write(
            "a/case.toml", 'library = "../lib/lib.txt"\noutput = "o.txt"\nthreads = 1\n'
        )
        r.write(
            "dvc.yaml",
            """\
            stages:
              toml-values:
                foreach:
                  - a/case.toml
                do:
                  cmd: python3 scripts/normalize_toml.py ${item} toml-values/${item}.json
                  deps: [scripts/normalize_toml.py, "${item}"]
                  outs:
                    - toml-values/${item}.json:
                        cache: false
            """,
        )
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              run:
                cmd: cat case.toml ../lib/lib.txt > o.txt
                deps: [../toml-values/a/case.toml.json, ../lib/lib.txt]
                params:
                  - case.toml:
                      - library
                      - threads
                outs: [o.txt]
            """,
        )
        r.repro("dvc.yaml")
        r.repro("a/dvc.yaml")
        r.commit()

    def move_case_deeper(self, extra=""):
        r = self.r
        r.mv("a/case.toml", "a/deep/case.toml")
        r.write(
            "a/deep/case.toml",
            'library = "../../lib/lib.txt"\noutput = "o.txt"\nthreads = 1\n' + extra,
        )
        r.mv("toml-values/a/case.toml.json", "toml-values/a/deep/case.toml.json")
        r.run(
            "python3",
            "scripts/normalize_toml.py",
            "a/deep/case.toml",
            "toml-values/a/deep/case.toml.json",
        )
        r.write(
            "dvc.yaml", r.read("dvc.yaml").replace("a/case.toml", "a/deep/case.toml")
        )
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace(
                "../toml-values/a/case.toml.json",
                "../toml-values/a/deep/case.toml.json",
            )
            .replace("cat case.toml", "cat deep/case.toml")
            .replace("- case.toml:", "- deep/case.toml:"),
        )

    def test_toml_whose_relative_paths_were_rewritten_gets_new_md5_and_values(self):
        r = self.r
        self.move_case_deeper()

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(r.status("dvc.yaml", "a/dvc.yaml"), {})
        params = r.lock("a/dvc.lock")["stages"]["run"]["params"]
        self.assertEqual(
            params, {"deep/case.toml": {"library": "../../lib/lib.txt", "threads": 1}}
        )
        self.assertIn("toml-values@a/deep/case.toml", r.lock()["stages"])

    def test_toml_with_another_value_changed_is_not_rewritten(self):
        r = self.r
        self.move_case_deeper()
        r.write(
            "a/deep/case.toml",
            r.read("a/deep/case.toml").replace("threads = 1", "threads = 2"),
        )
        r.run(
            "python3",
            "scripts/normalize_toml.py",
            "a/deep/case.toml",
            "toml-values/a/deep/case.toml.json",
        )

        result = r.tool()

        self.assertEqual(result.returncode, 1)
        self.assertIn("UNVERIFIED: a/dvc.yaml:run", result.stdout)

    def test_toml_naming_an_output_follows_the_moved_outs(self):
        r = self.r
        r.remove_outputs("a/o.txt")
        r.write(
            "a/case.toml",
            r.read("a/case.toml").replace('output = "o.txt"', 'output = "res/o.txt"'),
        )
        r.run(
            "python3",
            "scripts/normalize_toml.py",
            "a/case.toml",
            "toml-values/a/case.toml.json",
        )
        r.write(
            "a/dvc.yaml",
            r.read("a/dvc.yaml")
            .replace("> o.txt", "> res/o.txt")
            .replace("outs: [o.txt]", "outs: [res/o.txt]"),
        )

        result = r.tool()

        self.assert_ok(result)
        r.run("dvc", "checkout", "-q")
        self.assertEqual(r.status("dvc.yaml", "a/dvc.yaml"), {})


class ComparisonStage(Base):
    def setUp(self):
        super().setUp()
        r = self.r
        r.write("a/in.txt", "i\n")
        r.write(
            "a/compare.py",
            """\
            import hashlib, json, sys
            data = open(sys.argv[1], "rb").read()
            json.dump({"inputs": [{"path": sys.argv[1], "md5": hashlib.md5(data).hexdigest()}]}, open(sys.argv[2], "w"))
            """,
        )
        r.write(
            "a/dvc.yaml",
            """\
            stages:
              compare:
                cmd: python3 compare.py in.txt summary.json
                deps: [compare.py, in.txt]
                outs: [summary.json]
            """,
        )
        r.repro("a/dvc.yaml")
        r.commit()

    def test_comparison_stage_with_a_moved_input_is_listed_apart_and_not_an_error(self):
        r = self.r
        r.mv("a/in.txt", "a/inputs/in.txt")
        r.write("a/dvc.yaml", r.read("a/dvc.yaml").replace("in.txt", "inputs/in.txt"))
        before = r.lock("a/dvc.lock")

        result = r.tool()

        self.assert_ok(result)
        self.assertEqual(
            [l.split(":")[1] for l in lines(result, "COMPARISON:")], ["compare"]
        )
        self.assertEqual(r.lock("a/dvc.lock"), before)


SPEC = importlib.util.spec_from_file_location("rewrite_lock_paths", TOOL)
rlp = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rlp)


class Words(unittest.TestCase):
    def test_cmd_is_split_into_shell_words_keeping_flake_references(self):
        self.assertEqual(
            rlp.cmd_parts("nix develop ../..#references -c python3 a.py 'x y'>o 2>err"),
            [
                [
                    "nix",
                    "develop",
                    "../..#references",
                    "-c",
                    "python3",
                    "a.py",
                    "x y",
                    ">",
                    "o",
                    "2",
                    ">",
                    "err",
                ]
            ],
        )

    def test_word_is_split_into_option_path_and_fragment(self):
        self.assertEqual(
            rlp.split_word("../..#references"), ("", "../..", "#references")
        )
        self.assertEqual(
            rlp.split_word("--json=out/a.json"), ("--json=", "out/a.json", "")
        )
        self.assertEqual(rlp.split_word("a.toml"), ("", "a.toml", ""))


class TomlValues(unittest.TestCase):
    def same(self, a, b):
        return (a, b) == ("../x", "../../x")

    def test_only_strings_naming_the_same_file_may_differ(self):
        old = {"lib": "../x", "n": 1, "t": {"v": [1, "a"]}}
        self.assertTrue(
            rlp.toml_values_differ_only_by_paths(
                old, {"lib": "../../x", "n": 1, "t": {"v": [1, "a"]}}, self.same
            )
        )
        self.assertFalse(
            rlp.toml_values_differ_only_by_paths(
                old, {"lib": "../../x", "n": 2, "t": {"v": [1, "a"]}}, self.same
            )
        )
        self.assertFalse(
            rlp.toml_values_differ_only_by_paths(
                old, {"lib": "../y", "n": 1, "t": {"v": [1, "a"]}}, self.same
            )
        )

    def test_key_order_and_types_must_match(self):
        self.assertFalse(
            rlp.toml_values_differ_only_by_paths(
                {"a": 1, "b": 2}, {"b": 2, "a": 1}, self.same
            )
        )
        self.assertFalse(
            rlp.toml_values_differ_only_by_paths({"a": 1}, {"a": 1.0}, self.same)
        )
        self.assertFalse(
            rlp.toml_values_differ_only_by_paths({"a": True}, {"a": 1}, self.same)
        )


class Matching(unittest.TestCase):
    def test_content_matches_are_preferred_when_they_cover_everything(self):
        # 0 は右の 0 と中身が同じ（順位 0）、右の 1 とは動かす前の中身が同じ（順位 1）
        self.assertEqual(
            rlp.match_by_content(2, [[(1, 1), (0, 0)], [(0, 1)]]), {0: 0, 1: 1}
        )

    def test_falls_back_to_lower_ranks_and_fails_without_a_full_matching(self):
        self.assertEqual(
            rlp.match_by_content(2, [[(0, 0)], [(1, 0), (2, 1)]]), {0: 0, 1: 1}
        )
        self.assertIsNone(rlp.match_by_content(2, [[(0, 0)], [(0, 0)]]))

    def test_inputs_are_read_from_lists_and_dicts(self):
        self.assertEqual(
            rlp.inputs_paths({"inputs": [{"path": "a", "md5": "x"}]}), ["a"]
        )
        self.assertEqual(rlp.inputs_paths({"inputs": {"b": "y"}}), ["b"])
        self.assertIsNone(rlp.inputs_paths({"k": 1}))


if __name__ == "__main__":
    unittest.main()
