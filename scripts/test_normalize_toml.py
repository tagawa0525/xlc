"""normalize_toml.py の検査（DVC の deps に置く TOML の値の正規形）。

実行: python3 -m unittest discover -s scripts -p 'test_normalize_toml.py'
"""

import functools
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import tomllib
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import normalize_toml as norm

ROOT = Path(__file__).resolve().parent.parent


def form(text: str) -> str:
    return norm.normalize(textwrap.dedent(text))


BASE = """
title = "pin"
order = ["U-235", "U-238", "O-16"]

[geometry]
pitch = 1.26
rings = 3
reflective = true

[[material]]
name = "fuel"
density = 10.257

[[material]]
name = "water"
density = 0.743
"""


class InvariantTest(unittest.TestCase):
    def test_comment_changes_keep_the_form(self):
        changed = (
            "# 先頭のコメント\n"
            + BASE.replace("rings = 3", "rings = 3  # 行末のコメント")
            + "# 末尾のコメント\n"
        )
        self.assertEqual(form(changed), form(BASE))

    def test_whitespace_and_line_breaks_keep_the_form(self):
        changed = BASE.replace("pitch = 1.26", "pitch   =   1.26").replace(
            'order = ["U-235", "U-238", "O-16"]',
            'order = [\n  "U-235",\n  "U-238",\n  "O-16",\n]\n\n',
        )
        self.assertEqual(form(changed), form(BASE))

    def test_quote_styles_keep_the_form(self):
        self.assertEqual(
            form('a = "x\\\\y"\nb = "\\u00e9"\n'), form("a = 'x\\y'\nb = \"é\"\n")
        )

    def test_number_spellings_of_the_same_value_keep_the_form(self):
        self.assertEqual(
            form("a = 1e10\nb = 0.10\nc = 1_000\nd = 0x10\n"),
            form("a = 1.0e10\nb = 0.1\nc = 1000\nd = 16\n"),
        )

    def test_table_spellings_of_the_same_value_keep_the_form(self):
        self.assertEqual(
            form("[a]\nb = 1\nc.d = 2\n"), form("a = { b = 1, c = { d = 2 } }\n")
        )


class ValueChangeTest(unittest.TestCase):
    def test_value_change_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE.replace("1.26", "1.27")))

    def test_integer_and_float_differ(self):
        self.assertNotEqual(form("a = 1\n"), form("a = 1.0\n"))

    def test_bool_and_integer_differ(self):
        self.assertNotEqual(form("a = true\n"), form("a = 1\n"))

    def test_signed_zero_differs(self):
        self.assertNotEqual(form("a = 0.0\n"), form("a = -0.0\n"))

    def test_key_order_changes_the_form(self):
        # 核種の並びの順で出力を書くプログラムがあるので、表のキーの順も値に数える（ADR-0212 の (b) の条件 3）
        self.assertNotEqual(form("a = 1\nb = 2\n"), form("b = 2\na = 1\n"))

    def test_array_order_changes_the_form(self):
        self.assertNotEqual(
            form(BASE), form(BASE.replace('"U-235", "U-238"', '"U-238", "U-235"'))
        )

    def test_string_and_number_differ(self):
        self.assertNotEqual(form('a = "1"\n'), form("a = 1\n"))

    def test_added_key_changes_the_form(self):
        self.assertNotEqual(form(BASE), form(BASE + "extra = 0\n"))

    def test_form_reads_back_to_the_parsed_values_in_order(self):
        values = tomllib.loads(BASE)
        back = json.loads(norm.normalize(BASE))
        self.assertEqual(back, values)
        self.assertEqual(list(back), list(values))
        self.assertEqual(list(back["geometry"]), list(values["geometry"]))


class FailureTest(unittest.TestCase):
    def test_syntax_error_raises(self):
        with self.assertRaises(tomllib.TOMLDecodeError):
            form("a = \n")

    def test_datetime_raises(self):
        with self.assertRaises(TypeError):
            form("a = 2026-11-19\n")

    def test_nan_and_infinity_raise(self):
        for text in ["a = nan\n", "a = inf\n", "a = -inf\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                form(text)


class MainTest(unittest.TestCase):
    def test_writes_the_form_creating_parent_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.toml")
            src.write_text("# コメント\na = 1\n", encoding="utf-8")
            dst = Path(tmp, "out", "deep", "a.toml.json")
            norm.main([str(src), str(dst)])
            self.assertEqual(dst.read_text(encoding="utf-8"), norm.normalize("a = 1\n"))


class InputMd5Test(unittest.TestCase):
    def test_toml_input_is_the_md5_of_its_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.toml")
            src.write_text("a = 1  # コメント\n", encoding="utf-8")
            self.assertEqual(
                norm.input_md5(src),
                hashlib.md5(norm.normalize("a = 1\n").encode()).hexdigest(),
            )

    def test_other_input_is_the_md5_of_its_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp, "a.msgpack")
            src.write_bytes(b"\x81\xa1a\x01")
            self.assertEqual(
                norm.input_md5(src), hashlib.md5(b"\x81\xa1a\x01").hexdigest()
            )


def tracked(pattern: str) -> list[str]:
    return subprocess.run(
        ["git", "ls-files", "--", pattern],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()


def toml_values_items() -> list[str]:
    """リポジトリの直下の dvc.yaml の toml-values の foreach の並び"""
    stages = yaml.safe_load((ROOT / "dvc.yaml").read_text(encoding="utf-8"))["stages"]
    if "toml-values" not in stages:
        raise AssertionError("dvc.yaml has no toml-values stage")
    return stages["toml-values"]["foreach"]


def load_vars(base: Path, entries: list) -> dict:
    """dvc.yaml の vars（ファイルの名前か、その場の辞書）を 1 つの辞書にする"""
    merged = {}
    for entry in entries:
        if isinstance(entry, dict):
            merged.update(entry)
        elif entry.endswith(".toml"):
            merged.update(tomllib.loads((base / entry).read_text(encoding="utf-8")))
        elif entry.endswith((".yaml", ".yml")):
            merged.update(yaml.safe_load((base / entry).read_text(encoding="utf-8")))
        else:
            raise AssertionError(f"{base}: unsupported vars entry {entry!r}")
    return merged


def lookup(context: dict, name: str):
    value = context
    for part in name.split("."):
        value = value[part]
    return value


def resolve(text: str, context: dict) -> str:
    """`${name}` と `${item.attr}` を埋める（DVC の補間のうち、例の dvc.yaml が使う形だけ。知らない名前は KeyError）"""
    return re.sub(r"\$\{([^}]+)\}", lambda m: str(lookup(context, m.group(1))), text)


def stage_paths(base: Path, spec: dict, field: str, context: dict) -> list[str]:
    entries = spec.get(field) or []
    # outs の項目は、`cache: false` などを付けると {パス: 設定} の辞書になる
    names = [next(iter(e)) if isinstance(e, dict) else e for e in entries]
    return [
        os.path.relpath(os.path.normpath(base / resolve(n, context)), ROOT)
        for n in names
    ]


@functools.cache
def example_stages() -> list[dict]:
    """例の dvc.yaml のステージごとの、解決した deps と outs（リポジトリの直下からのパス）

    DVC の Python の API は python-tests の devShell に無いので、foreach と補間を例が使う形だけ解く。foreach の値は
    並び、`${名前}`（vars の並びか辞書）。辞書なら `${key}` と `${item.属性}`、並びなら `${item}` を埋める。
    """
    result = []
    for name in tracked(":(glob)examples/**/dvc.yaml"):
        base = (ROOT / name).parent
        doc = yaml.safe_load((ROOT / name).read_text(encoding="utf-8"))
        variables = load_vars(base, doc.get("vars", []))
        for stage, body in (doc.get("stages") or {}).items():
            if "foreach" in body:
                items = body["foreach"]
                if isinstance(items, str):
                    items = lookup(
                        variables, items.removeprefix("${").removesuffix("}")
                    )
                pairs = (
                    [(key, {"key": key, "item": value}) for key, value in items.items()]
                    if isinstance(items, dict)
                    else [(item, {"item": item}) for item in items]
                )
                expanded = [
                    (f"{stage}@{key}", body["do"], extra) for key, extra in pairs
                ]
            else:
                expanded = [(stage, body, {})]
            for label, spec, extra in expanded:
                context = {**variables, **extra}

                result.append(
                    {
                        "stage": f"{name}:{label}",
                        "deps": stage_paths(base, spec, "deps", context),
                        "outs": stage_paths(base, spec, "outs", context),
                    }
                )
    return result


def form_path(source: str) -> str:
    return f"toml-values/{source}.json"


class RepositoryDepsTest(unittest.TestCase):
    def test_example_stages_depend_on_the_form_not_the_toml_file(self):
        # ステージの出力の TOML（xs-range.toml）は、そのステージのプログラムが書くバイトなので、そのまま deps にする
        outputs = {out for stage in example_stages() for out in stage["outs"]}
        bare = [
            f"{stage['stage']}: {dep}"
            for stage in example_stages()
            for dep in stage["deps"]
            if dep.endswith(".toml") and dep not in outputs
        ]
        self.assertEqual(bare, [])

    def test_every_form_a_stage_depends_on_is_written_by_toml_values(self):
        written = {form_path(item) for item in toml_values_items()}
        unknown = [
            f"{stage['stage']}: {dep}"
            for stage in example_stages()
            for dep in stage["deps"]
            if dep.startswith("toml-values/") and dep not in written
        ]
        self.assertEqual(unknown, [])

    def test_toml_values_writes_only_forms_some_stage_depends_on(self):
        used = {dep for stage in example_stages() for dep in stage["deps"]}
        unused = [item for item in toml_values_items() if form_path(item) not in used]
        self.assertEqual(unused, [])

    def test_toml_values_reads_tracked_inputs(self):
        files = set(tracked("*.toml"))
        outputs = {out for stage in example_stages() for out in stage["outs"]}
        wrong = [
            item for item in toml_values_items() if item not in files or item in outputs
        ]
        self.assertEqual(wrong, [])

    def test_every_form_matches_its_toml_file(self):
        # tomllib と json の書き方は Python の版で変わらないので（.py の正規形と違い）、すべてを比べ直す
        stale = [
            item
            for item in toml_values_items()
            if (ROOT / form_path(item)).read_text(encoding="utf-8")
            != norm.normalize((ROOT / item).read_text(encoding="utf-8"))
        ]
        self.assertEqual(stale, [])


class JsonInputsTest(unittest.TestCase):
    def test_tracked_json_inputs_match_the_locked_deps(self):
        """比較のステージの JSON の `inputs` が、そのステージの dvc.lock の deps（TOML は正規形）と同じ md5 を持つ（#0295）"""
        files = set(tracked("*.json"))
        problems = []
        for lock in tracked(":(glob)examples/**/dvc.lock"):
            base = (ROOT / lock).parent
            stages = yaml.safe_load((ROOT / lock).read_text(encoding="utf-8"))["stages"]
            for stage, body in stages.items():
                deps = {dep["path"]: dep.get("md5") for dep in body.get("deps") or []}
                for out in body.get("outs") or []:
                    path = os.path.relpath(os.path.normpath(base / out["path"]), ROOT)
                    if not path.endswith(".json") or path not in files:
                        continue
                    summary = json.loads((ROOT / path).read_text(encoding="utf-8"))
                    if not isinstance(summary, dict) or "inputs" not in summary:
                        continue
                    for item in summary["inputs"]:
                        dep = item["path"]
                        if dep.endswith(".toml"):
                            source = os.path.relpath(os.path.normpath(base / dep), ROOT)
                            dep = os.path.relpath(ROOT / form_path(source), base)
                        if deps.get(dep) != item["md5"]:
                            problems.append(
                                f"{path}: input {item['path']} has {item['md5']}, "
                                f"dvc.lock {lock}:{stage} has {deps.get(dep)} for {dep}"
                            )
        self.assertEqual(problems, [])


if __name__ == "__main__":
    unittest.main()
