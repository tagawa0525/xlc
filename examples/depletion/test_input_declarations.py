"""リポジトリの xs の入力の全部が、評価の欠陥を読み替える宣言（ADR-0272、ADR-0267、ADR-0274）を、版の走査の宣言と同じに持つこと。

宣言の要る評価の一覧は、版ごとの走査の入力 `process-survey/<版>/xs.toml`（`process_survey.py input` が `declarations.toml` から書く）。
追跡している TOML のうち `endf-directory` と `[nuclide]` を持つもの（xs の入力）の核種ごとに、同じ版の同じ評価のファイルの宣言の
項目を持つことを確かめる。放射化の核種は遮蔽しないので、非分離共鳴の宣言（`unresolved-negative-background`）は持たない（xs は書くと
止まる）。宣言が無いと xs はその評価で止まるので、本番と検証の入力を止まったまま残さないための検査。

実行: python3 -m unittest discover -s examples/depletion -p 'test_input_declarations.py'
"""

import re
import subprocess
import sys
import unittest
from pathlib import Path

import tomllib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
from process_survey import DECLARATION_KEYS

VERSIONS = {"VII.0", "VII.1", "VIII.0", "VIII.1"}
DIRECTORY = re.compile(r"ENDF-B-(VII\.0|VII\.1|VIII\.0|VIII\.1)_neutrons")
# 放射化の核種が取らない宣言
SHIELDING_ONLY = ("unresolved-negative-background",)
# 宣言を持たないことに理由のある入力: 非分離共鳴の読み方の比較（#0319）は、宣言の前の xs と、採らなかった読み方の xs で流す入力
EXEMPT = {
    "examples/depletion/issue-0319/VII.0-add.toml",
    "examples/depletion/issue-0319/VII.0-factors.toml",
    "examples/depletion/issue-0319/VIII.0-add.toml",
    "examples/depletion/issue-0319/VIII.0-factors.toml",
}


def required(version: str) -> dict[str, dict]:
    """版の評価のファイル → 宣言の項目"""
    entries = tomllib.loads(
        (HERE / "process-survey" / version / "xs.toml").read_text()
    )["nuclide"]
    return {
        e["endf"]: {k: e[k] for k in DECLARATION_KEYS if k in e}
        for e in entries.values()
        if any(k in e for k in DECLARATION_KEYS)
    }


def xs_inputs() -> list[tuple[str, dict]]:
    files = subprocess.run(
        ["git", "ls-files", "*.toml"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split()
    out = []
    for f in files:
        data = tomllib.loads((ROOT / f).read_text())
        if "endf-directory" in data and isinstance(data.get("nuclide"), dict):
            out.append((f, data))
    return out


def missing(data: dict, needed: dict[str, dict]) -> list[str]:
    """入力 `data` の核種のうち、宣言が `needed` と違うもの（核種: 足りない項目か違う項目）"""
    out = []
    for name, entry in data["nuclide"].items():
        items = needed.get(Path(entry["endf"]).name, {})
        if entry.get("activation"):
            items = {k: v for k, v in items.items() if k not in SHIELDING_ONLY}
        wrong = sorted(k for k, v in items.items() if entry.get(k) != v)
        if wrong:
            out.append(f"{name}: {wrong}")
    return out


class InputDeclarationsTest(unittest.TestCase):
    def test_every_xs_input_carries_the_declarations_of_its_evaluations(self):
        needed = {v: required(v) for v in VERSIONS}
        problems = []
        for f, data in xs_inputs():
            if f in EXEMPT or f.startswith("examples/depletion/process-survey/"):
                continue
            m = DIRECTORY.search(data["endf-directory"])
            if m is None:
                continue
            found = missing(data, needed[m[1]])
            if found:
                problems.append(f"{f}: {', '.join(found)}")
        self.assertEqual(problems, [])

    def test_an_activation_nuclide_needs_no_unresolved_declaration(self):
        needed = {
            "x.endf": {
                "unresolved-negative-background": "scale",
                "corrections": [{"line": 1}],
            }
        }
        data = {
            "nuclide": {
                "X": {
                    "endf": "x.endf",
                    "activation": True,
                    "corrections": [{"line": 1}],
                }
            }
        }
        self.assertEqual(missing(data, needed), [])
        data["nuclide"]["X"]["activation"] = False
        self.assertEqual(
            missing(data, needed), ["X: ['unresolved-negative-background']"]
        )


if __name__ == "__main__":
    unittest.main()
