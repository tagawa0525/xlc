"""OpenMC の燃焼の 4 つのツールが、オペレーターと積分器に同じ設定を渡すこと（#0132）。

openmc_depletion.py の `solve` の設定（`TimedOperator` と `CECMIntegrator` の引数）は、openmc_depletion_resume.py、
openmc_depletion_assembly.py、openmc_depletion_mox.py にコピーしてある（共通の関数にすると openmc-1C の重いステージが古くなる）。
片方だけを変えると、比べる参照解どうしが違う設定で流れて比較が黙って崩れるので、ソースの呼び出しの引数を突き合わせる。
OpenMC を import せずに ast で読む。

実行方法: python3 -m unittest discover -s examples/depletion -p 'test_openmc_depletion_settings.py'
"""

import ast
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = [
    "openmc_depletion.py",
    "openmc_depletion_resume.py",
    "openmc_depletion_assembly.py",
    "openmc_depletion_mox.py",
]
# 問題ごとに違ってよい引数（モデル、チェーン、再開の前の結果、ステップ、出力）は比べない
OPERATOR_KEYWORDS = {
    "normalization_mode",
    "fission_yield_mode",
    "fission_yield_opts",
    "reaction_rate_mode",
}
INTEGRATOR_KEYWORDS = {"timestep_units", "solver"}


class _StripModule(ast.NodeTransformer):
    """`pin.YIELD_CUTOFF` のような、import したモジュール越しの名前を `YIELD_CUTOFF` にそろえる。"""

    def visit_Attribute(self, node: ast.Attribute) -> ast.AST:
        if isinstance(node.value, ast.Name) and node.attr.isupper():
            return ast.Name(id=node.attr, ctx=node.ctx)
        return self.generic_visit(node)


def calls_with(tree: ast.AST, keywords: set[str]) -> list[dict[str, str]]:
    """`keywords` のどれかを引数に持つ呼び出しごとに、その引数の名前と式（ソースの形）。"""
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        args = {k.arg: k.value for k in node.keywords if k.arg in keywords}
        if args:
            found.append(
                {
                    name: ast.unparse(_StripModule().visit(value))
                    for name, value in args.items()
                }
            )
    return found


def settings(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    tree = ast.parse(path.read_text(), filename=str(path))
    return calls_with(tree, OPERATOR_KEYWORDS), calls_with(tree, INTEGRATOR_KEYWORDS)


class SameSettings(unittest.TestCase):
    def test_each_tool_passes_the_settings_once(self):
        for name in TOOLS:
            with self.subTest(tool=name):
                operators, integrators = settings(HERE / name)
                self.assertEqual(len(operators), 1, operators)
                self.assertEqual(len(integrators), 1, integrators)
                self.assertEqual(set(operators[0]), OPERATOR_KEYWORDS)
                self.assertEqual(set(integrators[0]), INTEGRATOR_KEYWORDS)

    def test_copies_match_the_pincell_solve(self):
        expected = settings(HERE / TOOLS[0])
        for name in TOOLS[1:]:
            with self.subTest(tool=name):
                self.assertEqual(settings(HERE / name), expected)

    def test_a_changed_copy_is_caught(self):
        source = (
            (HERE / TOOLS[0]).read_text().replace('solver="cram48"', 'solver="cram16"')
        )
        self.assertNotEqual(
            calls_with(ast.parse(source), INTEGRATOR_KEYWORDS),
            settings(HERE / TOOLS[0])[1],
        )


if __name__ == "__main__":
    unittest.main()
