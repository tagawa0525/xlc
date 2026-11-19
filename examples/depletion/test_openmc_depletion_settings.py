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


class _ResolveDepletion(ast.NodeTransformer):
    """`pin.YIELD_CUTOFF`（`import openmc_depletion as pin`）を `YIELD_CUTOFF` にそろえ、式の中の裸の名前を `bare` に集める。"""

    def __init__(self, aliases: set[str]):
        self.aliases = aliases
        self.bare: set[str] = set()

    def visit_Attribute(self, node: ast.Attribute) -> ast.AST:
        if isinstance(node.value, ast.Name) and node.value.id in self.aliases:
            return ast.Name(id=node.attr, ctx=node.ctx)
        return self.generic_visit(node)

    def visit_Name(self, node: ast.Name) -> ast.AST:
        self.bare.add(node.id)
        return node


def depletion_aliases(tree: ast.AST) -> set[str]:
    """`import openmc_depletion as <名前>` の名前（関数の中の import を含む）"""
    return {
        alias.asname or alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
        if alias.name == "openmc_depletion"
    }


def imported_from_depletion(tree: ast.AST) -> set[str]:
    """`from openmc_depletion import ...` で取った名前"""
    return {
        alias.asname or alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module == "openmc_depletion"
        for alias in node.names
    }


def module_assigned(tree: ast.Module) -> set[str]:
    names = set()
    for node in tree.body:
        targets = node.targets if isinstance(node, ast.Assign) else []
        if isinstance(node, ast.AnnAssign):
            targets = [node.target]
        names |= {t.id for t in targets if isinstance(t, ast.Name)}
    return names


def integrator_class(tree: ast.AST, callee: str) -> str:
    """積分器の呼び出しの先のクラス。`f()(...)` の形は、関数 f の中で定義して返すクラスの基底"""
    if not callee.endswith("()"):
        return callee
    (function,) = [
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef) and n.name == callee[:-2]
    ]
    (cls,) = [n for n in ast.walk(function) if isinstance(n, ast.ClassDef)]
    (base,) = cls.bases
    return ast.unparse(base)


def calls_with(tree: ast.AST, keywords: set[str]) -> list[dict[str, str]]:
    """`keywords` のどれかを引数に持つ呼び出しごとに、その引数の名前と式（ソースの形）、呼び出しの先（`callee`）、
    式の中の openmc_depletion から取っていない裸の名前（`local`。並べて文字列にする）"""
    aliases = depletion_aliases(tree)
    imported = imported_from_depletion(tree)
    assigned = module_assigned(tree) if isinstance(tree, ast.Module) else set()
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        args = {k.arg: k.value for k in node.keywords if k.arg in keywords}
        if not args:
            continue
        resolve = _ResolveDepletion(aliases)
        call = {name: ast.unparse(resolve.visit(value)) for name, value in args.items()}
        names = resolve.bare
        resolve.bare = set()
        callee = ast.unparse(resolve.visit(node.func))
        call["callee"] = integrator_class(tree, callee)
        # `f()(...)` の f は同じファイルの関数で、比べるのはそれが返すクラス（`callee`）
        if not callee.endswith("()"):
            names |= resolve.bare
        # 呼び出しの先の `openmc.deplete...` の `openmc` は比べる対象でない
        local = {n for n in names - {"openmc"} if n not in imported or n in assigned}
        call["local"] = ",".join(sorted(local))
        found.append(call)
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
                self.assertEqual(
                    set(operators[0]), OPERATOR_KEYWORDS | {"callee", "local"}
                )
                self.assertEqual(
                    set(integrators[0]), INTEGRATOR_KEYWORDS | {"callee", "local"}
                )

    def test_copies_match_the_pincell_solve(self):
        # 名前（YIELD_CUTOFF、TimedOperator）は openmc_depletion.py のものを使い、コピーの側で定義し直さない
        operators, integrators = settings(HERE / TOOLS[0])
        expected = (
            [{**operators[0], "local": ""}],
            [{**integrators[0], "local": ""}],
        )
        for name in TOOLS[1:]:
            with self.subTest(tool=name):
                self.assertEqual(settings(HERE / name), expected)

    def test_the_pincell_solve_uses_cecm_and_the_timed_operator(self):
        operators, integrators = settings(HERE / TOOLS[0])
        self.assertEqual(operators[0]["callee"], "TimedOperator")
        self.assertEqual(integrators[0]["callee"], "openmc.deplete.CECMIntegrator")

    def test_a_copy_that_defines_its_own_cutoff_is_caught(self):
        source = (HERE / TOOLS[3]).read_text() + "\nYIELD_CUTOFF = 1.0e3\n"
        operators = calls_with(ast.parse(source), OPERATOR_KEYWORDS)
        self.assertEqual(operators[0]["local"], "YIELD_CUTOFF")

    def test_a_copy_with_another_integrator_is_caught(self):
        source = (
            (HERE / TOOLS[2])
            .read_text()
            .replace(
                "openmc.deplete.CECMIntegrator(", "openmc.deplete.PredictorIntegrator("
            )
        )
        (integrator,) = calls_with(ast.parse(source), INTEGRATOR_KEYWORDS)
        self.assertEqual(integrator["callee"], "openmc.deplete.PredictorIntegrator")

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
