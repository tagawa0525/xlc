#!/usr/bin/env python3
"""Python と TOML の deps の正規形を、ステージの deps に合わせて書き直す（ADR-0258 の「正規形の置き場」）。pre-commit の hook
`normalized-forms` の入口で、scripts/dvc-relock.py も `sync` を呼ぶ。

使い方（default の devShell で、リポジトリの直下から）: python3 scripts/normalized_forms.py

1. 直下の dvc.yaml の python-ast と toml-values の foreach を、追跡しているすべての dvc.yaml のステージが deps にする正規形の
   元のファイルの並び（並べたもの）にする。正規形を書くステージ自身の deps は数えない。どのステージも deps にしない正規形の
   ステージ（並びが空になるもの）は、何も書かずに止まる
2. 並びから外れた元のファイルの、直下の dvc.lock の項と、正規形のファイルを消す（dvc repro は消えた foreach の項の lock を残す）
3. 2 組を 1 つの `dvc repro -s` で実行する。dvc が失敗すれば止まる
4. 1〜3 が書き換えたファイル（直下の dvc.yaml、dvc.lock、python-ast/、toml-values/）を `REWRITTEN:` に出す。あれば終了コード 1
   （ほかの書式を直す hook と同じく、コミットを止める）、無ければ 0

例のステージの deps は、YAML を読み、foreach と補間を例が使う形だけ解いて得る（DVC の Python の API は、この hook と
テストを実行する Python に無い）。
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import tomllib
import yaml

# 正規形のステージ → 正規形のパスの前と後ろ（`python-ast/<元のパス>.ast`）。`dvc repro -s` は目標を挙げた順に実行し、DAG で
# 並べ直さない（dvc 3.67.1 の dvc/repo/reproduce.py）。toml-values は python-ast の normalize_toml.py.ast を読むので後に置く
FORMS = {
    "python-ast": ("python-ast/", ".ast"),
    "toml-values": ("toml-values/", ".json"),
}
REPRO = ["dvc", "repro", "-s", *(f"dvc.yaml:{stage}" for stage in FORMS)]


class Stop(Exception):
    pass


def tracked(root: Path, *patterns: str) -> list[str]:
    return subprocess.run(
        ["git", "ls-files", "--", *patterns],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()


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
            raise Stop(f"{base}: unsupported vars entry {entry!r}")
    return merged


def lookup(context: dict, name: str):
    value = context
    for part in name.split("."):
        value = value[part]
    return value


def resolve(text: str, context: dict) -> str:
    """`${name}` と `${item.attr}` を埋める（DVC の補間のうち、例の dvc.yaml が使う形だけ。知らない名前は KeyError）"""
    return re.sub(r"\$\{([^}]+)\}", lambda m: str(lookup(context, m.group(1))), text)


def stage_paths(
    root: Path, base: Path, spec: dict, field: str, context: dict
) -> list[str]:
    entries = spec.get(field) or []
    # outs の項目は、`cache: false` などを付けると {パス: 設定} の辞書になる
    names = [next(iter(e)) if isinstance(e, dict) else e for e in entries]
    return [
        os.path.relpath(os.path.normpath(base / resolve(n, context)), root)
        for n in names
    ]


def stages(root: Path, names: list[str]) -> list[dict]:
    """dvc.yaml ごとのステージの、解決した cmd、deps と outs（deps と outs はリポジトリの直下からのパス）

    foreach の値は並び、`${名前}`（vars の並びか辞書）。辞書なら `${key}` と `${item.属性}`、並びなら `${item}` を埋める。
    解けない dvc.yaml（知らない名前、do の無い foreach、壊れた YAML か TOML）は、その名前を挙げて Stop にする。
    """
    result = []
    for name in names:
        try:
            result += file_stages(root, name)
        except (
            KeyError,
            TypeError,
            AttributeError,
            yaml.YAMLError,
            tomllib.TOMLDecodeError,
        ) as e:
            raise Stop(
                f"{name}: cannot resolve the stages ({type(e).__name__}: {e})"
            ) from e
    return result


def file_stages(root: Path, name: str) -> list[dict]:
    result = []
    base = (root / name).parent
    doc = yaml.safe_load((root / name).read_text(encoding="utf-8"))
    variables = load_vars(base, doc.get("vars", []))
    for stage, body in (doc.get("stages") or {}).items():
        if "foreach" in body:
            items = body["foreach"]
            if isinstance(items, str):
                items = lookup(variables, items.removeprefix("${").removesuffix("}"))
            pairs = (
                [(key, {"key": key, "item": value}) for key, value in items.items()]
                if isinstance(items, dict)
                else [(item, {"item": item}) for item in items]
            )
            expanded = [(f"{stage}@{key}", body["do"], extra) for key, extra in pairs]
        else:
            expanded = [(stage, body, {})]
        for label, spec, extra in expanded:
            context = {**variables, **extra}
            cmd = spec.get("cmd", "")
            result.append(
                {
                    "stage": f"{name}:{label}",
                    "cmd": resolve(
                        " ".join(cmd) if isinstance(cmd, list) else cmd, context
                    ),
                    "deps": stage_paths(root, base, spec, "deps", context),
                    "outs": stage_paths(root, base, spec, "outs", context),
                }
            )
    return result


def form_items(root: Path) -> dict[str, list[str]]:
    """正規形のステージ → ステージが deps にする正規形の、元のファイルの並び（並べたもの）"""
    found = {stage: set() for stage in FORMS}
    for stage in stages(root, tracked(root, ":(glob)**/dvc.yaml")):
        for form, (prefix, suffix) in FORMS.items():
            # 正規形を書くステージは、自分の正規形を deps にできない（書くには先に自分が要る）
            if stage["stage"].startswith(f"dvc.yaml:{form}@"):
                continue
            found[form].update(
                dep.removeprefix(prefix).removesuffix(suffix)
                for dep in stage["deps"]
                if dep.startswith(prefix) and dep.endswith(suffix)
            )
    return {form: sorted(items) for form, items in found.items()}


def rewrite_foreach(text: str, form: str, items: list[str]) -> str:
    m = re.search(
        rf"^  {re.escape(form)}:\n    foreach:\n((?:      - \S+\n)+)    do:",
        text,
        re.MULTILINE,
    )
    if m is None:
        raise Stop(f"dvc.yaml has no {form} stage with a foreach list")
    listed = "".join(f"      - {item}\n" for item in items)
    return text[: m.start(1)] + listed + text[m.end(1) :]


def drop_lock_entries(text: str, form: str, keep: list[str]) -> str:
    """直下の dvc.lock から、keep に無い `<form>@<元のファイル>` の項を消す（ほかの行はバイトで残す）"""
    kept = set(keep)
    entry = re.compile(rf"^  {re.escape(form)}@(\S+):\n(?:    .*\n)*", re.MULTILINE)
    return entry.sub(lambda m: m.group(0) if m.group(1) in kept else "", text)


def sync(root: Path) -> None:
    """1 と 2: foreach を deps に合わせ、外れた項の lock と正規形を消す。並びが空になるなら、何も書かずに Stop"""
    items = form_items(root)
    empty = [form for form, listed in items.items() if not listed]
    if empty:
        raise Stop(
            f"no stage depends on a form written by {', '.join(empty)}; remove the stage from dvc.yaml "
            "or add the forms to the deps"
        )
    yaml_path, lock_path = root / "dvc.yaml", root / "dvc.lock"
    text = yaml_path.read_text(encoding="utf-8")
    lock = lock_path.read_text(encoding="utf-8")
    new_text, new_lock = text, lock
    for form, listed in items.items():
        new_text = rewrite_foreach(new_text, form, listed)
        new_lock = drop_lock_entries(new_lock, form, listed)
    if new_text != text:
        yaml_path.write_text(new_text, encoding="utf-8")
    if new_lock != lock:
        lock_path.write_text(new_lock, encoding="utf-8")
    for form, (prefix, suffix) in FORMS.items():
        keep = {f"{prefix}{item}{suffix}" for item in items[form]}
        for path in sorted((root / prefix).rglob(f"*{suffix}")):
            if path.relative_to(root).as_posix() not in keep:
                path.unlink()


def snapshot(root: Path) -> dict[str, bytes]:
    """hook が書き換えうるファイル → 中身"""
    paths = [root / "dvc.yaml", root / "dvc.lock"]
    for prefix, _ in FORMS.values():
        paths += (p for p in (root / prefix).rglob("*") if p.is_file())
    return {
        p.relative_to(root).as_posix(): p.read_bytes() for p in paths if p.is_file()
    }


def rewritten(before: dict[str, bytes], after: dict[str, bytes]) -> list[str]:
    return sorted(
        p for p in before.keys() | after.keys() if before.get(p) != after.get(p)
    )


def write_forms(root: Path) -> None:
    sync(root)
    # dvc は飛ばしたステージも 1 行ずつ出す（数百行）ので、失敗したときだけ出す
    result = subprocess.run(
        REPRO, cwd=root, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise Stop(f"{' '.join(REPRO)} failed with exit code {result.returncode}")


def main() -> int:
    root = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    )
    # dvc の目標（dvc.yaml:<ステージ>）は cwd からの相対なので、直下でだけ動かす
    if Path.cwd().resolve() != root.resolve():
        print(
            f"STOP: run from the root of the repository ({root}); the dvc targets are paths from there",
            file=sys.stderr,
        )
        return 1
    before = snapshot(root)
    try:
        write_forms(root)
        stop = None
    # dvc を起動できない（devShell の外）ときの OSError と、git ls-files の失敗も、止まる理由として出す
    except (Stop, OSError, subprocess.CalledProcessError) as e:
        stop = e
    changed = rewritten(before, snapshot(root))
    for path in changed:
        print(f"REWRITTEN: {path}")
    if stop is not None:
        print(f"STOP: {stop}", file=sys.stderr)
        return 1
    if changed:
        print(
            "NEXT: the Python and TOML forms were stale; git add the REWRITTEN: files and commit again",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
