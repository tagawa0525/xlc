#!/usr/bin/env python3
"""ファイルやディレクトリを動かした後に、同じファイルかを中身の md5 で判定して dvc.lock のパスだけを書き換える（ADR-0260）。

使い方（default の devShell で、リポジトリの直下から）:
    scripts/rewrite-lock-paths.py [--base <動かす前のコミット>]

`--base` の既定は HEAD（動かした変更をコミットする前に実行する）。移動をコミットした後に実行するときは、動かす前の
コミットを渡す。作業ツリーの dvc.lock は、どれも `--base` のどれかの dvc.lock と同じバイトでなければならない（動かした
lock は同じ中身で見つける）。書き換えた lock は作業ツリーに書き、コミットはしない。

- 判定の決まりは ADR-0260 の「決定と理由」にある。ステージごとに、deps、params、outs、cmd のすべてが確かめられたときだけ、
  lock の項の `path`、`params` のファイル、`cmd`、ステージの名前を書き換える。md5、size、nfiles は持ち越す（入力の TOML の
  中のパスだけを書き換えたときの、その TOML と正規形の項を除く）
- 出力の行: `REWRITTEN:`（書き換えたステージ）、`UNVERIFIED:`（確かめられず書き換えなかったステージと理由）、
  `COMPARISON:`（確かめられたが、出力の JSON の `inputs` が動かした入力を持つので書き換えなかった比較のステージ）、
  `NO-ENTRY:`（持ち越す lock の項が無いステージ）、`LEFT:`（どのステージにも対応しなかった、動かす前の lock の項）、
  `DVC-FILE:`（.dvc のファイルの確認で通らなかったもの）
- 終了コードは、`UNVERIFIED:` か `DVC-FILE:` があれば 1、ほかは 0。`COMPARISON:` は 1 の理由に数えない

DVC のライブラリ（ステージの展開、ファイルとディレクトリの md5、params の読み方）を使うので、PATH の dvc（flake.nix の
ラッパー）が使う Python で自分を実行し直す。
"""

import argparse
import copy
import importlib.util
import json
import os
import posixpath
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import tomllib

IN_DVC_ENV = "XLC_REWRITE_LOCK_PATHS_IN_DVC"
# .dvc-wrapped（nixpkgs の Python のラッパー）を __main__ でなく読むと、DVC の依存の site のディレクトリを足すだけで、
# dvc の CLI は動かない。その後でこの script を __main__ として実行する
BOOT = (
    "import runpy, sys; runpy.run_path(sys.argv[1], run_name='dvc_env'); "
    "sys.argv = sys.argv[2:]; runpy.run_path(sys.argv[0], run_name='__main__')"
)
LOCK = "dvc.lock"
# 計算のステージが deps に書く正規形と、cmd が読む元のファイル（ADR-0210、ADR-0254、ADR-0255）
NORMAL_FORMS = (
    ("python-ast/", ".ast", ""),
    ("toml-values/", ".json", ""),
    ("binary-forms/", ".md5", "bin/"),
)
TOML_FORM = ("toml-values/", ".json")
OPTION_WORD = re.compile(r"(-{1,2}[A-Za-z0-9_][A-Za-z0-9_-]*=)(.*)")
GLOB_CHARS = set("*?[")


class ToolError(Exception):
    pass


class Unverified(Exception):
    """ステージを確かめられない理由"""


def reexec_under_dvc(argv):
    if os.environ.get(IN_DVC_ENV):
        raise ToolError("DVC の Python でも dvc を import できない")
    dvc = shutil.which("dvc")
    if dvc is None:
        raise ToolError("dvc が PATH に無い（default の devShell で実行する）")
    m = re.search(
        r"^exec (\S+/bin/\.dvc-wrapped) ",
        Path(dvc).read_text(errors="replace"),
        re.MULTILINE,
    )
    if m is None:
        raise ToolError(
            f"{dvc} が flake.nix の dvc のラッパー（.dvc-wrapped を exec する）でない"
        )
    wrapped = m[1]
    with open(wrapped, encoding="utf-8") as f:
        first = f.readline()
    if not first.startswith("#!"):
        raise ToolError(f"{wrapped} の先頭に Python の #! が無い")
    interp = shlex.split(first[2:])
    env = dict(os.environ, PYTHONNOUSERSITE="true", **{IN_DVC_ENV: "1"})
    os.execve(
        interp[0],
        [*interp, "-c", BOOT, wrapped, str(Path(__file__).resolve()), *argv[1:]],
        env,
    )


# --- パス -----------------------------------------------------------------------------------------------------------


def join(dir_, path):
    """`dir_`（リポジトリの直下からのディレクトリ）から見た `path` を、リポジトリの直下からのパスにする"""
    if posixpath.isabs(path):
        return posixpath.normpath(path)
    joined = posixpath.normpath(posixpath.join(dir_, path))
    return "" if joined == "." else joined


def rel(path, dir_):
    return posixpath.relpath(path or ".", dir_ or ".")


def dirname(path):
    return posixpath.dirname(path)


def normal_form_source(path):
    """正規形のパスなら、cmd が読む元のファイルのパス。正規形でなければ None"""
    for prefix, suffix, source_prefix in NORMAL_FORMS:
        if path.startswith(prefix) and path.endswith(suffix):
            return source_prefix + path[len(prefix) : -len(suffix)]
    return None


def toml_source(path):
    """入力の TOML か、その正規形なら、元の TOML のパス。ほかは None"""
    prefix, suffix = TOML_FORM
    if path.startswith(prefix) and path.endswith(".toml" + suffix):
        return path[len(prefix) : -len(suffix)]
    if path.endswith(".toml"):
        return path
    return None


def cmd_parts(cmd):
    """cmd（文字列か文字列のリスト）を、要素ごとの shell の語のリストにする"""
    out = []
    for part in cmd if isinstance(cmd, list) else [cmd]:
        lex = shlex.shlex(part, posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        # bash は語の頭の # だけを注釈にする。`../..#references` の # を注釈にしない
        lex.commenters = ""
        out.append(list(lex))
    return out


def split_word(word):
    """語を (前, パスの部分, 後) に分ける。`<パス>#<属性>`（flake の参照）と `--<名前>=<パス>` を読む"""
    if "#" in word:
        path, frag = word.split("#", 1)
        return "", path, "#" + frag
    if m := OPTION_WORD.fullmatch(word):
        return m[1], m[2], ""
    return "", word, ""


# --- 対応 -----------------------------------------------------------------------------------------------------------


def bipartite(n_left, edges):
    """左の 0..n_left-1 をすべて右に対応させる。edges[i] は右の候補の並び（先のものを先に試す）。できなければ None"""
    match_right = {}

    def try_assign(i, seen):
        for j in edges[i]:
            if j in seen:
                continue
            seen.add(j)
            if j not in match_right or try_assign(match_right[j], seen):
                match_right[j] = i
                return True
        return False

    for i in range(n_left):
        if not try_assign(i, set()):
            return None
    return {i: j for j, i in match_right.items()}


def match_by_content(n_left, ranked_edges):
    """ranked_edges[i] は (順位, 右) の並び。順位 0 の辺だけで全部が対応すれば、それを使う（lock の md5 と今の中身が
    同じ項を、ほかの理由の対応より先にする）。できなければ、すべての辺で順位の順に試す"""
    best = [[j for rank, j in sorted(e) if rank == 0] for e in ranked_edges]
    found = bipartite(n_left, best)
    if found is not None:
        return found
    return bipartite(n_left, [[j for _, j in sorted(e)] for e in ranked_edges])


def toml_values_differ_only_by_paths(old, new, same_file):
    """2 つの TOML の値が、文字列の値のほかは同じで、違う文字列がどれも同じファイルを指すか。
    same_file(old の文字列, new の文字列) が文字列の対を判定する"""
    if isinstance(old, dict):
        if not isinstance(new, dict) or list(old) != list(new):
            return False
        return all(
            toml_values_differ_only_by_paths(old[k], new[k], same_file) for k in old
        )
    if isinstance(old, list):
        if not isinstance(new, list) or len(old) != len(new):
            return False
        return all(
            toml_values_differ_only_by_paths(a, b, same_file)
            for a, b in zip(old, new, strict=True)
        )
    if type(old) is not type(new):
        return False
    if isinstance(old, str) and old != new:
        return same_file(old, new)
    return old == new


def inputs_paths(summary):
    """比較のツールの出力の JSON の `inputs`（パスと md5 の並びか、パスから md5 への辞書）のパス"""
    if not isinstance(summary, dict) or "inputs" not in summary:
        return None
    inputs = summary["inputs"]
    if isinstance(inputs, dict):
        return list(inputs)
    if isinstance(inputs, list):
        return [i["path"] for i in inputs if isinstance(i, dict) and "path" in i]
    return None


def unwrap(value):
    """tomlkit の値を、YAML に書ける素の値にする"""
    if hasattr(value, "unwrap"):
        return value.unwrap()
    if isinstance(value, dict):
        return {k: unwrap(v) for k, v in value.items()}
    if isinstance(value, list):
        return [unwrap(v) for v in value]
    return value


# --- DVC と git の読み方 ----------------------------------------------------------------------------------------------


def git(root, *args, binary=False):
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, check=False
    )
    if result.returncode != 0:
        raise ToolError(
            f"git {' '.join(args)}: {result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout if binary else result.stdout.decode()


@dataclass
class Def:
    """展開したステージの定義。パスはリポジトリの直下から"""

    yaml: str
    name: str
    wdir: str
    cmd: object
    deps: list
    params: list  # (ファイル, キーの tuple)
    outs: list  # 定義の並び

    @property
    def lock(self):
        return join(dirname(self.yaml), LOCK)

    def key(self):
        return f"{self.yaml}:{self.name}"

    def same_definition(self, other):
        return (
            self.wdir,
            self.cmd,
            sorted(self.deps),
            sorted(self.params),
            self.outs,
        ) == (
            other.wdir,
            other.cmd,
            sorted(other.deps),
            sorted(other.params),
            other.outs,
        )


def read_defs(repo):
    from dvc.dependency.param import ParamsDependency
    from dvc.stage import PipelineStage

    defs = {}
    for st in repo.index.stages:
        if not isinstance(st, PipelineStage):
            continue
        wdir = posixpath.relpath(st.wdir, repo.root_dir)
        wdir = "" if wdir == "." else wdir
        deps, params = [], []
        for d in st.deps:
            if isinstance(d, ParamsDependency):
                keys = d.params if isinstance(d.params, list) else list(d.params or [])
                if not keys:
                    raise ToolError(
                        f"{st.path_in_repo}:{st.name}: キーを挙げない params（{d.def_path}）は扱わない"
                    )
                params.append((join(wdir, d.def_path), tuple(keys)))
            else:
                deps.append(join(wdir, d.def_path))
        yaml = st.path_in_repo
        defs[(yaml, st.name)] = Def(
            yaml=yaml,
            name=st.name,
            wdir=wdir,
            cmd=st.cmd,
            deps=deps,
            params=params,
            outs=[join(wdir, o.def_path) for o in st.outs],
        )
    return defs


class Side:
    """動かす前のコミットか、今の作業ツリー。md5 は DVC のライブラリで取る"""

    def __init__(self, root, repo):
        from dvc.stage import Stage

        self.root = root
        self.repo = repo
        self.defs = read_defs(repo)
        self._dummy = Stage(
            repo, path=repo.fs.join(repo.root_dir, "dvc.yaml"), wdir=repo.root_dir
        )
        self._md5 = {}
        # 追跡しない出力の md5（lock の outs と .dvc の outs）
        self.recorded = {}
        # ステージの定義が名指しするパス（作業ツリーに無い出力を、パスの語と読むため）
        self.known = {
            p
            for d in self.defs.values()
            for p in [*d.deps, *d.outs, *(f for f, _ in d.params)]
        }

    def file_md5(self, path):
        """ファイルかディレクトリの今の中身の md5（.dir を含む）。無ければ None"""
        from dvc.dependency import Dependency

        if path not in self._md5:
            try:
                self._md5[path] = Dependency(self._dummy, path).get_hash().value
            except FileNotFoundError:
                self._md5[path] = None
        return self._md5[path]

    def params(self, path, keys):
        from dvc.dependency.param import ParamsDependency

        try:
            return unwrap(ParamsDependency(self._dummy, path, list(keys)).read_params())
        except FileNotFoundError:
            return None

    def manifest_md5(self, dir_md5, relpath):
        """キャッシュの .dir の一覧から、ディレクトリの中のファイルの md5"""
        cached = self.repo.cache.local.oid_to_path(dir_md5)
        if not os.path.exists(cached):
            return None
        with open(cached, encoding="utf-8") as f:
            for entry in json.load(f):
                if entry.get("relpath") == relpath:
                    return entry.get("md5")
        return None

    def recorded_md5(self, path):
        if path in self.recorded:
            return self.recorded[path]
        for out, md5 in self.recorded.items():
            if md5 and md5.endswith(".dir") and path.startswith(out + "/"):
                return self.manifest_md5(md5, path[len(out) + 1 :])
        return None


class Base(Side):
    def __init__(self, root, rev):
        from dvc.repo import Repo

        self.rev = git(root, "rev-parse", "--verify", rev + "^{commit}").strip()
        super().__init__(root, Repo(str(root), rev=self.rev))
        self.files = set(
            git(root, "ls-tree", "-r", "-z", "--name-only", self.rev).split("\0")
        ) - {""}
        self.dirs = {""}
        for f in self.files:
            d = dirname(f)
            while d and d not in self.dirs:
                self.dirs.add(d)
                d = dirname(d)
        self.lock_texts = {
            f: self.blob(f) for f in self.files if posixpath.basename(f) == LOCK
        }

    def blob(self, path):
        if path not in self.files:
            return None
        return git(self.root, "show", f"{self.rev}:{path}", binary=True)

    def exists(self, path):
        return (
            path in self.files
            or path in self.dirs
            or path in self.recorded
            or path in self.known
        )

    def md5(self, path):
        if path in self.files or path in self.dirs:
            return self.file_md5(path)
        return self.recorded_md5(path)


class Current(Side):
    def __init__(self, root):
        from dvc.repo import Repo

        super().__init__(root, Repo(str(root)))

    def exists(self, path):
        return (
            os.path.lexists(self.root / path)
            or path in self.recorded
            or path in self.known
        )

    def md5(self, path):
        if os.path.lexists(self.root / path):
            return self.file_md5(path)
        return self.recorded_md5(path)

    def size(self, path):
        return os.path.getsize(self.root / path)


@dataclass
class Entry:
    """動かす前の lock の項"""

    lock: str
    name: str
    data: dict
    wdir: str
    base_def: object  # Def か None

    def key(self):
        return f"{self.lock}:{self.name}"

    def deps(self):
        return [(join(self.wdir, d["path"]), d) for d in self.data.get("deps", [])]

    def outs(self):
        return {join(self.wdir, o["path"]): o for o in self.data.get("outs", [])}

    def params(self):
        return {join(self.wdir, f): v for f, v in self.data.get("params", {}).items()}


# --- 判定 -----------------------------------------------------------------------------------------------------------


class Tool:
    def __init__(self, root, base_rev):
        from dvc.utils.serialize import parse_yaml

        self.root = root
        self.base = Base(root, base_rev)
        self.cur = Current(root)
        self.changed = set(
            git(root, "diff", "--name-only", "-z", self.base.rev).split("\0")
        ) | set(
            git(root, "ls-files", "--others", "--exclude-standard", "-z").split("\0")
        )
        self.changed.discard("")
        normalize = Path(__file__).resolve().parent / "normalize_toml.py"
        spec = importlib.util.spec_from_file_location("normalize_toml", normalize)
        self.normalize_toml = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.normalize_toml)

        self.base_locks = {
            p: parse_yaml(t.decode(), p) for p, t in self.base.lock_texts.items()
        }
        self.entries = {}
        for lock, data in self.base_locks.items():
            yaml = join(dirname(lock), "dvc.yaml")
            for name, entry in (data.get("stages") or {}).items():
                base_def = self.base.defs.get((yaml, name))
                wdir = base_def.wdir if base_def else dirname(lock)
                self.entries[(lock, name)] = Entry(lock, name, entry, wdir, base_def)
        for e in self.entries.values():
            for path, out in e.outs().items():
                self.base.recorded[path] = out.get("md5")
        self.base_dvc_files = {}
        for f in self.base.files:
            if f.endswith(".dvc") and posixpath.basename(f) != ".dvc":
                data = parse_yaml(self.base.blob(f).decode(), f)
                self.base_dvc_files[f] = data
                for out in data.get("outs", []):
                    self.base.recorded[join(dirname(f), out["path"])] = out.get("md5")
        self.cur_dvc_files = {}
        for f in self.tracked_or_new("*.dvc"):
            data = parse_yaml((root / f).read_text(), f)
            self.cur_dvc_files[f] = data
            for out in data.get("outs", []):
                self.cur.recorded[join(dirname(f), out["path"])] = out.get("md5")

        self.origin = self.match_locks()
        self.kept = set()
        for lock, origin in self.origin.items():
            yaml = join(dirname(lock), "dvc.yaml")
            for name in self.base_locks[origin].get("stages") or {}:
                if (yaml, name) in self.cur.defs:
                    self.kept.add((origin, name))
        self.consumed = set()
        self.pairs = set()  # 確かめたステージで組になった (動かす前, 今) のパス
        self.toml_ok = {}
        self.results = {}  # 今のステージの key → (動かす前の項, 新しい項)
        self.unverified = {}
        self.comparison = {}
        self.no_entry = []

    def tracked_or_new(self, pattern):
        out = git(
            self.root,
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
            "--",
            f":(glob)**/{pattern}",
            pattern,
        )
        return sorted(
            {f for f in out.split("\0") if f and os.path.lexists(self.root / f)}
        )

    def match_locks(self):
        """作業ツリーの dvc.lock ごとに、動かす前のどの lock か。同じパスで同じバイトか、ほかに同じバイトのものが 1 つ"""
        origin = {}
        for lock in self.tracked_or_new(LOCK):
            if posixpath.basename(lock) != LOCK:
                continue
            text = (self.root / lock).read_bytes()
            if lock in self.base.lock_texts:
                if self.base.lock_texts[lock] != text:
                    raise ToolError(
                        f"{lock} が動かす前のコミットと違う（ツールを 2 度実行したか、lock を手で変えた）。"
                        f"git checkout {self.base.rev[:9]} -- {lock} で戻してから実行する"
                    )
                origin[lock] = lock
                continue
            same = [p for p, t in self.base.lock_texts.items() if t == text]
            if len(same) != 1:
                raise ToolError(
                    f"{lock} と同じバイトの dvc.lock が動かす前のコミットに {len(same)} 個ある（1 個のときだけ、動かした lock と読む）"
                )
            origin[lock] = same[0]
        twice = [o for o in set(origin.values()) if list(origin.values()).count(o) > 1]
        if twice:
            raise ToolError(f"動かす前の {twice[0]} に当たる dvc.lock が 2 つある")
        return origin

    # 今の md5。作業ツリーに無い出力は、確かめたステージの新しい lock の md5
    def cur_md5(self, path):
        return self.cur.md5(path)

    def same_file(self, old, new, local_pairs):
        """old（動かす前）と new（今）が同じファイルか"""
        if old == new or (old, new) in local_pairs or (old, new) in self.pairs:
            return True
        base_md5 = self.base.md5(old)
        return base_md5 is not None and base_md5 == self.cur_md5(new)

    def toml_path_only(self, old, new, local_pairs):
        """動かす前の TOML（old）と今の TOML（new）の差が、同じファイルを指すパスの文字列だけか"""
        key = (old, new)
        if self.toml_ok.get(key):
            return True
        old_bytes = self.base.blob(old)
        new_path = self.root / new
        if old_bytes is None or not new_path.is_file():
            return False
        try:
            ov = tomllib.loads(old_bytes.decode())
            nv = tomllib.loads(new_path.read_text())
        except (tomllib.TOMLDecodeError, UnicodeDecodeError):
            return False

        def strings_same(a, b):
            return self.same_file(
                join(dirname(old), a), join(dirname(new), b), local_pairs
            )

        ok = toml_values_differ_only_by_paths(ov, nv, strings_same)
        if ok:
            self.toml_ok[key] = True
        return ok

    def normal_form_current(self, json_path, toml_path):
        """正規形のファイルが、今の TOML から書いたものと同じバイトか"""
        try:
            expected = self.normalize_toml.normalize(
                (self.root / toml_path).read_text()
            )
            return (self.root / json_path).read_text() == expected
        except (OSError, ValueError, TypeError):
            return False

    def toml_rewrite(self, old, new, local_pairs):
        """old から new への項が、入力の TOML（か正規形）の中のパスだけの書き換えで、確かめられたか"""
        src_old, src_new = toml_source(old), toml_source(new)
        if src_old is None or src_new is None or (old == src_old) != (new == src_new):
            return False
        if old == new and old not in self.changed and src_old not in self.changed:
            return False
        if self.base.md5(old) == self.cur_md5(new):
            return False
        if not self.toml_path_only(src_old, src_new, local_pairs):
            return False
        return new == src_new or self.normal_form_current(new, src_new)

    def rewritten_md5(self, item, old, new):
        """TOML の書き換えを確かめた項の新しい md5 と size。動かす前に lock が中身と合っていなかった項は持ち越す"""
        if item.get("md5") != self.base.md5(old):
            return item
        item = dict(item)
        item["md5"] = self.cur_md5(new)
        if "size" in item:
            item["size"] = self.cur.size(new)
        return item

    def match_deps(self, entry, s, out_pairs):
        lock_deps = entry.deps()
        if len(lock_deps) != len(s.deps):
            raise Unverified(
                f"deps の数が違う（lock {len(lock_deps)}、dvc.yaml {len(s.deps)}）"
            )
        result = []  # (old, item, new, how)
        new_left = list(s.deps)
        old_left = []
        for old, item in lock_deps:
            if old in new_left:
                new_left.remove(old)
                result.append((old, item, old, "path"))
            else:
                old_left.append((old, item))
        edges = []
        for old, item in old_left:
            ranked = []
            for j, new in enumerate(new_left):
                cm = self.cur_md5(new)
                if cm is not None and item.get("md5") == cm:
                    ranked.append((0, j))
                elif cm is not None and self.base.md5(old) == cm:
                    ranked.append((1, j))
                elif self.toml_rewrite(old, new, out_pairs):
                    ranked.append((2, j))
            edges.append(ranked)
        found = match_by_content(len(old_left), edges)
        if found is None:
            missing = ", ".join(o for o, _ in old_left)
            raise Unverified(
                f"deps の {missing} と中身の同じファイルが、今の deps に見つからない"
            )
        for i, (old, item) in enumerate(old_left):
            result.append((old, item, new_left[found[i]], "content"))
        out = []
        for old, item, new, how in result:
            if self.toml_rewrite(old, new, out_pairs):
                out.append((old, self.rewritten_md5(item, old, new), new, True))
            else:
                out.append((old, item, new, False))
        return out

    def match_params(self, entry, s, local_pairs):
        lock_params = entry.params()
        if len(lock_params) != len(s.params):
            raise Unverified("params のファイルの数が違う")
        cur = {f: keys for f, keys in s.params}
        result = []  # (old, values, new)
        old_left = []
        for old, values in lock_params.items():
            if old in cur:
                result.append((old, values, old))
            else:
                old_left.append((old, values))
        new_left = [f for f in cur if f not in lock_params]
        edges = []
        for old, values in old_left:
            ranked = []
            for j, new in enumerate(new_left):
                now = self.cur.params(new, cur[new])
                if now is not None and now == values:
                    ranked.append((0, j))
                elif now is not None and self.base.params(old, cur[new]) == now:
                    ranked.append((1, j))
                elif self.toml_rewrite(old, new, local_pairs):
                    ranked.append((2, j))
            edges.append(ranked)
        found = match_by_content(len(old_left), edges)
        if found is None:
            raise Unverified(
                "params の "
                + ", ".join(o for o, _ in old_left)
                + " と値の同じファイルが、今の params に見つからない"
            )
        for i, (old, values) in enumerate(old_left):
            result.append((old, values, new_left[found[i]]))
        out = []
        for old, values, new in result:
            keys = cur[new]
            if set(values) != set(keys):
                raise Unverified(f"params の {new} のキーが lock と違う")
            if self.toml_rewrite(old, new, local_pairs):
                before = self.base.params(old, keys) or {}
                now = self.cur.params(new, keys)
                values = {
                    k: (now[k] if before.get(k) == v and k in now else v)
                    for k, v in values.items()
                }
                out.append((old, values, new, True))
            else:
                out.append((old, values, new, False))
        return out

    def check_cmd(self, entry, s, pairs):
        base_cmd = entry.base_def.cmd if entry.base_def else entry.data.get("cmd")
        old_parts, new_parts = cmd_parts(base_cmd), cmd_parts(s.cmd)
        if len(old_parts) != len(new_parts):
            raise Unverified("cmd の要素の数が違う")
        for old_words, new_words in zip(old_parts, new_parts, strict=True):
            if len(old_words) != len(new_words):
                raise Unverified(
                    f"cmd の語の数が違う（{len(old_words)} と {len(new_words)}）"
                )
            for wo, wn in zip(old_words, new_words, strict=True):
                if not self.same_word(wo, wn, entry.wdir, s.wdir, pairs):
                    raise Unverified(
                        f"cmd の語 {wo!r} と {wn!r} が、同じファイルを指すと確かめられない"
                    )
        return base_cmd

    def same_word(self, wo, wn, old_wdir, new_wdir, pairs):
        if wo == wn and old_wdir == new_wdir:
            return True
        (po_pre, po, po_suf), (pn_pre, pn, pn_suf) = split_word(wo), split_word(wn)
        if po_pre != pn_pre or po_suf != pn_suf:
            return False
        if GLOB_CHARS & set(po + pn):
            return False
        old, new = join(old_wdir, po), join(new_wdir, pn)
        if old == new or (old, new) in pairs:
            return True
        if old in self.base.files and self.base.md5(old) == self.cur_md5(new):
            return True
        # 同じ文字列で、どちらの置き場にもファイルが無い語は、パスでない語（オプションや引数の値）
        return wo == wn and not self.base.exists(old) and not self.cur.exists(new)

    def read_json_out(self, item, paths):
        """lock の md5 の中身を、作業ツリー、動かす前のコミット（git で追跡する出力）、キャッシュの順に探す"""
        md5 = item.get("md5")
        for p in paths:
            if os.path.isfile(self.root / p) and self.cur.file_md5(p) == md5:
                return (self.root / p).read_text()
        if paths[0] in self.base.files and self.base.md5(paths[0]) == md5:
            return self.base.blob(paths[0]).decode()
        cached = self.cur.repo.cache.local.oid_to_path(md5) if md5 else None
        if cached and os.path.isfile(cached):
            with open(cached, encoding="utf-8") as f:
                return f.read()
        raise Unverified(
            f"出力の {paths[-1]} を、作業ツリーにもキャッシュにも見つけられない"
        )

    def comparison_reason(self, entry, s, out_pairs, deps):
        """出力の JSON の `inputs` が、動かした入力か中のパスを書き換えた入力を持てば、その説明"""
        dep_map = {old: (new, rewritten) for old, _, new, rewritten in deps}
        old_dir, new_dir = dirname(entry.lock), dirname(s.yaml)
        lock_outs = entry.outs()
        for old, new in out_pairs:
            item = lock_outs[old]
            if not old.endswith(".json") or str(item.get("md5", "")).endswith(".dir"):
                continue
            try:
                summary = json.loads(self.read_json_out(item, [old, new]))
            except json.JSONDecodeError:
                continue
            paths = inputs_paths(summary)
            for s_old in paths or []:
                p = join(old_dir, s_old)
                dep = p
                if p.endswith(".toml"):
                    dep = TOML_FORM[0] + p + TOML_FORM[1]
                if dep not in dep_map:
                    continue
                new_dep, rewritten = dep_map[dep]
                source = toml_source(new_dep) if p.endswith(".toml") else new_dep
                if rewritten or rel(source, new_dir) != s_old:
                    return f"{new} の inputs が、動かした入力 {s_old} を持つ"
        return None

    def verify(self, entry, s):
        """動かす前の lock の項を今のステージに対応させた新しい項。確かめられなければ Unverified"""
        lock_outs = entry.outs()
        if entry.base_def is not None:
            order = entry.base_def.outs
        elif len(lock_outs) <= 1:
            order = list(lock_outs)
        else:
            raise Unverified(
                "動かす前の dvc.yaml に定義が無く、outs の並びが分からない"
            )
        if sorted(order) != sorted(lock_outs):
            raise Unverified("lock の outs が、動かす前の定義の outs と違う")
        if len(order) != len(s.outs):
            raise Unverified("outs の数が違う")
        out_pairs = set(zip(order, s.outs, strict=True))
        deps = self.match_deps(entry, s, out_pairs)
        local = out_pairs | {(o, n) for o, _, n, _ in deps}
        params = self.match_params(entry, s, local)
        local |= {(o, n) for o, _, n, _ in params}
        for old, new in list(local):
            so, sn = normal_form_source(old), normal_form_source(new)
            if so is not None and sn is not None:
                local.add((so, sn))
        new_outs = []
        for old, new in zip(order, s.outs, strict=True):
            item = dict(lock_outs[old])
            if self.toml_rewrite(old, new, local):
                item = self.rewritten_md5(item, old, new)
            elif (
                old != new
                and os.path.lexists(self.root / new)
                and self.cur.file_md5(new) != item.get("md5")
            ):
                raise Unverified(
                    f"作業ツリーの {new} の md5 が、持ち越す lock の md5 と違う"
                )
            item["path"] = rel(new, s.wdir)
            new_outs.append((new, item))
        base_cmd = self.check_cmd(entry, s, local)

        data = copy.deepcopy(dict(entry.data))
        data["cmd"] = (
            s.cmd if entry.data.get("cmd") == base_cmd else entry.data.get("cmd")
        )
        if deps or "deps" in data:
            data["deps"] = sorted(
                ({**item, "path": rel(new, s.wdir)} for _, item, new, _ in deps),
                key=lambda d: d["path"],
            )
        if params or "params" in data:
            data["params"] = dict(
                sorted((rel(new, s.wdir), values) for _, values, new, _ in params)
            )
        if new_outs or "outs" in data:
            data["outs"] = sorted(
                (item for _, item in new_outs), key=lambda o: o["path"]
            )
        reason = self.comparison_reason(entry, s, sorted(out_pairs), deps)
        return data, local, new_outs, reason

    def candidates(self, s):
        """名前で対応の付かないステージの候補（どのステージにも残らず、まだ使っていない項）"""
        for key, entry in self.entries.items():
            if key in self.kept or key in self.consumed:
                continue
            if len(entry.data.get("deps", [])) != len(s.deps):
                continue
            if len(entry.data.get("outs", [])) != len(s.outs):
                continue
            if len(entry.data.get("params", {})) != len(s.params):
                continue
            yield entry

    def process(self, s):
        """1 つのステージを判定する。確かめたら True"""
        origin = self.origin.get(s.lock)
        if origin is not None and (origin, s.name) in self.entries:
            entry = self.entries[(origin, s.name)]
            if (
                entry.base_def is not None
                and entry.base_def.yaml == s.yaml
                and entry.base_def.same_definition(s)
                and not self.touches_changed_toml(s)
            ):
                self.accept(s, entry, entry.data, set(), None)
                return True
            try:
                data, local, new_outs, reason = self.verify(entry, s)
            except Unverified as e:
                self.unverified[s.key()] = str(e)
                return False
            self.accept(s, entry, data, local, reason, new_outs)
            return True
        unrun_yaml = join(dirname(origin), "dvc.yaml") if origin else s.yaml
        if (unrun_yaml, s.name) in self.base.defs:
            # 動かす前から lock の項が無い（流していない）ステージ
            return True
        found, reasons = [], []
        for entry in self.candidates(s):
            try:
                found.append((entry, *self.verify(entry, s)))
            except Unverified as e:
                reasons.append(f"{entry.key()}: {e}")
        if not found:
            if reasons or any(
                k not in self.kept and k not in self.consumed for k in self.entries
            ):
                self.unverified[s.key()] = (
                    "名前の同じ lock の項が無く、確認の通る項も無い"
                    + (f"（{'; '.join(reasons[:3])}）" if reasons else "")
                )
                return False
            self.no_entry.append(s.key())
            return True
        datas = {json.dumps(f[1], sort_keys=True, default=str) for f in found}
        if len(datas) > 1:
            self.unverified[s.key()] = (
                "確認の通る lock の項が複数あり、書き換えた結果が違う: "
                + ", ".join(f[0].key() for f in found)
            )
            return False
        entry, data, local, new_outs, reason = found[0]
        self.consumed.add((entry.lock, entry.name))
        self.accept(s, entry, data, local, reason, new_outs)
        return True

    def touches_changed_toml(self, s):
        for path in [*s.deps, *(f for f, _ in s.params), *s.outs]:
            source = toml_source(path)
            if source is not None and (path in self.changed or source in self.changed):
                return True
        return False

    def accept(self, s, entry, data, local, reason, new_outs=None):
        self.unverified.pop(s.key(), None)
        if reason is not None:
            self.comparison[s.key()] = reason
            return
        self.results[s.key()] = (entry, data, s)
        self.pairs |= local
        if new_outs is None:
            new_outs = [(join(s.wdir, o["path"]), o) for o in data.get("outs", [])]
        for path, item in new_outs:
            self.cur.recorded[path] = item.get("md5")

    def order(self):
        """今のステージを上流から並べる（出力のパスか、その上か下のディレクトリを deps に持つステージが下流）"""
        defs = {d.key(): d for d in self.cur.defs.values()}
        producer, under = {}, {}
        for key, d in defs.items():
            for out in d.outs:
                producer.setdefault(out, set()).add(key)
                parent = dirname(out)
                while parent:
                    under.setdefault(parent, set()).add(key)
                    parent = dirname(parent)
        upstream = {}
        for key, d in defs.items():
            ups = set()
            for dep in [*d.deps, *(f for f, _ in d.params)]:
                ups |= under.get(dep, set())
                path = dep
                while path:
                    ups |= producer.get(path, set())
                    path = dirname(path)
            upstream[key] = ups - {key}
        done, ordered = set(), []

        def visit(key, stack):
            if key in done:
                return
            if key in stack:
                raise ToolError(f"ステージの依存が循環する: {key}")
            stack.add(key)
            for up in sorted(upstream[key]):
                visit(up, stack)
            stack.discard(key)
            done.add(key)
            ordered.append(defs[key])

        for key in sorted(defs):
            visit(key, set())
        return ordered

    def run(self):
        """上流から判定する。確かめたステージで組になったパスが増えたら、始めから判定し直す（入力の TOML が名指しする
        出力の組は、その TOML を読む下流のステージの outs の対応で決まるので、上流の正規形のステージの判定に後から効く）"""
        stages = self.order()
        initial = dict(self.cur.recorded)
        while True:
            known = set(self.pairs)
            self.results, self.unverified, self.comparison = {}, {}, {}
            self.no_entry, self.consumed = [], set()
            self.cur.recorded = dict(initial)
            pending = stages
            while pending:
                failed = [s for s in pending if not self.process(s)]
                if len(failed) == len(pending):
                    break
                pending = failed
            if self.pairs == known:
                break
        self.check_dvc_files()

    def check_dvc_files(self):
        """動かした .dvc は、動かす前の .dvc のどれかと、outs の path のほかが同じであること"""
        self.dvc_file_problems = []
        unchanged = {
            f for f, d in self.cur_dvc_files.items() if self.base_dvc_files.get(f) == d
        }

        def without_paths(data):
            data = copy.deepcopy(dict(data))
            for out in data.get("outs", []):
                out.pop("path", None)
            return data

        for f, data in self.cur_dvc_files.items():
            if f in unchanged:
                continue
            same_md5 = [
                g
                for g, b in self.base_dvc_files.items()
                if g not in unchanged
                and [o.get("md5") for o in b.get("outs", [])]
                == [o.get("md5") for o in data.get("outs", [])]
            ]
            if not same_md5:
                continue
            if not any(
                without_paths(self.base_dvc_files[g]) == without_paths(data)
                for g in same_md5
            ):
                self.dvc_file_problems.append(
                    f"{f}: 中身の同じ動かす前の .dvc（{', '.join(same_md5)}）と、outs の path のほかが違う"
                )

    def write(self):
        """確かめたステージの新しい項を、今の dvc.yaml の隣の dvc.lock に書く。書いた lock のパスを返す"""
        from dvc.utils.serialize import dumps_yaml

        out = {}
        for lock, origin in self.origin.items():
            data = copy.deepcopy(dict(self.base_locks[origin]))
            data["stages"] = dict(data.get("stages") or {})
            out[lock] = data
        holder = {origin: lock for lock, origin in self.origin.items()}
        for key in sorted(self.results):
            entry, data, s = self.results[key]
            source = holder.get(entry.lock)
            if source == s.lock:
                # 同じ lock の中では、元の項の位置に置く
                out[s.lock]["stages"] = {
                    (s.name if name == entry.name else name): (
                        data if name == entry.name else value
                    )
                    for name, value in out[s.lock]["stages"].items()
                }
                continue
            if source is not None:
                out[source]["stages"].pop(entry.name, None)
            target = out.setdefault(s.lock, {"schema": "2.0", "stages": {}})
            target["stages"][s.name] = data
        written = []
        for lock, data in sorted(out.items()):
            path = self.root / lock
            if not data["stages"]:
                if path.exists():
                    path.unlink()
                    written.append(lock)
                continue
            text = dumps_yaml(data)
            if not path.exists() or path.read_text() != text:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                written.append(lock)
        return written

    def rewritten(self):
        return sorted(
            key
            for key, (entry, data, s) in self.results.items()
            if data != entry.data
            or entry.name != s.name
            or self.origin.get(s.lock) != entry.lock
        )

    def left(self):
        used = {(e.lock, e.name) for e, _, _ in self.results.values()}
        return sorted(
            f"{lock}:{name}"
            for lock, name in self.entries
            if (lock, name) not in self.kept | self.consumed | used
        )


def main(argv):
    parser = argparse.ArgumentParser(
        prog="scripts/rewrite-lock-paths.py",
        description=__doc__.splitlines()[0],
    )
    parser.add_argument(
        "--base", default="HEAD", help="動かす前のコミット（既定は HEAD）"
    )
    args = parser.parse_args(argv[1:])
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel").strip())
    tool = Tool(root, args.base)
    tool.run()
    rewritten = tool.rewritten()
    tool.write()
    for key in rewritten:
        print(f"REWRITTEN: {key}（← {tool.results[key][0].key()}）")
    for key in sorted(tool.comparison):
        print(f"COMPARISON: {key}: {tool.comparison[key]}")
    for key in sorted(tool.no_entry):
        print(f"NO-ENTRY: {key}")
    for key in tool.left():
        print(f"LEFT: {key}")
    for key in sorted(tool.unverified):
        print(f"UNVERIFIED: {key}: {tool.unverified[key]}")
    for problem in tool.dvc_file_problems:
        print(f"DVC-FILE: {problem}")
    return 1 if tool.unverified or tool.dvc_file_problems else 0


if __name__ == "__main__":
    try:
        import dvc.repo  # noqa: F401
    except ModuleNotFoundError:
        try:
            reexec_under_dvc(sys.argv)
        except ToolError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)
    try:
        sys.exit(main(sys.argv))
    except ToolError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
