#!/usr/bin/env python3
"""形式を変える修正の検査（ADR-0258）。xlc の実行ファイルを呼ぶステージが deps にする上流の出力を、ブランチのコードの
読み手（crates/outputs の read-outputs）で読み、読めない出力とそれを書いたステージを出す。

使い方（default の devShell で、リポジトリの中から）:
    scripts/check-readable-outputs.py [--base <ref>] [--all] [--cache <dir>] [--reader <実行ファイル>]

この検査は、再実行が要るかを決めない（DVC が changed で示す）。答えるのは、修正を入れた後の実行ファイルが、
dvc.lock の指す古い上流の出力をまだ読めるか（changed をためたまま下流の作業を続けられるか）だけである。

- 読むのは、deps に `binary-forms/` の実行ファイルの正規形を持つステージ（xlc の実行ファイルを呼ぶステージ）の deps の
  うち、例の dvc.lock のステージの outs（ディレクトリの outs の中のファイルと、deps のディレクトリの下にある outs を含む）と
  `.dvc` の生成物であるもののうち、xlc の形式のもの（`.msgpack` と、`xs range` の書く範囲の TOML）。
  直下の dvc.lock（build と正規形）の出力、md5 の無い（流していない）出力、xlc の形式でない出力は読まない
- 中身は作業ツリーでなく DVC のキャッシュから md5 で読む。出力を戻していない worktree でも同じ結果になり、走っている
  ジョブが書いている途中の出力には当たらない。キャッシュに無く、git の index が追う出力（`cache: false` の出力。DVC は
  git が追う出力をキャッシュに入れない（`dvc repro` が OutputAlreadyTrackedError で止まる））は、index の blob の md5 が lock の
  md5 と同じならその blob を読み、違えば `MISMATCH:` で止める（lock が古いか、出力の変更をステージしていない）。HEAD でなく
  index を読むのは、マージの hook の中では HEAD が既定ブランチのままで、index と lock がマージの結果になるため。どちらにも
  無い出力は `MISSING:` で止める
- 実行ファイルを変えうるブランチ（`--base` との差分（移動の元と先、追跡していないファイルを含む）が直下の dvc.yaml の `build` の deps、この script、下の一覧に当たる）
  だけを読み、ほかは `SKIP:` で終わる（`--all` は差分によらず読む）。`--base` の既定は HEAD と main の分岐点。
  一覧の形の確かめは、差分によらず行う
- 読めない出力を残してマージするときは、`scripts/unreadable-outputs.toml` の `[[output]]` に `md5`（ファイルの md5 か、
  ディレクトリの outs の `.dir` の md5）、`issue`（上流の作り直しの見積もりと時期を書いた open の issue。`#NNNN`）、
  `path`（結び付けたときのパス。照合には使わない）を書く。open の issue に結び付けた出力は `ALLOWED:` で通し、
  issue が無いか閉じていれば `STOP:` で止める
- 出力の行: `UNREADABLE:`、`ALLOWED:`、`MISSING:`、`MISMATCH:`、`STOP:`、`SKIP:`、`READ:`（読んだファイルの数、時間、xlc の形式でないので読まなかった出力の数）。終了コードは、
  `UNREADABLE:`、`MISSING:`、`MISMATCH:`、`STOP:` のどれかがあれば 1、ほかは 0
- 読み手の既定は `cargo run -p xlc-outputs`（dev のプロファイル。run-checks.sh の cargo test が作ったものを再利用する）
"""

import argparse
import dataclasses
import hashlib
import json
import posixpath
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import tomllib
import yaml

EXEMPTIONS = "scripts/unreadable-outputs.toml"
# これらを変えたブランチも読む（検査の中身と、通す出力が変わるため）
CHECK_FILES = ("scripts/check-readable-outputs.py", EXEMPTIONS)
BINARY_FORMS = "binary-forms/"
PACKAGE = ["--quiet", "--locked", "-p", "xlc-outputs", "--bin", "read-outputs"]
ISSUE = re.compile(r"^#(\d{4})$")
MD5 = re.compile(r"^[0-9a-f]{32}(\.dir)?$")
# xs の範囲のファイルを書くコマンド（`xs [--threads <数>] range <入力> <範囲>`）
XS_RANGE = re.compile(r"(^|[\s/])bin/xs(\s+--threads\s+\S+)?\s+range(\s|$)")


class Stop(Exception):
    pass


@dataclass(frozen=True)
class Output:
    """dvc.lock か .dvc の 1 つの outs"""

    path: str
    md5: str
    # dvc.lock の `hash: md5`（キャッシュの files/md5 の下）。無ければ DVC 2 の置き場（キャッシュの直下）
    hash: str | None
    writer: str
    # 書いたのが `xs range` のステージか
    by_xs_range: bool


@dataclass(frozen=True)
class Target:
    """読む 1 つのファイル"""

    path: str
    md5: str
    # ディレクトリの outs の中のファイルなら、そのディレクトリの md5
    dir_md5: str | None
    writer: str
    file: Path
    # 読み手に渡す種類。None は xlc の形式でない（読まない）
    kind: str | None


def git(root, *args):
    return subprocess.run(
        ["git", "-C", str(root), *args], check=True, capture_output=True, text=True
    ).stdout


def join(dir_, path):
    return posixpath.normpath(posixpath.join(dir_, path))


def load_yaml(root, path):
    # 766 の .dvc と 23 の lock を読むので、C の読み手にする（Python の読み手の数倍速い）
    with open(root / path, encoding="utf-8") as f:
        return yaml.load(f, Loader=yaml.CSafeLoader) or {}


def tracked(root, pattern):
    return git(root, "ls-files", "-z", "--", pattern).split("\0")[:-1]


@dataclass(frozen=True)
class Stage:
    lock: str
    name: str
    cmd: str
    # deps のパス（リポジトリの直下から）
    deps: tuple
    outs: tuple

    @property
    def by_xlc(self):
        """xlc の実行ファイルを呼ぶか（deps に実行ファイルの正規形を持つ）"""
        return any(p.startswith(BINARY_FORMS) for p in self.deps)


def stages(root):
    """例の dvc.lock のステージ（直下の dvc.lock の build と正規形は計算の出力を書かない）"""
    found = []
    for lock in tracked(root, "*dvc.lock"):
        if lock == "dvc.lock":
            continue
        d = posixpath.dirname(lock)
        for name, stage in (load_yaml(root, lock).get("stages") or {}).items():
            deps = tuple(join(d, dep["path"]) for dep in stage.get("deps") or [])
            cmd = stage.get("cmd") or ""
            cmd = " && ".join(cmd) if isinstance(cmd, list) else cmd
            found.append(Stage(lock, name, cmd, deps, tuple(stage.get("outs") or [])))
    return found


def locked(out, where, writer, by_xs_range, found):
    """md5 のある outs を found に足す"""
    md5 = out.get("md5")
    if md5 is None:
        return
    if not isinstance(md5, str) or not MD5.match(md5):
        raise Stop(f"{where} の {out.get('path')} の md5 を読めない: {md5!r}")
    path = join(posixpath.dirname(where), out["path"])
    found[path] = Output(path, md5, out.get("hash"), writer, by_xs_range)


def outputs(root, all_stages):
    """例の dvc.lock と .dvc の、md5 のある outs（パスはリポジトリの直下から）"""
    found = {}
    for stage in all_stages:
        writer = f"{posixpath.dirname(stage.lock) or '.'}/dvc.yaml:{stage.name}"
        for out in stage.outs:
            locked(out, stage.lock, writer, bool(XS_RANGE.search(stage.cmd)), found)
    for dvc_file in tracked(root, "*.dvc"):
        for out in load_yaml(root, dvc_file).get("outs") or []:
            locked(out, dvc_file, dvc_file, False, found)
    return found


def kind_of(path, by_xs_range):
    """読み手に渡す種類。xlc の形式は、形式の名前を持つ MessagePack と、`xs range` の書く範囲の TOML だけ。

    MessagePack は書き手（Python のステージを含む）によらず読む: xlc の読み手は xlc の形式の MessagePack しか受けないので、
    xlc のステージが deps にする MessagePack は xlc の形式のはずである。ほかの TOML（xlc の書く要約、Python が書く入力）に
    Rust の形式の読み手は無い
    """
    if path.endswith(".msgpack"):
        return "msgpack"
    if path.endswith(".toml") and by_xs_range:
        return "xs-range"
    return None


def cache_file(cache, md5, hash_):
    base = cache / "files" / "md5" if hash_ == "md5" else cache
    return base / md5[:2] / md5[2:]


def expand(out, cache, inside, targets, missing):
    """`out` のファイル（ディレクトリなら、その中で `inside` の下にあるもの）を targets に足す"""
    if not out.md5.endswith(".dir"):
        file = cache_file(cache, out.md5, out.hash)
        kind = kind_of(out.path, out.by_xs_range)
        targets[out.path] = Target(out.path, out.md5, None, out.writer, file, kind)
        return
    manifest = cache_file(cache, out.md5, out.hash)
    if not manifest.exists():
        missing.append((out.path, out.md5, f"DVC のキャッシュに無い（{manifest}）"))
        return
    for entry in json.loads(manifest.read_text()):
        relpath = entry["relpath"]
        if inside is None or relpath == inside or relpath.startswith(inside + "/"):
            path = f"{out.path}/{relpath}"
            file = cache_file(cache, entry["md5"], out.hash)
            kind = kind_of(path, out.by_xs_range)
            targets[path] = Target(path, entry["md5"], out.md5, out.writer, file, kind)


def collect(root, cache):
    """読むファイル（パスごと）と、キャッシュに無いディレクトリ（パス、md5、理由）の一覧"""
    all_stages = stages(root)
    found = outputs(root, all_stages)
    deps = {p for stage in all_stages if stage.by_xlc for p in stage.deps}
    targets, missing = {}, []
    for dep in sorted(deps):
        if dep in found:
            expand(found[dep], cache, None, targets, missing)
            continue
        container = next(
            (
                p
                for p in found
                if dep.startswith(p + "/") and found[p].md5.endswith(".dir")
            ),
            None,
        )
        if container is not None:
            expand(found[container], cache, dep[len(container) + 1 :], targets, missing)
            continue
        for path in found:
            if path.startswith(dep + "/"):
                expand(found[path], cache, None, targets, missing)
    return targets, missing


def issue_status(root, number):
    for path in sorted((root / "docs" / "issues").glob(f"{number}-*.md")):
        match = re.search(
            r"^status:\s*(\S+)\s*$", path.read_text(encoding="utf-8"), re.MULTILINE
        )
        return match[1] if match else None
    return None


def exemptions(root):
    """一覧の md5 から issue。一覧の誤りは Stop"""
    path = root / EXEMPTIONS
    if not path.exists():
        return {}
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as e:
        raise Stop(f"{EXEMPTIONS} を読めない: {e}") from e
    if set(data) - {"output"}:
        raise Stop(
            f"{EXEMPTIONS} の項目は [[output]] だけ: {sorted(set(data) - {'output'})}"
        )
    entries = data.get("output", [])
    if not isinstance(entries, list) or not all(isinstance(e, dict) for e in entries):
        raise Stop(f"{EXEMPTIONS} の output は [[output]] の表の並び")
    allowed = {}
    for i, entry in enumerate(entries, start=1):
        where = f"{EXEMPTIONS} の {i} 番目の [[output]]"
        if set(entry) != {"md5", "issue", "path"}:
            raise Stop(f"{where} の項目は md5、issue、path: {sorted(entry)}")
        if not all(isinstance(v, str) for v in entry.values()):
            raise Stop(f"{where} の md5、issue、path は文字列: {entry}")
        if not MD5.match(entry["md5"]):
            raise Stop(f"{where} の md5 が md5 の形でない: {entry['md5']!r}")
        match = ISSUE.match(entry["issue"])
        if match is None:
            raise Stop(f"{where} の issue が #NNNN の形でない: {entry['issue']!r}")
        status = issue_status(root, match[1])
        if status != "open":
            state = "無い" if status is None else f"閉じている（status: {status}）"
            raise Stop(
                f"{where} の {entry['issue']} が {state}。上流を作り直したなら、その項を消す"
            )
        if entry["md5"] in allowed:
            raise Stop(
                f"{where} の md5 {entry['md5']} は前の項にもある（issue は 1 つにする）"
            )
        allowed[entry["md5"]] = entry["issue"]
    return allowed


def git_bytes(root, *args):
    """git の標準出力のバイト。失敗は Stop"""
    result = subprocess.run(
        ["git", "-C", str(root), "--literal-pathspecs", *args],
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise Stop(
            f"git {' '.join(args)} が失敗した: {result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


def from_index(root, xlc, tmp):
    """キャッシュに無いファイルを、git の index が追っていれば、lock の md5 と同じ blob を tmp に書いたものに替える。
    読むファイル、どちらにも無いファイル（パス、md5、理由）、blob の md5 が lock と違うファイルの組。

    HEAD でなく index を読むのは、マージの hook の中では HEAD が既定ブランチのままで、index と作業ツリーの lock が
    マージの結果になるため（ステージしていない作業ツリーの変更は読まない）"""
    present, missing, mismatched = [], [], []
    for t in xlc:
        if t.file.exists():
            present.append(t)
            continue
        entries = git_bytes(root, "ls-files", "-s", "-z", "--", t.path).split(b"\0")[
            :-1
        ]
        if not entries:
            missing.append(
                (
                    t.path,
                    t.md5,
                    f"DVC のキャッシュに無く（{t.file}）、git の index も追っていない",
                )
            )
            continue
        mode, oid, stage = entries[0].split(b"\t", 1)[0].decode().split()
        if len(entries) != 1 or stage != "0":
            raise Stop(f"{t.path} は git の index で衝突している（マージを解いてから）")
        if mode not in ("100644", "100755"):
            raise Stop(
                f"{t.path} は git の index で通常のファイルでない（mode {mode}）ので、中身を lock の md5 と照らせない"
            )
        blob = git_bytes(root, "cat-file", "blob", oid)
        actual = hashlib.md5(blob).hexdigest()
        if actual != t.md5:
            mismatched.append((t.path, t.md5, actual))
            continue
        file = tmp / t.md5
        file.write_bytes(blob)
        present.append(dataclasses.replace(t, file=file))
    return present, missing, mismatched


def needs_reading(root, base):
    """--base との差分（作業ツリー）が、実行ファイルを変えうるパスに当たるか。当たったパスか None"""
    build = load_yaml(root, "dvc.yaml").get("stages", {}).get("build")
    if build is None:
        raise Stop("直下の dvc.yaml に build のステージが無い")
    watched = [*build.get("deps", []), *CHECK_FILES]
    # 移動は元と先の両方を見る（crates/ の外へ動かしても実行ファイルは変わる）。追跡していないファイルもビルドに入る
    changed = git(root, "diff", "--name-only", "--no-renames", "-z", base)
    untracked = git(root, "ls-files", "--others", "--exclude-standard", "-z")
    for path in (changed + untracked).split("\0")[:-1]:
        if any(path == w or path.startswith(w.rstrip("/") + "/") for w in watched):
            return path
    return None


def read(root, reader, targets):
    """読み手の UNREADABLE の行から、読めなかったファイルと理由"""
    lines = "".join(
        f"{t.kind}\t{t.file}\n" for t in {t.file: t for t in targets}.values()
    )
    result = subprocess.run(
        reader, cwd=root, input=lines, capture_output=True, text=True, check=False
    )
    failed = {}
    for line in result.stdout.splitlines():
        kind, sep, rest = line.partition("\t")
        file, sep2, reason = rest.partition("\t")
        if kind != "UNREADABLE" or not sep or not sep2:
            raise Stop(f"読み手の出力の行を読めない: {line!r}")
        failed[file] = reason
    if result.returncode not in (0, 1) or (result.returncode == 1) != bool(failed):
        raise Stop(
            f"読み手（{' '.join(reader)}）が終了コード {result.returncode} で終わった:\n"
            f"{result.stderr.strip()}"
        )
    return failed


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--base", help="比べる ref（既定は HEAD と main の分岐点）")
    parser.add_argument("--all", action="store_true", help="差分によらず読む")
    parser.add_argument(
        "--cache", type=Path, help="DVC のキャッシュ（既定は dvc cache dir）"
    )
    parser.add_argument(
        "--reader", help="読み手の実行ファイル（既定は cargo の xlc-outputs）"
    )
    args = parser.parse_args(argv[1:])
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel").strip())
    try:
        allowed = exemptions(root)
        if not args.all:
            base = args.base or git(root, "merge-base", "HEAD", "main").strip()
            hit = needs_reading(root, base)
            if hit is None:
                print(
                    f"SKIP: {base} からの差分は実行ファイルを変えない（直下の dvc.yaml の build の deps）"
                )
                return 0
            print(f"RUN: {hit} を変えたので、上流の出力を読む")
        cache = args.cache
        if cache is None:
            out = subprocess.run(
                ["dvc", "cache", "dir"],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
            )
            if out.returncode != 0:
                raise Stop(f"dvc cache dir が失敗した: {out.stderr.strip()}")
            cache = Path(out.stdout.strip())
        targets, missing_dirs = collect(root, cache)
        xlc = [t for t in targets.values() if t.kind is not None]
        if args.reader:
            reader = [args.reader]
        else:
            # 読む時間にビルドを混ぜないため、先に作る
            built = subprocess.run(
                ["cargo", "build", *PACKAGE],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
            )
            if built.returncode != 0:
                raise Stop(f"読み手をビルドできない:\n{built.stderr.strip()}")
            reader = ["cargo", "run", *PACKAGE]
        # index の blob を書いたファイルは、読み終えるまで残す
        with tempfile.TemporaryDirectory() as tmp:
            present, missing_files, mismatched = from_index(root, xlc, Path(tmp))
            start = time.monotonic()
            failed = read(root, reader, present) if present else {}
            seconds = time.monotonic() - start
        missing = missing_dirs + missing_files
    except Stop as e:
        print(f"STOP: {e}")
        return 1
    bad = bool(missing or mismatched)
    for path, md5, reason in sorted(missing):
        print(f"MISSING: {path} (md5 {md5}): {reason}")
    for path, md5, actual in sorted(mismatched):
        print(
            f"MISMATCH: {path} (md5 {md5}): DVC のキャッシュに無く、git の index の版の md5 は {actual}。"
            "lock が古いか、出力の変更をステージしていない"
        )
    for t in sorted(present, key=lambda t: t.path):
        if str(t.file) not in failed:
            continue
        reason = failed[str(t.file)]
        issue = allowed.get(t.md5) or (t.dir_md5 and allowed.get(t.dir_md5))
        if issue:
            print(
                f"ALLOWED: {t.path} (md5 {t.md5}, written by {t.writer}, {issue}): {reason}"
            )
        else:
            bad = True
            print(
                f"UNREADABLE: {t.path} (md5 {t.md5}, written by {t.writer}): {reason}"
            )
    files = len({t.file for t in present})
    others = len(targets) - len(xlc)
    print(
        f"READ: {files} files ({len(present)} outputs) in {seconds:.1f} s. "
        f"Not read: {others} outputs in no xlc format"
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
