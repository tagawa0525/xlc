#!/usr/bin/env bash
# check-references-shell.sh — references の devShell の store のパスと、追跡している references-shell.txt を比べる
# （.pre-commit-config.yaml の references-shell が呼ぶ）
#
# Usage: scripts/check-references-shell.sh（リポジトリのどこからでもよい）
#
# references-shell.txt は NJOY と OpenMC のステージの deps で、git で追跡する（cache: false。ADR-0130）。書き出す
# references-shell のステージは always_changed で、dvc repro が走るまで書き直されないので、flake.nix の変更で devShell の
# store のパスが変わっても dvc status は参照解を changed と示さない。食い違ったまま main に入ると、次の dvc repro が
# 参照解（数時間）までを再実行する（#0183）。比べるのは git が追跡するファイルだけ（追跡しない試しのファイルで止めない）。
# 食い違いは、すべてを挙げてから止める（最初の 1 つで止めると、直す範囲が分からない）。

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

mapfile -t records < <(git ls-files -- ':(glob)**/references-shell.txt')
if (( ${#records[@]} == 0 )); then
  echo "ERROR: 追跡している references-shell.txt が無い（比べるものが無いまま通すと、何も確かめずに成功する）" >&2
  exit 2
fi

# 評価で flake.lock を書き換えない（コミットに入らない lock で比べることになる）。入力が lock に無ければ止まる
shell=$(nix eval --no-update-lock-file --raw .#devShells.x86_64-linux.references.outPath)

stale=()
for record in "${records[@]}"; do
  if [[ ! -f "$record" ]]; then
    echo "FAILED: $record は追跡しているが、作業ツリーに無い" >&2
    stale+=("$record")
    continue
  fi
  # DVC はバイトで md5 を取り、nix eval --raw は末尾に改行を書かないので、文字列でなくバイトで比べる（$(<file) は末尾の
  # 改行を落とす）
  if ! cmp -s -- "$record" <(printf '%s' "$shell"); then
    echo "FAILED: $record は $(<"$record")（バイトで比べる。末尾の改行も違いになる）、references の devShell は $shell" >&2
    stale+=("$record")
  fi
done
(( ${#stale[@]} == 0 )) && exit 0

cat >&2 <<MSG

references の devShell の store のパスが、NJOY と OpenMC のステージの deps（上の references-shell.txt）と食い違う。
このまま入ると、dvc status に出ないまま、次の dvc repro が参照解（数時間）までを再実行する（ADR-0212、#0183）。
修正方法（ADR-0212 の「(a)」）:
- 足した Python のパッケージを参照のステージが使わない（検証用のツールやテストだけが使う）なら、flake.nix の
  referencesPythonPackages でなく pythonPackages に足す。references の devShell が変わらなくなり、この検査が通る
- 参照のステージのツール（NJOY、OpenMC、referencesPython）を変えたなら、default の devShell で次を実行して
  references-shell.txt と、同じ実行が書き換える dvc.lock の references-shell の項をコミットする。参照のステージは changed になり、結果が要るときに再実行する（ADR-0130）。
  dvc commit で lock を付け直さない（ADR-0212 の「(b)」の対象外）
MSG
for record in "${stale[@]}"; do
  echo "    dvc repro $(dirname -- "$record")/dvc.yaml:references-shell" >&2
done
exit 1
