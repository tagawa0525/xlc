#!/usr/bin/env bash
# check-references-shell.sh — references の devShell の store のパスと、追跡している references-shell.txt を比べる
# （.pre-commit-config.yaml の references-shell が、flake.nix、flake.lock、references-shell.txt を変えるコミットで呼ぶ）
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

shell=$(nix eval --raw .#devShells.x86_64-linux.references.outPath)

stale=()
for record in "${records[@]}"; do
  recorded=$(<"$record")
  if [[ "$recorded" != "$shell" ]]; then
    echo "FAILED: $record は $recorded、references の devShell は $shell" >&2
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
  references-shell.txt をコミットする。参照のステージは changed になり、結果が要るときに再実行する（ADR-0130）。
  dvc commit で lock を付け直さない（ADR-0212 の「(b)」の対象外）
MSG
for record in "${stale[@]}"; do
  echo "    dvc repro $(dirname -- "$record")/dvc.yaml:references-shell" >&2
done
exit 1
