#!/usr/bin/env bash
# check-dvc-lock.sh — 渡された dvc.lock を dvc が読めるかを確かめる（.pre-commit-config.yaml の dvc-lock が呼ぶ）
#
# Usage: check-dvc-lock.sh <examples/<問題>/dvc.lock>...
#
# dvc stage list は dvc.lock を読み込んで検証する。先頭の schema: '2.0' が無いと終了コード 255 で止まる。
# dvc status は出力のハッシュを計算するので vera で約 27 秒かかるため使わない。stage list は約 2 秒。
# 読めない lock は、残りも見てから止める（最初の失敗で止めると、直す順が分からなくなる）。stage list の失敗の原因は
# lock の schema 行だけとは限らない（dvc.yaml の vars の参照先、dvc が PATH に無い、など）。

set -euo pipefail

if (( $# == 0 )); then
  echo "ERROR: 確かめる dvc.lock が無い（引数に lock を渡す。空のまま通すと、何も確かめずに成功する）" >&2
  exit 2
fi

status=0
for lock in "$@"; do
  dvc stage list "$(dirname -- "$lock")/dvc.yaml" >/dev/null || {
    echo "FAILED: $lock を dvc stage list が読めない（先頭の schema: '2.0' の行、同じ場所の dvc.yaml、dvc.yaml の vars、dvc の有無を確かめる）" >&2
    status=1
  }
done
exit "$status"
