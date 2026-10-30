#!/usr/bin/env bash
# worktree-restore-outputs.sh — worktree に git で追跡しない出力（bin/ の実行ファイル、各例の出力、記録が引く生成物）を戻す
#
# Usage: scripts/worktree-restore-outputs.sh（worktree のどこからでもよい。dvc が PATH にあること）
#
# 新しい worktree は DVC の出力が空なので、dvc.lock と .dvc ファイルの hash に合わせてキャッシュ（references/dvc-cache。
# main と共有で reflink。ADR-0126）から戻す。引数なしの dvc checkout はすべての dvc.yaml と、記録が引く生成物を追う
# .dvc ファイル（ADR-0222）を対象にする。キャッシュの置き方は .dvc/config の cache.type（reflink だけ）で、中身をコピーしない。
#
# 終了コード: dvc.lock か .dvc ファイルに hash がある出力を戻せなければ 0 でない。checkout は scripts/dvc-checkout.py で行い、hash の
# まだ無いステージ（未実行。No file hash info）の出力は数えるだけにする（#0184）。キャッシュに無い出力（OpenMC の結果
# など、作った機械から持ってきていないもの）は MISSING に並ぶので、作業に要るものが入っていなければ dvc repro で作る。
# 黙って捨てない。

set -uo pipefail

cd "$(git rev-parse --show-toplevel)" || exit

# references が main への symlink でないと、dvc は worktree の中に空のキャッシュを黙って作る（ADR-0126）
if [[ ! -L references ]]; then
  echo "ERROR: references が main への symlink でない。worktree-add.sh の SHARED: references を確かめる" >&2
  exit 1
fi

exec python3 scripts/dvc-checkout.py
