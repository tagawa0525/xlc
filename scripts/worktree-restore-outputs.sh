#!/usr/bin/env bash
# worktree-restore-outputs.sh — worktree に git で追跡しない出力（bin/ の実行ファイル、各例の出力、参照解）を戻す
#
# Usage: scripts/worktree-restore-outputs.sh（worktree のどこからでもよい。dvc が PATH にあること）
#
# 新しい worktree は DVC の出力が空なので、dvc.lock の hash に合わせてキャッシュ（references/dvc-cache。main と共有で
# reflink。ADR-0126）から戻す。引数なしの dvc checkout はすべての dvc.yaml を対象にする。手元（btrfs、reflink）で
# 約 10 秒、追加の使用量は約 14 MB（2026-10-21、379e81d6a）。reflink が効かない場所では全体のコピーになる。
# 続けて、DVC の出力でない参照解 examples/vera/p2/openmc-2K.msgpack（ADR-0121）を main からコピーする。
#
# 終了コード: dvc.lock に hash がある出力を戻せなければ 0 でない。checkout は scripts/dvc-checkout.py で行い、hash の
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

status=0

python3 scripts/dvc-checkout.py || status=$?

main_root=$(git rev-parse --path-format=absolute --git-common-dir)/..
pin=examples/vera/p2/openmc-2K.msgpack
if [[ -e "$pin" ]]; then
  :
elif [[ -e "$main_root/$pin" ]]; then
  cp --reflink=auto "$main_root/$pin" "$pin"
else
  echo "ERROR: $pin が main にも無い。参照解を作る（スキル openmc-reference）か、持っている worktree からコピーする" >&2
  status=1
fi

exit "$status"
