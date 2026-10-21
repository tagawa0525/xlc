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
# 終了コード: 戻せない出力があれば 0 でない。dvc checkout は、dvc.lock にまだ hash が無いステージ（No file hash
# info）と、キャッシュに無い出力（OpenMC の結果など、作った機械から持ってきていないもの）を戻せず、255 で終わる。
# 戻せなかったものは dvc の出力に並ぶので、作業に要るものが入っていなければ dvc repro で作る。黙って捨てない。

set -uo pipefail

cd "$(git rev-parse --show-toplevel)" || exit

# references が main への symlink でないと、dvc は worktree の中に空のキャッシュを黙って作る（ADR-0126）
if [[ ! -L references ]]; then
  echo "ERROR: references が main への symlink でない。worktree-add.sh の SHARED: references を確かめる" >&2
  exit 1
fi

status=0

log=$(mktemp)
trap 'rm -f "$log"' EXIT
dvc checkout 2>&1 | tee "$log"
rc=${PIPESTATUS[0]}
if (( rc != 0 )); then
  unlocked=$(grep -c 'No file hash info found' "$log" || true)
  echo "WARNING: dvc checkout が終了コード $rc で終わった。戻っていない出力がある" \
    "（dvc.lock に hash が無いもの $unlocked 行。ほかはキャッシュに無い出力など。上の一覧を見る）" >&2
  status=$rc
fi

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
