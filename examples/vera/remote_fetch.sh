#!/usr/bin/env bash
# VM で流した参照解の結果を手元に戻す（docs/environment.md の「外部の計算機での参照解」）。
#
# 使い方: examples/vera/remote_fetch.sh [--dry-run] <ssh の接続先> <VM の作業ディレクトリ> <手元の置き場>
#
# results/<入力の名前>/ の k.json、pins.msgpack、layers.msgpack などと、最後の世代の statepoint（extract をやり直すため）、
# log、done を戻す。途中の statepoint（再開用）は戻さない。手元で extract し直すには:
#   python3 examples/vera/openmc_model_3d.py extract <入力> <手元の置き場>/<名前>/run <結果のディレクトリ>
set -euo pipefail
flags=(-a --partial --info=stats1 --prune-empty-dirs)
if [[ ${1:-} == --dry-run ]]; then flags+=(--dry-run --itemize-changes); shift; fi
if [[ $# -ne 3 || $1 == --help || $1 == -h ]]; then
  sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
# 最後の世代の statepoint は、世代数が入力ごとに違うので、最大の番号のものを VM で選ぶ
rsync "${flags[@]}" --exclude='run/statepoint.*.h5' --exclude='run/summary.h5' "$1:$2/results/" "$3/"
for name in $(ssh "$1" "ls $2/results"); do
  last=$(ssh "$1" "ls $2/results/$name/run | grep -E '^statepoint\\.[0-9]+\\.h5$' | sort -t. -k2 -n | tail -1" || true)
  [[ -n $last ]] && rsync "${flags[@]}" "$1:$2/results/$name/run/$last" "$3/$name/run/"
done
