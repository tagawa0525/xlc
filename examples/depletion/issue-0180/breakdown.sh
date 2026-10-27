#!/usr/bin/env bash
# #0180 の自己遮蔽の時間の内訳: 1C の燃焼のピン（1C-1step.toml）と 2A の集合体の燃焼（2A-1step.toml）の最初の自己遮蔽を、
# 与えた lattice の実行ファイル（名前=パス の組。例: before=/tmp/lattice-before）ごとに解き、自己遮蔽の区間（ログの
# `shielding` の span）が閉じたところで止める。止めるので燃焼の出力は書かれない。debug のログは作業ディレクトリ <dir> の
# <入力>-<名前>.log（途中のチェックポイントは入力の隣に書かれ、消す）、要約（ピンごとの UFG の内訳の行と自己遮蔽の区間の行）は logs/<入力>-<名前>.txt
#
# Usage: breakdown.sh <dir> <名前>=<lattice> ...（Slurm の --timed のジョブで、排他で実行する）
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
work=$(realpath -m "$1")
shift
mkdir -p "$work"
for input in 1C-1step 2A-1step; do
  for pair in "$@"; do
    name=${pair%%=*}
    bin=$(realpath "${pair#*=}")
    log="$work/$input-$name.log"
    rm -f "$here/$input.checkpoint"
    (cd "$here" && RUST_LOG=xlc_lattice=debug exec "$bin" "$input.toml" > /dev/null 2> "$log") &
    pid=$!
    until grep -q "INFO shielding: xlc_lattice::shielding: close" "$log" 2>/dev/null; do
      sleep 1
      kill -0 "$pid" 2>/dev/null || { echo "ERROR: $input $name ended before the first shielding" >&2; exit 1; }
    done
    kill "$pid"
    wait "$pid" || true
    rm -f "$here/$input.checkpoint"
    grep -E "UFG slowing down|UFG free-gas sweep|INFO shielding: xlc_lattice::shielding: (close|fuel pins)" "$log" > "$here/logs/$input-$name.txt"
    echo "$input $name: $(grep 'INFO shielding: xlc_lattice::shielding: close' "$log")"
  done
done
