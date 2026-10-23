#!/usr/bin/env bash
# #0163 の直列の部分の測定: 2A の 1 輸送（run.sh と同じ入力）を、実行ファイル $BIN（既定は target/release/lattice）、スレッド数
# $THREADS（既定は 1、4、16）で順に実行し、ログを $LOGS/<名前>-t<n>.log に、少数群定数の出力を $LOGS/<名前>-t<n>.msgpack に書く
# （実行ファイルどうしの出力のビットの一致を比べるため）。<名前> は $NAME（既定は burned-1e-4）。MOC の反復の段と輸送の外の段の区間
# （trace と debug）を有効にする。スレッドは物理コアの 0 番から n 個に固定する（r995 は CPU 0〜15 が別々の物理コア）。
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
bin=$(realpath "${BIN:-$here/../../../target/release/lattice}")
LOGS=$(realpath -m "${LOGS:?}")
name=${NAME:-burned-1e-4}
mkdir -p "$LOGS"
spans=(setup equilibrium sweep tallies update homogenize coarse scale converge finish)
RUST_LOG=debug$(printf ',[%s]=trace' "${spans[@]}")
export RUST_LOG
for n in ${THREADS:-1 4 16}; do
  run=$name-t$n
  # 出力の置き場所だけを書き換えた入力（入力の中のパスは絶対パス）
  input=$LOGS/$run.toml
  sed "s|output = \"[^\"]*\"|output = \"$LOGS/$run.msgpack\"|" "$here/lattice-2A-burned-1e-4.toml" > "$input"
  start=$(date +%s.%N)
  RAYON_NUM_THREADS=$n taskset -c 0-$((n - 1)) "$bin" "$input" 2> "$LOGS/$run.log" > "$LOGS/$run.out"
  end=$(date +%s.%N)
  echo "$run wall $(echo "$end - $start" | bc)" | tee -a "$LOGS/wall.txt"
done
