#!/usr/bin/env bash
# 2A の 1 輸送（60 MWd/kg の組成、ufg-tolerance = 1e-4。燃焼の本番と同じ入力）を、実行ファイル <bin> とスレッド数 <n>... で
# 順に実行し、ログを $LOGS/<名前>-t<n>.log に、少数群定数の出力を $LOGS/<名前>-t<n>.msgpack に書く（実行ファイルどうしの
# 出力のビットの一致を比べるため）。プロセスは `perf stat`（ユーザー空間のサイクル、命令、DRAM からの L2 の充填、バックエンドで止まったディスパッチのスロット）の下で実行し、その
# 数を $LOGS/<名前>-t<n>.perf に書く。
# MOC の sweep などの区間は trace なので、その区間だけ有効にする。スレッドは物理コアの 0 番から n 個に固定する
# （r995 は CPU 0〜15 が別々の物理コア。16〜31 はその SMT の相方）。
#
# 使い方: LOGS=<ディレクトリ> run.sh <名前> <bin> <n>...
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
LOGS=$(realpath -m "${LOGS:?}")
name=${1:?}
bin=$(realpath "${2:?}")
shift 2
mkdir -p "$LOGS"
export RUST_LOG='debug,[sweep]=trace,[tallies]=trace,[homogenize]=trace,[coarse]=trace'
input=$LOGS/$name.toml
for n in "$@"; do
  run=$name-t$n
  # 出力の置き場所だけを書き換えた入力（入力の中の相対パスは入力のディレクトリから読むので、LOGS は絶対パスにしてある）
  sed "s|output = \"[^\"]*\"|output = \"$LOGS/$run.msgpack\"|" "$here/lattice-2A-burned-1e-4.toml" > "$input"
  start=$(date +%s.%N)
  RAYON_NUM_THREADS=$n perf stat -x, -o "$LOGS/$run.perf" -e cycles:u,instructions:u,l2_fill_rsp_src.dram_io_near:u,de_no_dispatch_per_slot.backend_stalls:u \
    taskset -c 0-$((n - 1)) "$bin" "$input" 2> "$LOGS/$run.log" > "$LOGS/$run.out"
  end=$(date +%s.%N)
  echo "$run wall $(echo "$end - $start" | bc)" | tee -a "$LOGS/wall.txt"
done
