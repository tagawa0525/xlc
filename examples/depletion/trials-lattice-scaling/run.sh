#!/usr/bin/env bash
# 2A の 1 輸送（60 MWd/kg の組成、ufg-tolerance = 1e-4。燃焼の本番と同じ入力）をスレッド数 1、4、16 で実行し、ログを $LOGS に書く。
# MOC の sweep などの区間は trace なので、その区間だけ有効にする。スレッドは物理コアの 0 番から n 個に固定する
# （r995 は CPU 0〜15 が別々の物理コア。16〜31 はその SMT の相方）。
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
bin=$here/../../../target/release/lattice
LOGS=${LOGS:?}
mkdir -p "$LOGS"
export RUST_LOG='debug,[sweep]=trace,[tallies]=trace,[homogenize]=trace,[coarse]=trace'
for n in 1 4 16; do
  name=burned-1e-4-t$n
  start=$(date +%s.%N)
  RAYON_NUM_THREADS=$n taskset -c 0-$((n - 1)) "$bin" "$here/lattice-2A-burned-1e-4.toml" 2> "$LOGS/$name.log" > "$LOGS/$name.out"
  end=$(date +%s.%N)
  echo "$name wall $(echo "$end - $start" | bc)" | tee -a "$LOGS/wall.txt"
done
