#!/usr/bin/env bash
# 2A の 1 輸送を、燃料（新しい／60 MWd/kg）× ufg-tolerance（なし／1e-4）× スレッド数 1、4、16 で実行し、ログを $LOGS に書く。
# スレッドは物理コアの 0 番から n 個に固定する（r995 は CPU 0〜15 が別々の物理コア。16〜31 はその SMT の相方）。
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
bin=$here/../../../target/release/lattice
LOGS=${LOGS:?}
mkdir -p "$LOGS"
for fuel in fresh burned; do
  for tol in none 1e-4; do
    for n in 1 4 16; do
      name=$fuel-$tol-t$n
      start=$(date +%s.%N)
      RUST_LOG=debug RAYON_NUM_THREADS=$n taskset -c 0-$((n - 1)) "$bin" "$here/lattice-2A-$fuel-$tol.toml" 2> "$LOGS/$name.log" > "$LOGS/$name.out"
      end=$(date +%s.%N)
      echo "$name wall $(echo "$end - $start" | bc)" | tee -a "$LOGS/wall.txt"
    done
  done
done
