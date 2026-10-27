#!/usr/bin/env bash
# ADR-0207 の測定のジョブ: 4 つの状態の run.sh を同時に流し、そのあと 1C-0 と 2A-40 を P0 の散乱だけで解いた k を logs/p0.txt に書く。
# リポジトリの直下で、`cargo build --release -p xlc-lattice --example swap_transfer --example balance --example fine_flux` の後に実行する。
# 状態の OpenMC の行と k は 1C-0、1C-5、2A-40 が Slurm 534、MOX1 が 558（細かい入射の実行で、核も同じ実行）。核は 2A-40 が 557。
set -euo pipefail
R=references/openmc/pin-transfer-1C
K="2A-40=$R/2A-40-fine-seed2/transfer.msgpack,MOX1=$R/MOX1-fine/transfer.msgpack"
S=water:27-27,water:14-26,fuel:27-27,fuel:14-26
W=water:27-27,water:14-26
run=examples/depletion/issue-0133/ufg-transfer/run.sh
pids=()
$run 1C-0 $R/1C-0/transfer.msgpack references/ufg-sab/1C-0-sab-fast3-fine.msgpack "$K" $S $W &
pids+=($!)
$run 1C-5 $R/1C-5/transfer.msgpack references/ufg-transfer/1C-5-fine.msgpack "$K" $S $W &
pids+=($!)
$run 2A-40 $R/2A-40/transfer.msgpack references/ufg-sab/2A-40-sab-fast3-fine.msgpack "2A-40=$R/2A-40-fine-seed2/transfer.msgpack" $S $S &
pids+=($!)
$run MOX1 $R/MOX1-fine/transfer.msgpack references/ufg-sab/MOX1-sab-fast3-fine.msgpack "MOX1=$R/MOX1-fine/transfer.msgpack" $S $S &
pids+=($!)
for p in "${pids[@]}"; do
  wait "$p"
done
p0=examples/depletion/issue-0133/ufg-transfer/logs/p0.txt
: >"$p0"
for s in 1C-0 2A-40; do
  for m in dilute ring; do
    k=$(target/release/lattice "examples/depletion/issue-0133/ufg-transfer/pin-$s-$m-p0.toml" 2>/dev/null)
    echo "$s $m $k" >>"$p0"
  done
done
