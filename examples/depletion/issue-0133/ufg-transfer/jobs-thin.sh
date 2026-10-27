#!/usr/bin/env bash
# ADR-0208 の測定のジョブ: 1C-0 と 2A-40 の状態で、ADR-0207 の測定（jobs.sh）の 2 つの移行の形の UFG の格子を許容差 1e-4 で間引いた
# `dilute-thin` と `ring-thin` の run.sh を同時に流し、そのあと 2 つの状態の間引かない `dilute` と `ring` の k を logs/thin-k.txt に書く
# （間引かない形の置き換えと行は jobs.sh のログ。間引かない格子の計算は ADR-0208 で変わらないことを、この k で確かめる）。
# リポジトリの直下で、`cargo build --release -p xlc-lattice --bin lattice --example swap_transfer --example balance --example fine_flux` の後に実行する。
# 状態の OpenMC の行と k、核、xlc の格子点ごとの中性子束は jobs.sh と同じ。
set -euo pipefail
R=references/openmc/pin-transfer-1C
S=water:27-27,water:14-26,fuel:27-27,fuel:14-26
W=water:27-27,water:14-26
here=examples/depletion/issue-0133/ufg-transfer
run=$here/run.sh
pids=()
$run 1C-0 $R/1C-0/transfer.msgpack references/ufg-sab/1C-0-sab-fast3-fine.msgpack \
  "2A-40=$R/2A-40-fine-seed2/transfer.msgpack,MOX1=$R/MOX1-fine/transfer.msgpack" $S $W dilute-thin,ring-thin &
pids+=($!)
$run 2A-40 $R/2A-40/transfer.msgpack references/ufg-sab/2A-40-sab-fast3-fine.msgpack \
  "2A-40=$R/2A-40-fine-seed2/transfer.msgpack" $S $S dilute-thin,ring-thin &
pids+=($!)
for p in "${pids[@]}"; do
  wait "$p"
done
k=$here/logs/thin-k.txt
: >"$k"
for s in 1C-0 2A-40; do
  for m in dilute ring; do
    case $m in
      dilute) i=$here/../ufg-sab/pin-$s-sab-fast3.toml ;;
      ring) i=$here/pin-$s-ring.toml ;;
    esac
    echo "$s $m $(target/release/lattice "$i" 2>"${k%.txt}-$s-$m.err")" >>"$k"
  done
done
