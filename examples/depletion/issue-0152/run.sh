#!/usr/bin/env bash
# #0152 の切り分け: 0 MWd/kg の MOX-1 と 1C の格子の k を設定ごとに解き（bin/lattice）、1C の形状のピンセルで燃料の吸収と物質の移行を
# OpenMC の値に置き換える（swap_absorption）。リポジトリの直下で、`cargo build --release -p xlc-lattice --example swap_absorption` の後に実行する。
#
# 使い方: run.sh lattice <入力の名前>...   （examples/depletion/issue-0152/<名前>.toml を解き、logs/<名前>.{out,log} に書く）
#         run.sh swap                      （logs/swap-{MOX1,1C-0}.{txt,err} に書く）
#
# 置き換えの OpenMC は、MOX-1 が Slurm 558（references/openmc/pin-transfer-1C/MOX1-fine）、1C の 0 MWd/kg が Slurm 534
# （references/openmc/pin-transfer-1C/1C-0）の openmc_pin_transfer.py の出力
set -euo pipefail
here=examples/depletion/issue-0152
mkdir -p "$here/logs"
case $1 in
  lattice)
    shift
    for name in "$@"; do
      t=$SECONDS
      bin/lattice "$here/$name.toml" >"$here/logs/$name.out" 2>"$here/logs/$name.log"
      echo "elapsed $((SECONDS - t)) s" >>"$here/logs/$name.log"
    done
    ;;
  swap)
    t=$SECONDS
    target/release/examples/swap_absorption "$here/pin-MOX1-ring.toml" \
      references/openmc/pin-transfer-1C/MOX1-fine/transfer.msgpack "$here/swaps.txt" \
      --transfer references/openmc/pin-transfer-1C/MOX1-fine/transfer.msgpack \
      >"$here/logs/swap-MOX1.txt" 2>"$here/logs/swap-MOX1.err"
    echo "elapsed $((SECONDS - t)) s" >>"$here/logs/swap-MOX1.err"
    t=$SECONDS
    target/release/examples/swap_absorption examples/depletion/issue-0133/ufg-transfer/pin-1C-0-ring.toml \
      references/openmc/pin-transfer-1C/1C-0/transfer.msgpack "$here/swaps-1c.txt" \
      --transfer references/openmc/pin-transfer-1C/1C-0/transfer.msgpack \
      >"$here/logs/swap-1C-0.txt" 2>"$here/logs/swap-1C-0.err"
    echo "elapsed $((SECONDS - t)) s" >>"$here/logs/swap-1C-0.err"
    ;;
  *) echo "unknown command: $1" >&2; exit 2 ;;
esac
