#!/usr/bin/env bash
# ADR-0205 の測定: 1 つの状態の UFG の散乱核ごとに、swap_absorption（吸収の置き換え）、swap_transfer（移行の置き換え）、fine_flux
# （格子点ごとの中性子束）を実行する。リポジトリの直下で、
# `cargo build --release -p xlc-lattice --example swap_absorption --example swap_transfer --example fine_flux` の後に実行する。
#
# 使い方: run.sh <状態> <OpenMC の openmc-rates.msgpack> <置き換えの一覧> <OpenMC の transfer.msgpack>
#             <swap_absorption の散乱核,...> <swap_transfer の散乱核,...> <fine_flux の散乱核,...>
#
# 散乱核の入力: static は ../pin-<状態>.toml、heavy、heavy-cap、heavy-fast3 は ADR-0204 の ../ufg-heavy-free-gas/pin-<状態>-<散乱核>.toml、
# ほかは pin-<状態>-<散乱核>.toml。移行の置き換えは、減速材と燃料の 27 群（4〜9.877 eV）と 14〜26 群（9.877 eV〜15 keV）。
# ログは logs/<状態>-<散乱核>-{absorption,transfer,fine}.{txt,err}（transfer.txt は置き換えごとに `material <物質>` の行を前に
# 置く）、細かい中性子束は references/ufg-sab/<状態>-<散乱核>-fine.msgpack。
set -euo pipefail
state=$1 rates=$2 swaps=$3 transfer=$4 absorption=$5 transfers=$6 fine=$7
here=examples/depletion/issue-0133/ufg-sab
out=references/ufg-sab
mkdir -p "$here/logs" "$out"
input() {
  case $1 in
    static) echo "$here/../pin-$state.toml" ;;
    heavy | heavy-cap | heavy-fast3) echo "$here/../ufg-heavy-free-gas/pin-$state-$1.toml" ;;
    *) echo "$here/pin-$state-$1.toml" ;;
  esac
}
IFS=, read -ra modes <<<"$absorption"
for m in "${modes[@]}"; do
  log=$here/logs/$state-$m-absorption
  t=$SECONDS
  target/release/examples/swap_absorption "$(input "$m")" "$rates" "$swaps" >"$log.txt" 2>"$log.err"
  echo "elapsed $((SECONDS - t)) s" >>"$log.err"
done
IFS=, read -ra modes <<<"$transfers"
for m in "${modes[@]}"; do
  log=$here/logs/$state-$m-transfer
  : >"$log.txt"
  : >"$log.err"
  t=$SECONDS
  for spec in "water 27 27" "water 14 26" "fuel 27 27" "fuel 14 26"; do
    read -r material first last <<<"$spec"
    echo "material $material" >>"$log.txt"
    target/release/examples/swap_transfer "$(input "$m")" "$transfer" "$material" "$first" "$last" >>"$log.txt" 2>>"$log.err"
  done
  echo "elapsed $((SECONDS - t)) s" >>"$log.err"
done
IFS=, read -ra modes <<<"$fine"
for m in "${modes[@]}"; do
  log=$here/logs/$state-$m-fine
  t=$SECONDS
  target/release/examples/fine_flux "$(input "$m")" "$out/$state-$m-fine.msgpack" >"$log.err" 2>&1
  echo "elapsed $((SECONDS - t)) s" >>"$log.err"
done
