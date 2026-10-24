#!/usr/bin/env bash
# ADR-0204 の測定: 1 つの状態の UFG の散乱核ごとに、swap_absorption（吸収の置き換え）と fine_flux（格子点ごとの中性子束）を実行する。
# リポジトリの直下で、`cargo build --release -p xlc-lattice --example swap_absorption --example fine_flux` の後に実行する。
#
# 使い方: run.sh <状態> <OpenMC の openmc-rates.msgpack> <置き換えの一覧> <swap_absorption の散乱核,...> <fine_flux の散乱核,...>
#
# 散乱核の入力: static は ../pin-<状態>.toml、light、light-cap、light-fast は ADR-0203 の ../ufg-free-gas/pin-<状態>-{exact,cap400,fast}.toml、
# ほかは pin-<状態>-<散乱核>.toml。ログは logs/<状態>-<散乱核>-{absorption,fine}.{txt,err}、細かい中性子束は
# references/ufg-heavy-free-gas/<状態>-<散乱核>-fine.msgpack。
set -euo pipefail
state=$1 rates=$2 swaps=$3 absorption=$4 fine=$5
here=examples/depletion/issue-0133/ufg-heavy-free-gas
out=references/ufg-heavy-free-gas
mkdir -p "$here/logs" "$out"
input() {
  case $1 in
    static) echo "$here/../pin-$state.toml" ;;
    light) echo "$here/../ufg-free-gas/pin-$state-exact.toml" ;;
    light-cap) echo "$here/../ufg-free-gas/pin-$state-cap400.toml" ;;
    light-fast) echo "$here/../ufg-free-gas/pin-$state-fast.toml" ;;
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
IFS=, read -ra modes <<<"$fine"
for m in "${modes[@]}"; do
  log=$here/logs/$state-$m-fine
  t=$SECONDS
  target/release/examples/fine_flux "$(input "$m")" "$out/$state-$m-fine.msgpack" >"$log.err" 2>&1
  echo "elapsed $((SECONDS - t)) s" >>"$log.err"
done
