#!/usr/bin/env bash
# ADR-0207 の測定: 1 つの状態で、UFG の範囲の群の散乱行列の移行の形を無限希釈（`dilute`。../ufg-sab/pin-<状態>-sab-fast3.toml）と
# リングの中性子束で数えた移行（`ring`。pin-<状態>-ring.toml）にした xlc の、OpenMC の状態の行への swap_transfer（SPH 因子は置き換える
# 前のまま）と、物質の行（balance の例）を OpenMC の行と核の行に比べた group_rows.py rows を実行する。リポジトリの直下で、
# `cargo build --release -p xlc-lattice --example swap_transfer --example balance --example fine_flux` の後に実行する。
#
# 使い方: run.sh <状態> <状態の transfer.msgpack> <xlc の fine.msgpack> <核の名前>=<核の transfer.msgpack>,... <置き換え>,... <rows の置き換え>,... [<形>,...]
#
# 形は省くと dilute と ring。ADR-0208 の測定は、それぞれの格子を間引いた `dilute-thin` と `ring-thin`（pin-<状態>-<形>.toml）を選ぶ。
#
# 置き換えは `<物質>:<最初の群>-<最後の群>`。rows の置き換えは、核の行と比べる物質と群（燃料の核は同じ状態の実行にしかないので、
# 状態ごとに選ぶ。空なら rows を実行しない）。<xlc の fine.msgpack> が無ければ、無限希釈の入力の fine_flux で作る（リングの中性子束は
# dilute と ring で同じ。間引いた形でも、核は間引かない格子のこの中性子束でまとめる）。ログは logs/<状態>-<形>-transfer.txt（置き換えごとに `material <物質> <群>` の行を前に置く）と
# logs/<状態>-<形>-rows.txt（核と置き換えごとに `kernel <名前> material <物質> <群>` の行を前に置く）、balance の出力は
# references/ufg-transfer/<状態>-<形>-balance.msgpack（git で追跡しない）。
set -euo pipefail
state=$1 transfer=$2 fine=$3 kernels=$4 specs=$5 row_specs=$6 modes=${7:-dilute,ring}
here=examples/depletion/issue-0133/ufg-transfer
out=references/ufg-transfer
mkdir -p "$here/logs" "$out"
input() {
  case $1 in
    dilute) echo "$here/../ufg-sab/pin-$state-sab-fast3.toml" ;;
    ring | dilute-thin | ring-thin) echo "$here/pin-$state-$1.toml" ;;
    *)
      echo "unknown mode $1" >&2
      exit 1
      ;;
  esac
}
if [[ ! -e $fine ]]; then
  target/release/examples/fine_flux "$(input dilute)" "$fine" >"$here/logs/$state-fine.err" 2>&1
fi
IFS=, read -ra spec_list <<<"$specs"
IFS=, read -ra row_list <<<"$row_specs"
IFS=, read -ra kernel_list <<<"$kernels"
IFS=, read -ra mode_list <<<"$modes"
for mode in "${mode_list[@]}"; do
  input "$mode" >/dev/null
done
for mode in "${mode_list[@]}"; do
  log=$here/logs/$state-$mode
  : >"$log-transfer.txt"
  : >"$log-transfer.err"
  for spec in "${spec_list[@]}"; do
    material=${spec%%:*} groups=${spec#*:}
    echo "material $material $groups" >>"$log-transfer.txt"
    target/release/examples/swap_transfer "$(input "$mode")" "$transfer" "$material" "${groups%-*}" "${groups#*-}" \
      >>"$log-transfer.txt" 2>>"$log-transfer.err"
  done
  balance=$out/$state-$mode-balance.msgpack
  target/release/examples/balance "$(input "$mode")" "$balance" >"$log-balance.txt"
  : >"$log-rows.txt"
  for entry in "${kernel_list[@]}"; do
    name=${entry%%=*} kernel=${entry#*=}
    for spec in "${row_list[@]}"; do
      material=${spec%%:*} groups=${spec#*:}
      echo "kernel $name material $material $groups" >>"$log-rows.txt"
      python3 examples/depletion/issue-0133/group_rows.py rows "$transfer" "$kernel" "$balance" "$fine" \
        "$material" "${groups%-*}" "${groups#*-}" >>"$log-rows.txt"
    done
  done
done
