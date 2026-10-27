#!/usr/bin/env bash
# #0152 の測定: 1 つの状態で、物質の群の移行の行を、OpenMC の細かい入射の区間の核を xlc の群の中の中性子束（`xlc`）か 1/E（`lethargy`）で
# まとめた行に置き換えた swap_transfer（group_rows.py write）と、OpenMC の状態の行そのものに置き換えた swap_transfer を実行する。
# リポジトリの直下で、`cargo build --release -p xlc-lattice --example swap_transfer --example balance` の後に実行する。
#
# 使い方: run.sh <状態> <xlc の入力> <状態の transfer.msgpack> <xlc の fine.msgpack> <核の名前>=<核の transfer.msgpack>,... <置き換え>,...
#
# 置き換えは `<物質>:<最初の群>-<最後の群>`。ログは logs/<状態>-transfer.txt（OpenMC の行。置き換えごとに `material <物質> <群>` の行を
# 前に置く）、logs/<状態>-<核の名前>-<重み>-transfer.txt（核の行）と logs/<状態>-rows.txt（group_rows.py rows）。合成した transfer.msgpack は
# references/group-rows/ に置く（git で追跡しない）。
set -euo pipefail
state=$1 input=$2 transfer=$3 fine=$4 kernels=$5 specs=$6
here=examples/depletion/issue-0133/group-rows
out=references/group-rows
mkdir -p "$here/logs" "$out"
target/release/examples/balance "$input" "$out/$state-balance.msgpack" >"$here/logs/$state-balance.txt"
IFS=, read -ra kernel_list <<<"$kernels"
IFS=, read -ra spec_list <<<"$specs"
swap() { # <transfer.msgpack> <ログ>
  : >"$2.txt"
  : >"$2.err"
  for spec in "${spec_list[@]}"; do
    material=${spec%%:*} groups=${spec#*:}
    echo "material $material $groups" >>"$2.txt"
    target/release/examples/swap_transfer "$input" "$1" "$material" "${groups%-*}" "${groups#*-}" >>"$2.txt" 2>>"$2.err"
  done
}
swap "$transfer" "$here/logs/$state-transfer"
: >"$here/logs/$state-rows.txt"
for entry in "${kernel_list[@]}"; do
  name=${entry%%=*} kernel=${entry#*=}
  for spec in "${spec_list[@]}"; do
    material=${spec%%:*} groups=${spec#*:}
    echo "kernel $name material $material $groups" >>"$here/logs/$state-rows.txt"
    python3 examples/depletion/issue-0133/group_rows.py rows "$transfer" "$kernel" "$out/$state-balance.msgpack" "$fine" \
      "$material" "${groups%-*}" "${groups#*-}" >>"$here/logs/$state-rows.txt"
  done
  for weight in xlc lethargy; do
    synthetic=$out/$state-$name-$weight.msgpack
    cp "$transfer" "$synthetic"
    for spec in "${spec_list[@]}"; do
      material=${spec%%:*} groups=${spec#*:}
      python3 examples/depletion/issue-0133/group_rows.py write "$synthetic" "$kernel" "$fine" "$material" \
        "${groups%-*}" "${groups#*-}" "$weight" "$synthetic"
    done
    swap "$synthetic" "$here/logs/$state-$name-$weight-transfer"
  done
done
