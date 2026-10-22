#!/usr/bin/env bash
# #0133 の試行: 燃料の核種の吸収と核分裂の断面積を、核種の組とエネルギーの範囲の組ごとに OpenMC の値に置き換えて（`swap_absorption` の例）
# k の変化を測る。k の差への寄与はその符号を変えたもの。
#
# 使い方（references の devShell の外で。リポジトリの直下で `cargo build --release -p xlc-lattice --example swap_absorption` のあと）:
#   swap_categories.sh <作業ディレクトリ（lattice.toml、openmc-rates.msgpack）>
#
# 群は WIMS 69 群の番号: 1〜13 は 15.03 keV より上、14〜27 は UFG の範囲（4 eV〜15.03 keV）、28〜45 は 0.625〜4 eV、46〜56 は 0.1〜0.625 eV、
# 57〜69 は 0.1 eV より下。高速の群は非弾性散乱と (n,2n) を置き換えないので、重い核種だけにする（O-16 の (n,α) は OpenMC の absorption に入る）。
# 「ほかの核種」は U-235、U-238、Pu-239、Pu-241、Xe-135、Sm-149、O-16 を除く燃料の核種。
set -euo pipefail
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
d=$1
others=$(python3 -c "
import msgpack
o = msgpack.unpackb(open('$d/openmc-rates.msgpack', 'rb').read())
main = {'U-235', 'U-238', 'Pu-239', 'Pu-241', 'Xe-135', 'Sm-149', 'O-16'}
print(' '.join(n for n in o['rates'] if n not in main))")
join() {
  local r=$1
  shift
  local out=()
  for n in "$@"; do out+=("$n:$r"); done
  (
    IFS=,
    echo "${out[*]}"
  )
}
fast=U-235:1-13
for n in U-234 U-236 Pu-239 Pu-240 Pu-241 Pu-242; do
  if python3 -c "import msgpack,sys; sys.exit('$n' not in msgpack.unpackb(open('$d/openmc-rates.msgpack','rb').read())['rates'])"; then
    fast+=",$n:1-13"
  fi
done
specs=("U-238:1-13" "$fast" "U-238:14-27" "U-238:28-45" "U-238:46-69" "U-235:14-27" "U-235:28-45" "U-235:46-69")
if [[ " $others " == *" Pu-240 "* ]]; then
  specs+=("Pu-239:14-27,Pu-241:14-27" "Pu-239:28-45,Pu-241:28-45" "Pu-239:46-56,Pu-241:46-56" "Pu-239:57-69,Pu-241:57-69"
    "Xe-135:14-69,Sm-149:14-69")
fi
specs+=("$(join 14-27 $others)" "$(join 28-45 $others)" "$(join 46-69 $others)" "*:14-69,$fast,U-238:1-13")
cd "$d"
"$root/target/release/examples/swap_absorption" lattice.toml openmc-rates.msgpack "${specs[@]}"
