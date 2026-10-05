#!/usr/bin/env bash
# M0: ENDF → 多群ライブラリ → ピンセルの2群定数 → 炉心の k_eff を順に計算する。
#
# 使い方: examples/m0/run.sh
# ENDF/B-VIII.0 の中性子サブライブラリのディレクトリは環境変数 ENDF で指定する
# （既定は references/endf/ENDF-B-VIII.0_neutrons）。
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
endf="${ENDF:-$root/references/endf/ENDF-B-VIII.0_neutrons}"
here="$root/examples/m0"
bin="$root/target/release"

cargo build --release --quiet --manifest-path "$root/Cargo.toml"

echo "== xs: ENDF/B-VIII.0 -> WIMS 69-group library"
"$bin/xs" "$here/wims69.lib" \
  H1="$endf/n-001_H_001.endf" \
  O16="$endf/n-008_O_016.endf" \
  U235="$endf/n-092_U_235.endf" \
  U238="$endf/n-092_U_238.endf"

echo "== lattice: pin cell MOC -> 2-group constants"
"$bin/lattice" "$here/pincell.inp"

echo "== core: 3D NEM"
"$bin/core" "$here/core.inp"
