#!/usr/bin/env bash
# M0: ENDF → 多群ライブラリ → ピンセルの2群定数 → 炉心の k_eff を順に計算する。
#
# 使い方: examples/m0/run.sh
# ENDF/B-VIII.0 の中性子サブライブラリは references/endf/ENDF-B-VIII.0_neutrons に置く
# （別の場所なら xs.toml の endf-directory を書き換える）。
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/examples/m0"
bin="$root/target/release"

cargo build --release --quiet --manifest-path "$root/Cargo.toml"

echo "== xs: ENDF/B-VIII.0 -> WIMS 69-group library"
"$bin/xs" "$here/xs.toml"

echo "== lattice: pin cell MOC -> 2-group constants"
"$bin/lattice" "$here/pincell.toml"

echo "== core: 3D NEM"
"$bin/core" "$here/core.inp"
