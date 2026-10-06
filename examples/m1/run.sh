#!/usr/bin/env bash
# M1: ENDF → 温度と σ0 に依存する多群ライブラリ → ピンセルの2群定数を順に計算する。
#
# 使い方: examples/m1/run.sh
# ENDF/B-VIII.0 の中性子と熱中性子散乱のサブライブラリは references/endf/ に置く
# （別の場所なら xs.toml の endf-directory と thermal-directory を書き換える）。
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
here="$root/examples/m1"
bin="$root/target/release"

cargo build --release --quiet --manifest-path "$root/Cargo.toml"

echo "== xs: ENDF/B-VIII.0 -> WIMS 69-group library (temperatures, sigma0)"
"$bin/xs" "$here/xs.toml"

echo "== lattice: pin cell MOC -> 2-group constants"
"$bin/lattice" "$here/pincell.inp"
