#!/usr/bin/env bash
# 区間に分けた燃焼が分けない燃焼と統計の範囲で同じ結果を出すことの確認（ADR-0174）を、種 1〜8 で流して compare_split.py で比べる。
# リポジトリの直下から、references の devShell の中で走らせる（OpenMC と examples/depletion/openmc_depletion.py、ライブラリ、chain が要る）。
# 使い方: nix develop .#references -c examples/depletion/aws-check/split_check.sh <作業ディレクトリ>
# 1 本は約 2 分（8 スレッド）。合計 16 本。重い計算ではないが、排他の測定が走っている間は流さない。
set -euo pipefail
work=$(realpath "${1:?作業ディレクトリを渡す}")
dir=$(realpath "$(dirname "$0")")
mkdir -p "$work"
cd "$work"
whole=() split=()
for seed in 1 2 3 4 5 6 7 8; do
  sed "s/^seed = .*/seed = $seed/; s#\"\\.\\./#\"$dir/../#; s/^threads = .*/threads = 8/" "$dir/medium.toml" > "m$seed.toml"
  OMP_NUM_THREADS=8 python3 "$dir/../openmc_depletion.py" run "m$seed.toml" "w$seed" > "w$seed.log" 2>&1
  python3 "$dir/../openmc_depletion.py" extract "m$seed.toml" "w$seed" "w$seed.msgpack" > /dev/null
  OMP_NUM_THREADS=8 python3 "$dir/../openmc_depletion_resume.py" run "m$seed.toml" "s$seed" > "s$seed.log" 2>&1
  python3 "$dir/../openmc_depletion_resume.py" extract "m$seed.toml" "s$seed" "so$seed" > /dev/null
  whole+=("w$seed.msgpack") split+=("so$seed/m$seed.msgpack")
done
python3 "$dir/compare_split.py" "${whole[@]}" -- "${split[@]}"
