#!/usr/bin/env bash
# flake の NJOY と OpenMC が、references/ に手でビルドした NJOY と OpenMC と同じ結果を出すかを比べる（ADR-0124 の確認）。
#
# 使い方（worktree の直下で）: trials/pin-reference-assets/compare_tools.sh <main の checkout> <作業ディレクトリ>
#   old: 手でビルドしたツール（<main>/references/NJOY2016/build/njoy、<main>/references/openmc/{install/bin,venv}）を、
#        main の flake の references の devShell で動かす（OpenMC が GCC 13 の libgomp にリンクしているため）
#   new: この worktree の flake の references の devShell のツール
# NJOY は references/njoy-runs の入力（B10-acer は ACER のパッチを通る）、OpenMC は C5G7（多群）と VERA 1B（連続
# エネルギー、references/openmc/endfb80）を、同じ粒子数・乱数のシード・スレッド数で解く。同じツールを 2 回実行した差と
# 比べるため、old と new をそれぞれ 2 回実行する。openmc.data（Python API）が NJOY を通して作る H-1 の HDF5 も比べる。
set -euo pipefail

main=$(realpath "$1")
work=$(realpath -m "$2")
here=$(pwd)
mkdir -p "$work"

old() { nix develop "$main#references" -c env PATH="$main/references/openmc/install/bin:$main/references/NJOY2016/build:$PATH" "$@"; }
new() { nix develop "$here#references" -c "$@"; }
python_of() { if [[ $1 == old ]]; then echo "$main/references/openmc/venv/bin/python"; else echo python3; fi; }

for side in old new; do
  for run in 1 2; do
    for case in B10-acer O17-groupr; do
      d=$work/njoy-$side$run/$case
      mkdir -p "$d"
      cp "$main/references/njoy-runs/$case/input" "$main/references/njoy-runs/$case/tape20" "$d/"
      (cd "$d" && $side sh -c 'njoy < input > log')
    done
    $side sh -c "cd $work && OMP_NUM_THREADS=16 $(python_of $side) \
      $here/examples/c5g7/openmc_model.py c5g7-$side$run 20000 120 20 1 > c5g7-$side$run.log"
    $side sh -c "cd $work && OMP_NUM_THREADS=16 $(python_of $side) \
      $here/examples/vera/openmc_model.py $main/references/openmc/endfb80/cross_sections.xml 1B vera-1B-$side$run \
      20000 120 20 1 > vera-1B-$side$run.log"
    $side sh -c "cd $work && $(python_of $side) $here/trials/pin-reference-assets/openmc_data_h1.py \
      $main/references/endf/ENDF-B-VIII.0_neutrons/n-001_H_001.endf h1-$side$run.h5 > h1-$side$run.log"
  done
done

# 違うバイトの位置（1 始まり）を「最初-最後 (個数)」で書く。cmp は違うと 1 を返す
bytes() {
  local listing
  listing=$(cmp -l "$1" "$2" || [[ $? == 1 ]])
  awk 'NF { if (!n) f = $1; l = $1; n++ } END { print n ? f "-" l " (" n " bytes)" : "identical" }' <<<"$listing"
}

echo "== NJOY: tapes (old1 vs new1; old1 vs old2)"
for case in B10-acer O17-groupr; do
  for t in "$work/njoy-old1/$case"/tape*; do
    n=$(basename "$t")
    echo "$case $n: $(bytes "$t" "$work/njoy-new1/$case/$n"); $(bytes "$t" "$work/njoy-old2/$case/$n")"
  done
done
echo "== openmc.data + NJOY: H-1 at 600 K (old1 vs new1; old1 vs old2)"
nix develop "$here#default" -c python3 "$here/trials/pin-reference-assets/diff_hdf5.py" "$work/h1-old1.h5" "$work/h1-new1.h5"
nix develop "$here#default" -c python3 "$here/trials/pin-reference-assets/diff_hdf5.py" "$work/h1-old1.h5" "$work/h1-old2.h5"
echo "== OpenMC: version"
old openmc --version | head -3
new openmc --version | head -3
echo "== OpenMC: k"
for p in c5g7 vera-1B; do
  for r in old1 old2 new1 new2; do echo "$p $r: $(tail -1 "$work/$p-$r.log")"; done
  for pair in "old1 new1" "old1 old2" "new1 new2"; do
    read -r a b <<<"$pair"
    echo "$p $a vs $b:"
    nix develop "$here#default" -c python3 "$here/trials/pin-reference-assets/diff_statepoints.py" \
      "$work/$p-$a/statepoint.120.h5" "$work/$p-$b/statepoint.120.h5"
  done
done
