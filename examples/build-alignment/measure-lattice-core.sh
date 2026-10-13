#!/usr/bin/env bash
# ADR-0157 の lattice と core の測定（ジョブ 305）: 6b39ad13 の 4 つの整列の設定の lattice と core を、M0 のピンセルと炉心、
# VERA 1B のピンセル、2A の集合体で 2 回ずつ走らせ、出力と時間を比べる。先に build.sh と measure-xs.sh（M0 の
# ライブラリを使う）。examples/vera/vera.msgpack が要る。Slurm の排他のジョブで走らせる
set -euo pipefail
cd "$(dirname "$0")"
root=$(pwd)
# examples/m0 の DVC の出力（wims69.msgpack、pincell-fg.msgpack）を上書きしないよう、あれば止める
for f in ../m0/wims69.msgpack ../m0/pincell-fg.msgpack; do
  if [ -e "$f" ]; then echo "$f exists; run this in a worktree without the M0 outputs" >&2; exit 1; fi
done
c=6b39ad13
rm -rf out/lc; mkdir -p out/lc
t() { # <out> <tag> <dir> <command...>
  local out=$1 tag=$2 dir=$3; shift 3
  local s e
  s=$(date +%s.%N)
  (cd "$dir" && "$@" 2>/dev/null) > "$root/$out"
  e=$(date +%s.%N)
  echo "TIME $tag $(echo "$e - $s" | bc)" >&2
}
cp out/xs/$c-v0-1.msgpack ../m0/wims69.msgpack
trap 'rm -f ../m0/wims69.msgpack ../m0/pincell-fg.msgpack' EXIT
for round in 1 2; do
  for v in v0 v1 v2 v3; do
    t "out/lc/m0-lattice-$v-$round.out" "m0-lattice $v $round" ../m0 "$root/bin/lattice-$c-$v" pincell.toml
    cp ../m0/pincell-fg.msgpack "out/lc/m0-fg-$v-$round.msgpack"
    t "out/lc/m0-core-$v-$round.out" "m0-core $v $round" ../m0 "$root/bin/core-$c-$v" core.toml
    t "out/lc/1B-$v-$round.out" "1B $v $round" . "$root/bin/lattice-$c-$v" pincell-1B.toml
    mv pincell-1B-fg.msgpack "out/lc/1B-fg-$v-$round.msgpack"
    t "out/lc/2A-$v-$round.out" "2A $v $round" . "$root/bin/lattice-$c-$v" assembly-2A.toml
    mv assembly-2A-fg.msgpack "out/lc/2A-fg-$v-$round.msgpack"
    mv pins-2A.msgpack "out/lc/2A-pins-$v-$round.msgpack"
  done
done
for p in m0-fg 1B-fg 2A-fg 2A-pins; do
  for f in out/lc/$p-*.msgpack; do cmp -s out/lc/$p-v0-1.msgpack $f && echo "SAME $f" || echo "DIFFER $f"; done
done
for p in m0-lattice m0-core 1B 2A; do
  for f in out/lc/$p-v*.out; do cmp -s out/lc/$p-v0-1.out $f && echo "SAME $f" || echo "DIFFER $f"; done
done
