#!/usr/bin/env bash
# ADR-0157 の M0 の炉心と VERA 1B の繰り返しの測定: 6b39ad13 の整列なし（v0）と -align-loops=64（v3）を、回ごとに順を
# 入れ替えて 10 回ずつ走らせ、出力と時間を比べる。先に build.sh と measure-xs.sh（M0 のライブラリを使う）。
# Slurm の排他のジョブで走らせる
set -euo pipefail
cd "$(dirname "$0")"
root=$(pwd)
c=6b39ad13
# examples/m0 の DVC の出力（wims69.msgpack、pincell-fg.msgpack）を上書きしないよう、あれば止める
for f in ../m0/wims69.msgpack ../m0/pincell-fg.msgpack; do
  if [ -e "$f" ]; then echo "$f exists; run this in a worktree without the M0 outputs" >&2; exit 1; fi
done
rm -rf out/repeat; mkdir -p out/repeat
trap 'rm -f ../m0/wims69.msgpack ../m0/pincell-fg.msgpack' EXIT
cp out/xs/$c-v0-1.msgpack ../m0/wims69.msgpack
(cd ../m0 && "$root/bin/lattice-$c-v0" pincell.toml > /dev/null 2> "$root/out/repeat/m0-lattice.log")
for round in 1 2 3 4 5 6 7 8 9 10; do
  if (( round % 2 )); then order=(v0 v3); else order=(v3 v0); fi
  for v in "${order[@]}"; do
    s=$(date +%s.%N)
    (cd ../m0 && "$root/bin/core-$c-$v" core.toml) > "out/repeat/core-$v-$round.out" 2> "out/repeat/core-$v-$round.log"
    e=$(date +%s.%N)
    dt=$(echo "$e - $s" | bc)
    echo "TIME m0-core $v $round $dt"
    s=$(date +%s.%N)
    "$root/bin/lattice-$c-$v" pincell-1B.toml > "out/repeat/1B-$v-$round.out" 2> "out/repeat/1B-$v-$round.log"
    e=$(date +%s.%N)
    dt=$(echo "$e - $s" | bc)
    echo "TIME 1B $v $round $dt"
    mv pincell-1B-fg.msgpack "out/repeat/1B-fg-$v-$round.msgpack"
  done
done
for p in core 1B; do
  for f in out/repeat/$p-v*.out; do cmp -s "out/repeat/$p-v0-1.out" "$f" && echo "SAME $f" || echo "DIFFER $f"; done
done
for f in out/repeat/1B-fg-*.msgpack; do cmp -s out/repeat/1B-fg-v0-1.msgpack "$f" && echo "SAME $f" || echo "DIFFER $f"; done
