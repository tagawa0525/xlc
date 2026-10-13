#!/usr/bin/env bash
# ADR-0157 の 2A の追加の測定（ジョブ 307）: 6b39ad13 の整列なし（v0）と -align-loops=64（v3）の lattice を、VERA 2A の集合体で
# 交互に 6 回ずつ走らせ、出力と時間を比べる。先に build.sh。Slurm の排他のジョブで走らせる
set -euo pipefail
cd "$(dirname "$0")"
root=$(pwd)
c=6b39ad13
rm -rf out/2A; mkdir -p out/2A
for round in 1 2 3 4 5 6; do
  for v in v0 v3; do
    s=$(date +%s.%N)
    "$root/bin/lattice-$c-$v" assembly-2A.toml > "out/2A/$v-$round.out" 2>/dev/null
    e=$(date +%s.%N)
    echo "TIME 2A $v $round $(echo "$e - $s" | bc)"
    mv assembly-2A-fg.msgpack "out/2A/fg-$v-$round.msgpack"
    mv pins-2A.msgpack "out/2A/pins-$v-$round.msgpack"
  done
done
for p in fg pins; do
  for f in out/2A/$p-*.msgpack; do cmp -s out/2A/$p-v0-1.msgpack $f && echo "SAME $f" || echo "DIFFER $f"; done
done
for f in out/2A/v*.out; do cmp -s out/2A/v0-1.out $f && echo "SAME $f" || echo "DIFFER $f"; done
