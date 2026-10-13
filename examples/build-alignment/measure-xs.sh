#!/usr/bin/env bash
# ADR-0157 の xs の測定（ジョブ 302）: 3 つのコミットと 4 つの整列の設定の xs を examples/m0/xs.toml で 2 回ずつ実行し、
# 出力と時間を比べる。先に build.sh <作業ディレクトリ> b873a967 ebe7c100 d39c5db3。Slurm の排他のジョブで実行する
set -euo pipefail
cd "$(dirname "$0")"
root=$(pwd)
# examples/m0 の DVC の出力（wims69.msgpack、pincell-fg.msgpack）を上書きしないよう、あれば止める
for f in ../m0/wims69.msgpack ../m0/pincell-fg.msgpack; do
  if [ -e "$f" ]; then echo "$f exists; run this in a worktree without the M0 outputs" >&2; exit 1; fi
done
rm -rf out/xs; mkdir -p out/xs
trap 'rm -f ../m0/wims69.msgpack' EXIT
for round in 1 2; do
  for v in v0 v1 v2 v3; do
    for c in b873a967 ebe7c100 d39c5db3; do
      tag=$c-$v-$round
      s=$(date +%s.%N)
      (cd ../m0 && "$root/bin/xs-$c-$v" xs.toml 2>/dev/null)
      e=$(date +%s.%N)
      mv ../m0/wims69.msgpack "out/xs/$tag.msgpack"
      echo "TIME $c $v $round $(echo "$e - $s" | bc)"
    done
  done
done
ref=out/xs/b873a967-v0-1.msgpack
n=0
for f in out/xs/*.msgpack; do
  if cmp -s $ref $f; then n=$((n+1)); else echo "DIFFER $f"; fi
done
echo "IDENTICAL $n of $(ls out/xs/*.msgpack | wc -l) to $ref ($(md5sum < $ref))"
