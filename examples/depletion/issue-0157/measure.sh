#!/usr/bin/env bash
# ADR-0189 の測定（Slurm 479 で実行した手順）: 全コアで 1 スレッドずつ同じ計算を実行したときの遅れ、核種ごとのステージの形
# （1 核種ずつ順に、既定のスレッドの数）、全核種を 1 プロセスで（案 B の形）、CPU 時間の順に同時に投入する形（案 C の試行）。
# 入力は ../xs.toml の 65 核種（このディレクトリの xs.toml はそのコピー）と ../xs-range.toml。DVC を通さない。Slurm 479 は
# 同じ手順を trials/ に置いて実行した（パスだけが違う）。
# Usage: examples/depletion/issue-0157/measure.sh（cargo build --release -p xlc-xs の後、排他の Slurm のジョブで）
set -euo pipefail
cd "$(dirname "$0")"
xs=$(readlink -f ../../../target/release/xs)
T=/run/current-system/sw/bin/time
rm -rf out && mkdir -p out/one out/stages out/runner
echo "== slowdown: Cm-244 565 K, --threads 1"
$T -f "single elapsed=%e user=%U" taskset -c 0 "$xs" --threads 1 nuclide one.toml ../xs-range.toml Cm-244 out/one/single.msgpack 2>&1 | grep -v " INFO \| WARN " || true
for k in $(seq 0 15); do
  $T -f "copy$k elapsed=%e user=%U" taskset -c "$k" "$xs" --threads 1 nuclide one.toml ../xs-range.toml Cm-244 "out/one/copy$k.msgpack" 2>&1 | grep -v " INFO \| WARN " &
done
wait
md5sum out/one/*.msgpack
echo "== stages: one nuclide per process, default threads"
start=$(date +%s.%N)
sed -nE 's/^\[nuclide\.([^]]+)\]$/\1/p' xs.toml | while read -r n; do
  $T -f "nuclide $n elapsed=%e user=%U sys=%S maxrss=%M" "$xs" nuclide xs.toml ../xs-range.toml "$n" "out/stages/$n.msgpack" 2>&1 | grep -v " INFO \| WARN "
done | tee out/stages-times.txt
echo "stages total elapsed=$(echo "$(date +%s.%N) - $start" | bc)"
md5sum out/stages/*.msgpack
echo "== whole run"
$T -f "whole elapsed=%e user=%U sys=%S maxrss=%M" "$xs" xs.toml 2>&1 | grep -v " INFO \| WARN "
md5sum all.msgpack
echo "== runner"
sed -nE 's/^nuclide (\S+) elapsed=\S+ user=(\S+) sys=(\S+).*/\1 \2 \3/p' out/stages-times.txt | awk '{print $1, $2 + $3}' > out/cpu.txt
python3 runner_trial.py "$xs" . out/cpu.txt "$(pwd)/out/runner"
md5sum out/runner/*.msgpack
