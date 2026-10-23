---
status: open
labels: [xs, performance]
opened: 2026-10-23
closed:
milestone:
blocked-by: []
resolution:
---

# xs の既定のスレッドの数が、許された CPU の中の物理コアでなく、全体の物理コアを許された CPU の数で抑えたものになっている

## 問題

`xs` の既定のスレッドの数は `num_cpus::get_physical().min(available_parallelism)`（crates/xs/src/parallel.rs。ADR-0189）で、
ノード全体の物理コアの数を、許された CPU（ハードウェアスレッド）の数で抑える。r995（物理 16、SMT で 32）で Slurm の `-c 16` のジョブは
許された CPU が 16 のハードウェアスレッド（8 コアにあたる）なので、既定は 16 スレッドになり、8 コアに 2 倍で載る（ADR-0067 の
「物理コアの数」の意図と違う）。排他のジョブと taskset -c 0 では違いは出ない。

examples/stage_runner.py（ADR-0193）は、許された CPU の中の物理コアを数えて予算にし、各ステージに `--threads` を渡すので、この影響を
受けない。`xs` を直接、`-c` のジョブで既定のまま実行したときだけ効く。確かめていない（`-c 16` での測定はない）。

関係する ADR: なし（ADR-0067、ADR-0189 の既定の数え方）

## 完了の条件

許された CPU の中の物理コアを数えるように直すか、今の数え方で足りる理由を ADR に書く。

## 経緯

- 2026-10-23: perf/0157-nuclide-stage-runner の /code-review（2 周目）が runner の同じ数え方を指摘し、runner は直した。xs の側をここに分けた
