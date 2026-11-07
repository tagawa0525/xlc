---
status: open
labels: [dvc, slurm]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# Slurm でジョブを取り消しても stage_runner.py が残り、commit を続ける

## 問題

- Slurm 724（`python3 examples/stage_runner.py examples/depletion/p2/dvc.yaml xs-nuclide` から始まるジョブ）を 2026-11-07 14:33:02 に `scancel` で
  取り消したあと、runner（PID 4030250）が残った。親は PID 1、cgroup は `0::/system.slice/slurmstepd.scope/system` で、SIGTERM を保留した
  （/proc の SigBlk に SIGTERM）まま、取り消しの時刻に始めた 16 ステージの `dvc commit -f -q`（PID 4142128）を 14:43:59 まで続けて止まった
  （14:35 に手で SIGTERM を送った）。commit を途中で切らない保留は runner の設計（ADR-0193 の実装）だが、Slurm の SIGKILL が届かず、
  取り消した後もジョブの外で 1 コアを使い、worktree の dvc.lock を書き換え続けた
- runner がジョブの cgroup の外にいた理由（slurmstepd の cgroup の扱いか、runner のセッションの作り方か）は確かめていない

関係する ADR: なし

## 完了の条件

- 取り消しで runner と子が止まる（保留中の commit を終えるなら、その間もジョブの中に残る）ことを確かめ、残るなら原因を直す

## 経緯

- 2026-11-07: ADR-0238 のライブラリの作り直し（Slurm 724）を、見積もりを外れたので取り消して分かった
