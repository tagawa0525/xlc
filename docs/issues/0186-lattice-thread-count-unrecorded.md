---
status: open
labels: [lattice, xs, docs]
opened: 2026-10-28
closed:
milestone:
blocked-by: []
resolution:
---

# xs と lattice のスレッドの数がログに出ず、記録の「32 スレッド」を確かめられない

## 問題

xs と lattice は、`RAYON_NUM_THREADS` が無ければ、物理コアの数をプロセスに許された CPU の数で抑えた数のスレッドで解く
（ADR-0067。crates/lattice/src/main.rs の `use_physical_cores`、crates/xs/src/parallel.rs）。r995 は 1 ソケット 16 コアで、
1 コアに 2 スレッド（`lscpu`）なので、Slurm の排他のジョブ（32 CPU）では 16 になる。ログ（`lattice-*.log` など）には
スレッドの数が出ないので、記録を書く側は Slurm の割り当ての 32 をスレッドの数と取り違えうる。次の記述は lattice の実行を
「32 スレッド」と書いているが、ログからは確かめられない:

- ADR-0195 の背景: 2A の燃焼（examples/depletion/p2/lattice-2A.log、0f97cea6、2026-10-23）の 1 輸送の MOC が「7.4 秒（32 スレッド）」。
  examples/depletion/p2/dvc.yaml のステージは `RAYON_NUM_THREADS` を付けない（`git log -S RAYON_NUM_THREADS` で履歴にも無い）。
  Slurm に投げたコマンドに付けたかは確かめていない
- examples/vera/verification.md の「格子計算の並列の効率（2A）」: 「燃焼の本番（2A、…、排他のジョブで 32 スレッド）」
- #0108 の比較の行（「32 スレッドで 7 秒」。同じ行が `RAYON_NUM_THREADS` を設定したかは確かめていないと断っている）

取り違えが確かな例は 1 つある。2026-10-28 に、親のセッションが #0181 の実装の依頼に「2A は 32 CPU で 39 本の解くピンを並列に
解き（ADR-0198）」と書き、ユーザーへの答えにも「今は rayon の既定で 32 スレッド使っています」と書いた。ADR-0198 は 16 スレッドと
書いている。担当が main.rs で確かめ、#0181 の一連の 2A と 1C の測定（Slurm 573〜578、580）は 16 スレッドで、2A の燃焼を SMT の
32 スレッドで測ったのは Slurm 579 が初めてと分かった（ADR-0213 の確認、#0181 の経緯。16 と 32 の比較は ADR-0067 が 2B の P3 で
先に測っている）。

関係する ADR: なし

## 完了の条件

- xs と lattice が始めにスレッドの数（rayon のプールの数と、それを決めたもの: 既定か `RAYON_NUM_THREADS`）をログに出す。
  bin が変わるので、ステージは changed になる（結果が要るときに再実行する。ADR-0130）
- 上の記述のうち、ログ、Slurm のジョブのスクリプト（`sacct` の SubmitLine）で数が確かめられるものは直し、確かめられないものは
  「確かめていない」を添える（accepted の ADR-0195 は本文を書き換えず、ADR の運用の追記の形で）
- CLAUDE.md の「計算」の、スレッドの数の行（「ログに出ない」）を、ログを読む形に直す

## 経緯

- 2026-10-28: /retro 11 回目（#0181 の実装の区間）で、依頼の前提の誤り（32 CPU を 32 スレッドとした）の原因を探して見つけた
- 2026-10-28: ADR-0214 で lattice の既定を論理 CPU の数にし、始めにスレッドの数をログに出すようにした（`thread pool threads=`）。
  xs の数のログ、上の記述の直し、CLAUDE.md の xs の行は残る
