---
status: open
labels: [dvc, performance]
opened: 2026-11-15
closed:
milestone:
blocked-by: []
resolution:
---

# stage_runner と dvc repro の固定費が p2 の計算より長い

## 問題

- Slurm 887（`python3 examples/stage_runner.py examples/depletion/p2/dvc.yaml xs-nuclide`、545 核種、排他、xlc d6cc90d96、
  2026-11-15）の 1 時間 49 分のうち、16 コアで計算していたのは 31〜58 分の約 27 分だけだった。ログ
  （`~/github/slurm-logs/xlc-p2-library-defaults-887.out`）の runner の時刻で見た内訳:
  - 0〜547 秒: 上流の changed の確認（`dvc repro examples/depletion/p2/dvc.yaml:xs-range` の前）
  - 547〜1864 秒: `xs-range` の実行と、核種の投入の前の処理（deps の指紋など。何に時間がかかったかは確かめていない）
  - 1864〜約 3500 秒: 核種の計算（小さい核種は 1 本 0.05 秒ほど）
  - 約 3500〜6500 秒より後: `dvc commit` が 1 コアで 1 回約 100 秒。まとめた核種の数によらない（64 回で約 1 時間 40 分）
- `dvc commit` は 1 回ごとに 558 ステージの `dvc.yaml` 全体を読み直しているとみられる（確かめていない）。後に続く p2 の
  `dvc repro examples/depletion/p2/dvc.yaml:compare@<履歴>`（Slurm 888〜895）も、計算の前に同じ確認で 10〜20 分かかる見込み
  （見込みで、まだ測っていない）
- CLAUDE.md の方針の「問題の大きさに釣り合わない時間」に当たる。ライブラリを作り直すたびに、計算より長い固定費がかかる

## やること

- 小さい入力で、`dvc commit` 1 回の時間を、ステージの数と deps の数を変えて測り、何に比例するかを決める（cProfile など）
- runner の commit を、計算が終わった核種をまとめて少ない回数で行う形にできるかを見る（ADR-0193 の「commit の前後で deps の
  指紋を比べる」の制約と合わせる）
- 547〜1864 秒の処理の中身を測る

関係する ADR: ADR-0193、ADR-0184
