---
status: open
labels: [depletion, verification]
opened: 2026-10-24
closed:
milestone: M4
blocked-by: []
resolution:
---

# 2M の集合体の燃焼を今の履歴の形式で再実行する

## 問題

履歴の形式に `shieldings` が足された（ADR-0197）ので、`examples/depletion/p2/lattice-2M.msgpack`（前の本番。2ff103e1、Slurm 480）は
今の lattice の `[history]` で読めない。ほかの 7 本（2A、2F、2N、2O、2P、2A-density、2A-rods）は 2026-10-24 に今の `bin/lattice` で
再実行した（ブランチ chore/p2-rerun-after-speedups）が、2M は流していない（約 2 時間。ユーザーの判断で、要るときまで後回し）。
`lattice@2M` と `compare@2M` は `dvc status` で changed（`bin/lattice`）のままで、`dvc repro -P` などで下流を流すと、2M の
再実行が走る。

## やること

- 2M の履歴を読む作業（履歴からのテーブル、分岐、warm-start の比較など）を始める前に、
  `dvc repro examples/depletion/p2/dvc.yaml:compare@2M` を Slurm の timed で流し、`dvc.lock` と `compare-2M.json` をコミットする
- examples/depletion/verification.md の「集合体の燃焼」の「計算」「時間」と「2M」の節を書き直す（今は前の本番の数値）

## 完了の条件

- `lattice@2M` と `compare@2M` が今の `bin/lattice` で relock され、記録の 2M の節が新しい結果になっている

## 経緯

- 2026-10-24: p2 の再実行で、2P の後に続けて投げた 2M（Slurm 527）を、ステージを始める前に取り消した
