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

関係する ADR: なし

## やること

- 2M の履歴を読む作業（履歴からのテーブル、分岐、warm-start の比較など）を始める前に、
  `dvc repro examples/depletion/p2/dvc.yaml:compare@2M` を Slurm の timed で流し、`dvc.lock` と `compare-2M.json` をコミットする
- examples/depletion/verification.md の「集合体の燃焼」の「計算」「時間」と「2M」の節を書き直す（今は前の本番の数値）

## 完了の条件

- `lattice@2M` と `compare@2M` が今の `bin/lattice` で relock され、記録の 2M の節が新しい結果になっている

## 経緯

- 2026-10-24: p2 の再実行で、2P の後に続けて投げた 2M（Slurm 527）を、ステージを始める前に取り消した
- 2026-11-05: 2M は ADR-0226 の後の 8 本の再実行で流された（a1ae3b0a7。ジョブを実行したコミットはブランチ chore/rerun-after-adr-0226 の 9069e1cb、
  Slurm 668。マージ 029d2958c）。`compare-2M.json` と examples/depletion/p2/dvc.lock が更新され、記録（examples/depletion/verification.md）の
  「集合体の燃焼」の「計算」と「2M」の節がその実行を引く。ただし、その後に `bin/lattice` が作り直された（7db11a9d6、d77dcf60a。直下の dvc.lock の
  `build` の md5 7e959698 と、p2 の dvc.lock の `lattice@*` の md5 48387d65）ので、`lattice@2M` はほかの `lattice@*` とともに再び deps の changed で、完了の条件の
  「今の `bin/lattice` で relock」は満たさない。668 は割り込まれうるジョブで、2M の時間は記録の「時間」に無い（記録の「時間」の節の本文）。閉じずに残す
  （ブランチ docs/m4-plan-progress）
