---
status: open
labels: [dvc, depletion]
opened: 2026-11-07
closed:
milestone:
blocked-by: []
resolution:
---

# 2A の OpenMC との比較が時間を記録しない履歴を deps に持ち、作り直すと記録の時間の内訳を上書きしうる

## 問題

- examples/depletion/p2/dvc.lock の compare@2A は、時間を記録しない Slurm 726（2026-11-07）の比較で timed の記録（Slurm 664）を上書きしない
  ように前の版に戻したので、deps の 2A の履歴の md5 が前の履歴（a2cc87df）のままで、今の履歴（c028d773。k と数密度は同じで時間だけが違う）に
  対して `dvc status` で常に changed になる
- `dvc repro` で比べ直すと、examples/depletion/verification.md の「集合体の燃焼」の時間の内訳の表が描く `p2/compare-2A.json` を、時間を記録しない
  履歴の値で上書きしうる

関係する ADR: なし

## 完了の条件

- 次のどちらかにする: 履歴の時間（`seconds`）を比較の JSON から別のステージ（timed のジョブだけで作るもの）に分ける、または 2A の履歴を
  timed のジョブで作り直して比べ直す

## 経緯

- 2026-11-07: ADR-0238 の 2A の作り直し（Slurm 726）で、比較を前の版に戻したときに分けた
