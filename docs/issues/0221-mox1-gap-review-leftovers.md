---
status: open
labels: [verification]
opened: 2026-11-03
closed:
milestone: M4
blocked-by: []
resolution:
---

# MOX-1 の 0 の点の切り分けの記録の周で残った指摘

## 問題

ブランチ feat/0152-mox1-gap（#0152 の 2026-11-03 の MOX-1 の 0 の点の切り分けの項、#0219）のレビューは 3 周を回し、3 周目で直した
head（9cd0d1772）の記録の周に出た指摘のうち、文書の事実の誤り（ADR-0205 の差と比べる符号、Slurm 665 の一時停止の語、引いた値のコミットと日付）は
同じブランチで直した。残りは直していない:

- #0001 に、MOX-1 の 0 の点の差の主な寄与（Pu-242 の 0.625〜4 eV −178.9 pcm、29 群の +11.29%）を書いていない（このブランチの変える範囲の外）
- `examples/depletion/issue-0152/fuel_xs_compare.py` は `issue-0133/groupxs.py` と同じ実効断面積の比べ方を持つ（標準偏差のキーだけが違う）
- `examples/depletion/issue-0152/run.sh`: リポジトリの直下からしか動かない、失敗した再実行が追跡しているログを空にする、`lattice` に入力の名前が
  無くても 0 で終わる、`bin/lattice`（DVC の build）と `target/release/examples`（cargo）の版がそろうことを確かめない

関係する ADR: なし

## 完了の条件

issue #0001 の経緯に #0152 の結果を引き、run.sh の 4 件を直すか、run.sh を消して手順を記録に書き、fuel_xs_compare.py を groupxs.py に
まとめるか残すかを決めて閉じる。

## 経緯

- 2026-11-03: ブランチ feat/0152-mox1-gap の記録の周（9cd0d1772）で出た。文書の事実の誤りの 3 件は同じブランチで直した
