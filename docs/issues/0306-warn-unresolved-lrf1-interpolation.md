---
type: Issue
status: open
labels: [xs]
opened: 2026-11-20
closed:
milestone:
blocked-by: []
resolution:
---

# 非分離共鳴の LRF=1、LFW=1 の断面積を、xs は lin-lin、NJOY RECONR は log-log で補間し、差を warn も記録もしていない

## 問題

- xs は、非分離共鳴の LRF=1 の平均断面積を格子点（核分裂幅のエネルギー ES と NJOY と同じ補助格子）で求め、その間を lin-lin で
  補間する（crates/xs/src/endf/mf2.rs の `read_unresolved_lrf1` が系列の `int` を LinLin にし、crates/xs/src/pointwise.rs の
  `replace_with_unresolved` がその補間則を使う）。ENDF-102（BNL-224854-2023-INRE）の 2.4.1 の Case B（LRF=1、LFW=1）は INT を持たず
  lin-lin とし、2.4.2 は断面積を補間するとする
- NJOY2016 の RECONR は、LRF=1 の平均断面積の表を `intunr=5`（log-log）で補間する（`rdfil2`、`sigunr`）。LRF=2 はファイルの INT
- 差は格子点の間に出る。ENDF/B-VII.0 の Pu-238（200 eV〜10 keV）の 0 K の再構成で、格子点の上では NJOY と 1e-5 以下で合い、間の点で
  捕獲が最大 1.0%、核分裂が 0.36%、全断面積が 0.40% 大きい（examples/depletion/verification.md の「版ごとのライブラリの分離共鳴の
  SLBW」。212.5 eV の全断面積は、200 eV と 250 eV の NJOY の値の lin-lin の補間 16.16549 b と log-log の補間 16.0541 b が、xs と
  NJOY の値に 1e-5 で合う）
- NJOY と違う扱いなのに `tracing::warn!` を出さず、ADR に理由と差の大きさが無い（ADR-0058 は LRF=1 を LRF=2 と同じ形に読み替える
  ことだけを書く）。リポジトリの CLAUDE.md の「方針」の、参照コードと違う扱いの項に合わない

関係する ADR: ADR-0058

## 完了の条件

- [ ] lin-lin を保つか決め、理由と NJOY との差の大きさを ADR に書く（ADR-0058 を置き換えるか、新しい ADR）
- [ ] NJOY と違う補間をする範囲で warn を出す

## 経緯

- 2026-11-20: SLBW の再構成を NJOY と比べたとき（#0305）、Pu-238 の非分離共鳴の範囲の格子点の間で差が見つかり、起こした
