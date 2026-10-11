---
status: accepted
date: 2026-10-07
supersedes: [ADR-0051]
superseded-by:
issues: []
requires: []
---

# ADR-0072: ライブラリの σ0 ごとの群定数の重みは、分離共鳴の範囲で減速計算から作る

## 背景

σ0 ごとの群定数は Bondarenko の重み 1/(σt + σ0)（NR 近似）で作っていた（[ADR-0051](0051-weight-shielded-data-by-nr.md)）。この重みは、吸収体自身の共鳴散乱を
中性子束を減らすものとして扱うため、U-238 の低い共鳴（散乱でのエネルギーの損失が共鳴の幅と同じ程度）で遮蔽を強く
見積もる（均質な媒質で 21〜25 群の吸収が −9〜−30%。examples/vera/verification.md、docs/plans/m2-assembly.md、42f67de、
2026-10-07）。σ0 ごとの群定数は、UFG の範囲の外の等価原理（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）と、核データの検証に使う。

- 評価済みデータだけから決まること（手で調整する数表を持たない）
- 群構造を変えてもライブラリを作り直すだけで済むこと

## 検討した案

- Bondarenko の重み（[ADR-0051](0051-weight-shielded-data-by-nr.md)）: 上のとおり遮蔽を強く見積もる
- NJOY GROUPR の flux calculator と同じ減速計算（採用）

## 決定と理由

ライブラリの σ0 ごとの群定数の重みは、分離共鳴の範囲で UFG の減速計算から作る。NJOY GROUPR の flux calculator と
同じく、吸収体と背景の減速材の無限の均質な混合物で、背景の散乱は重み関数 w(E) どおりの滑らかなスペクトルを与えると
し、吸収体の弾性散乱による減速を積分形の減速方程式 φ(E)(σt(E) + σ0) = σ0 w(E) + ∫ σs(E') φ(E') / ((1 − α)E') dE' で
解く。手で調整する補正の数表を持たず、評価済みデータだけから決まるので、群構造を変えてもライブラリを作り直すだけで
済む。

- 計算する範囲は、重み関数の熱中性子の区切り（0.1 eV）以上で最も低い群の境界から、分離共鳴の範囲の上端以下で最も高い
  群の境界まで。群の中で重みの式が切り替わらないので、両者の規格化は群平均に効かない。熱中性子の領域では標的の熱運動が
  効き、静止した標的の減速計算は成り立たないので、その下の中性子束は σ0 によらない重み関数とする（NJOY と同じ。燃料の
  中の熱中性子束のへこみは輸送計算で扱う）。上は Bondarenko の重みのまま
- 分離共鳴の範囲を持たない核種は Bondarenko の重みのまま

## 確認

NJOY GROUPR の flux calculator（重み関数の番号を負にする）と 1e-3 以内で一致（docs/plans/m2-assembly.md、42f67de、
2026-10-07）。OpenMC の均質な無限媒質と比べた（examples/vera/verification.md）。テストは
`narrow_resonance_flux_matches_analytic_average`、`slowing_down_covers_the_resolved_range_in_whole_groups`、
`slowing_down_weight_changes_only_groups_in_its_range`、`below_the_slowing_down_range_the_flux_is_the_weight`、
`finer_grids_do_not_change_the_flux`。

## 帰結

- 散乱核の近似（標的は静止）は [ADR-0070](0070-treat-targets-at-rest.md)
- 格子の細かさは [ADR-0004](0004-keep-solver-settings-per-stage.md)
