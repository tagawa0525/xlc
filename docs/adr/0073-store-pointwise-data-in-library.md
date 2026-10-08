---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: ["#0001"]
---

# ADR-0073: ライブラリに、ポイントワイズの断面積、分離共鳴の上端、ポテンシャル散乱を持たせる

## 背景

lattice のピンごとの減速計算（[ADR-0068](0068-solve-pinwise-slowing-down.md)）はポイントワイズの断面積を使う。lattice は ENDF を読まない（[ADR-0002](0002-keep-stage-crates-independent.md)）ので、xs がライブラリに
入れる必要がある。

## 検討した案

- 共鳴核種の、UFG の範囲の全断面積・弾性散乱・捕獲・核分裂だけを入れる。計画で一時採った（2026-10-07、cc96628）。
  実装するときに下の形に変えた（実装 3361df0、計画 0ebdc01、2026-10-07）
- 全核種の、ライブラリ全体で共通の範囲のポイントワイズの断面積を入れる（採用）

## 決定と理由

- ライブラリには、全核種の温度ごとのドップラー拡がりのポイントワイズの断面積（全断面積・弾性散乱・核分裂・νσf。それ以外の反応は
  減速計算では吸収とみなす）を、ライブラリ全体で共通の範囲で加える。範囲は、熱中性子散乱の打ち切り（4 eV）以上で
  最も低い群の境界から、全核種の分離共鳴の上端のうち最も高いもの以下で最も高い群の境界まで
- 下限をライブラリの UFG の重みの下限（0.1 eV）より高くするのは、H in H2O が 4 eV 未満で S(α,β)（束縛原子）を使い、
  静止した自由原子の散乱核が成り立たないため。セルの UFG もこの下限から計算し、その下は範囲の外の規則（[ADR-0074](0074-shield-outside-range-by-equivalence.md)）に
  従う
- ライブラリの各核種は、分離共鳴の上端（セルの UFG の範囲を決める）と、s 波の散乱半径 a からポテンシャル散乱 4π a² を
  持つ（Reich-Moore は l=0 の APL、0 なら AP。MLBW は AP。R-matrix limited は弾性散乱の l=0 のチャネルの APE）。
  M2 の計画は「Reich-Moore と MLBW は l=0 の APL」と書いていたが、ENDF の MLBW（LRF=2）の l ごとのレコードには
  APL がないので、MLBW は AP を使う（`crates/xs/src/process.rs` の `s_wave_radius`）

## 確認

`keeps_the_pointwise_data`、`pointwise_data_covers_the_range_on_the_total_grid`、`resolved_range_has_its_top_and_potential_scattering`、
`rejects_a_zero_s_wave_radius`、`rejects_differing_s_wave_radii`。

## 結果

4 eV より下の共鳴（WIMS-69 の 2.6〜3.3 eV の群に入る Pu-242 の 2.67 eV など）は減速計算で遮蔽しない（#0001）。
