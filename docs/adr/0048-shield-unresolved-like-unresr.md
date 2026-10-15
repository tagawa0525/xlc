---
status: superseded
date: 2026-10-06
supersedes: [ADR-0038]
superseded-by: ADR-0177
issues: []
requires: []
---

# ADR-0048: 非分離共鳴の自己遮蔽は NJOY UNRESR と同じ式で求める

## 背景

M0 では、非分離共鳴の範囲は MF3 の無限希釈の平均断面積を使っていた（[ADR-0038](0038-use-zero-kelvin-cross-sections-in-m0.md)）。U-238 の非分離共鳴（20〜149 keV）の
自己遮蔽は無視できない大きさがある。

## 検討した案

計画には、ほかの方式を比べた記録はない。

## 決定と理由

共鳴パラメータ（LRF=2）から、バックグラウンド断面積と温度ごとの実効断面積を NJOY UNRESR と同じ式（ETOX 法）で求める。U-238 の
非分離共鳴（20〜149 keV）の自己遮蔽は無視できない大きさがあるため。LSSF=0 の核種（U-234, U-236）の無限希釈の断面積も、
NJOY RECONR と同じ式で共鳴パラメータから求める。

## 確認

`unresolved_shielding_lowers_capture_at_small_background`、`lssf0_adds_shielding_difference`、
`lssf1_scales_by_shielding_factor`。NJOY UNRESR との比較（`crates/xs/examples/compare_njoy_shielding.rs`、
examples/m1/verification.md）。

## 帰結

- 幅のゆらぎの期待値の求め方は [ADR-0049](0049-integrate-width-fluctuations-adaptively.md)
- LRF=1 は [ADR-0058](0058-read-lrf1-as-lrf2-with-l3.md)
