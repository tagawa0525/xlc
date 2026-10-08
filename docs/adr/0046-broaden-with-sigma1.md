---
status: accepted
date: 2026-10-06
supersedes: [ADR-0038]
superseded-by:
issues: []
---

# ADR-0046: ドップラー拡がりは SIGMA1 法で求める

## 背景

M0 の断面積は 0 K だった（[ADR-0038](0038-use-zero-kelvin-cross-sections-in-m0.md)）。M1 では、温度に依存する多群ライブラリを作る。

- NJOY と直接比べられること
- 速さ

## 検討した案

- 共鳴を温度ごとに再構成し直す: 遅い
- 0 K の線形化した断面積を、区間ごとの解析的な積分（NJOY BROADR と同じ SIGMA1 法）で広げる（採用）

## 決定と理由

0 K の線形化した断面積を、区間ごとの解析的な積分（NJOY BROADR と同じ SIGMA1 法）で広げる。共鳴を温度ごとに再構成し
直すより速く、NJOY と直接比べられる。

M1 の計画は、この項目に温度の組（293.6, 600, 900, 1200 K）も書いていた。温度の選び方は別の決定として [ADR-0047](0047-choose-temperatures-per-nuclide.md) に
分けた（のちに [ADR-0080](0080-interpolate-sab-in-temperature.md) で置き換え）。

## 確認

`constant_cross_section_matches_analytic_broadening`、`one_over_v_is_unchanged`。NJOY BROADR との比較
（`crates/xs/examples/compare_njoy_broadr.rs`、examples/m1/README.md）。

## 結果

温度の選び方は [ADR-0047](0047-choose-temperatures-per-nuclide.md)（M1）と、それを置き換えた [ADR-0080](0080-interpolate-sab-in-temperature.md)。
