---
status: superseded
date: 2026-10-22
supersedes: []
superseded-by: [ADR-0186]
issues: []
requires: []
---

# ADR-0042: ν は MF1 MT452、χ は MF5 MT18 から取り、χ は核分裂源で重み付けして1本にする

## 背景

多群ライブラリの ν と核分裂スペクトル χ の作り方を決める必要があった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

ν は MF1 MT452、核分裂スペクトル χ は MF5 MT18 を使い、χ は核分裂源（νσf と重み関数の積）で重み付けして1本の
ベクトルにする。

## 確認

`chi_is_normalized`、`fissile_nuclide_needs_normalized_chi`。NJOY GROUPR との比較（examples/m1/verification.md）。

## 帰結

遅発中性子の ν とスペクトルは M4 で加える（[ADR-0033](0033-store-kinetics-and-cooling-data-in-m4-tables.md)）。
