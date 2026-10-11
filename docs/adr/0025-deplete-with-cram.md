---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0025: 燃焼計算は CRAM で解き、予測子・修正子法で時間を進める

## 背景

格子計算の燃焼計算の解法を決める必要があった（M4）。

## 検討した案

計画には、ほかの解法を比べた記録はない。

## 決定と理由

燃焼計算は CRAM（Chebyshev 有理近似）で解き、予測子・修正子法で時間を進める。燃焼チェーンは xs が ENDF の崩壊・
収率データから作ったもの（全核種）を、輸送計算に必要な範囲に縮約して使う（[ADR-0034](0034-build-full-depletion-chain.md)）。燃焼のソルバは共通のクレートに
置く（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）。

## 確認

まだ実装していない（M4）。M4 で格子燃焼のベンチマークと比べる。

## 結果

CRAM の疎行列の LU 分解は faer で解く（[ADR-0007](0007-use-faer-for-linear-algebra.md)）。
