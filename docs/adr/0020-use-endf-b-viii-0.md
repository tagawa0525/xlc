---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0020: 核データは ENDF/B-VIII.0 とし、リポジトリに含めない

## 背景

核データ処理の対象とする評価済み核データを決める必要があった。

## 検討した案

計画には、ほかの評価済みライブラリを比べた記録はない。

## 決定と理由

対象データは ENDF/B-VIII.0 とする（中性子、熱中性子散乱、崩壊、核分裂収率の各サブライブラリ）。ENDF ファイルは
リポジトリに含めず、パスを設定して参照する。

## 確認

約束事で、確かめるテストはない。

## 結果

- 参照解の OpenMC も同じ ENDF/B-VIII.0 から作る（[ADR-0031](0031-verify-against-openmc-with-same-data.md)）
- ENDF/B-VIII.0 にないもの（O in H2O、天然の C など）は、それぞれの ADR で扱う（[ADR-0052](0052-scatter-h-in-h2o-by-sab.md)、[ADR-0079](0079-split-natural-carbon.md)）
