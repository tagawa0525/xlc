---
status: accepted
date: 2026-10-08
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0079: VERA の天然の C は、天然の存在比で C-12 と C-13 に分ける

## 背景

ENDF/B-VIII.0 には天然の C がなく、C-12 と C-13 だけがある。VERA の仕様は天然の C の個数密度で与える。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

VERA の天然の C の個数密度は、天然の存在比で C-12 と C-13 に分ける。

## 確認

約束事で、確かめるテストはない。

## 結果

OpenMC の参照解も同じ分け方にする（[ADR-0031](0031-verify-against-openmc-with-same-data.md)）。
