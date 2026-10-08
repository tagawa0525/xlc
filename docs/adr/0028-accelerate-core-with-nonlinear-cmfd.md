---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0028: 炉心計算の固有値は非線形 CMFD で加速する

## 背景

炉心計算の固有値の反復を加速する方法を決める必要があった。M2 の NEM は CMFD を持たない。

## 検討した案

計画には、ほかの加速法を比べた記録はない。

## 決定と理由

固有値は非線形 CMFD で加速する。粗メッシュのソルバは lattice と共通のクレートのものを使う（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）。

## 確認

M3 で実装したときに、反復の回数と結果の一致で確かめる。

## 結果

演算子と固有値の反復は分けておく（[ADR-0037](0037-separate-operators-from-iterations.md)）。
