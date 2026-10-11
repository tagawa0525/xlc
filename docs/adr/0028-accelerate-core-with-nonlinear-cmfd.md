---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0028: 炉心計算の固有値は非線形 CMFD で加速する

## 背景

炉心計算の固有値の反復を加速する方法を決める必要があった。M2 の NEM は CMFD を持たない。

## 検討した案

計画には、ほかの加速法を比べた記録はない。

## 決定と理由

固有値は非線形 CMFD で加速する。粗メッシュのソルバは lattice と共通のクレートのものを使う（[ADR-0003](0003-extract-shared-computation-on-second-use.md)）。

## 確認

まだ実装していない（M3）。実装したら、加速の有無で固有値と中性子束が一致し、反復が減ることをテストで確かめる。

## 帰結

演算子と固有値の反復は分けておく（[ADR-0037](0037-separate-operators-from-iterations.md)）。
