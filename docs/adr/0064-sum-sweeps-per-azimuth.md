---
status: accepted
date: 2026-10-07
supersedes: []
superseded-by:
issues: []
---

# ADR-0064: トラックのスイープは rayon で並列化し、領域ごとの和を方位角ごとに持って方位角の順に足す

## 背景

トラックのスイープをスレッドで並列化する。浮動小数点の和の順序が変わると、結果がスレッドの数で変わる。

- 結果の再現性（ビット単位）
- 仕事の量の釣り合い

## 検討した案

計画には、ほかの分け方を比べた記録はない。

## 決定と理由

トラックのスイープはスレッドで並列化する（rayon）。領域ごとの和は方位角ごとに持ち、方位角の順に足す。結果がスレッドの数や
仕事の分け方によらず、ビット単位で同じになるようにするため。方位角ごとのトラックの長さの和はほぼ等しいので、仕事の
量も揃う。

## 確認

`result_does_not_depend_on_the_thread_count`、`anisotropic_result_does_not_depend_on_the_thread_count`。

## 結果

方位角ごとの集計のメモリが要る（非等方散乱では [ADR-0075](0075-solve-moc-with-anisotropic-scattering.md)）。
