---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
---

# ADR-0060: 集合体全体を1組の cyclic tracking で解く

## 背景

ピンセルでは cyclic tracking を使っている（[ADR-0043](0043-track-pin-cells-cyclically.md)）。集合体と C5G7 に広げる方法を決める必要があった。

## 検討した案

- モジュール単位のトラッキング: メモリが少ない。メモリが問題になるまで使わない
- 集合体全体で1組のトラック（採用）

## 決定と理由

トラッキングは、ピンセルで使っている cyclic tracking を集合体全体に広げる（全体で1組のトラック）。モジュール単位の
トラッキングは、メモリが問題になるまで使わない。

## 確認

`lattice_tracks_reflect_at_the_outer_boundary`、`links_join_tracks_at_reflected_points`、`vacuum_sides_have_no_links`、
`vacuum_lattice_balances_production_with_absorption_and_leakage`。C5G7 との比較（examples/c5g7/verification.md）。

## 結果

- M3 で、長方形の範囲の cyclic tracking にする（[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）。全体で1組のトラックであることは変わらない
- 外周の境界条件（鏡面反射、真空、アルベド）は [ADR-0093](0093-choose-lattice-boundary-conditions.md)
