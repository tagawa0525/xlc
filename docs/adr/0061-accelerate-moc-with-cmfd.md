---
status: superseded
date: 2026-10-08
supersedes: []
superseded-by: ADR-0083
issues: []
---

# ADR-0061: MOC の固有値計算を、ピンセル単位の粗メッシュによる CMFD で加速する

## 背景

ピンセルではべき乗法に 400 回程度かかる（docs/plans/m2-assembly.md、d59903e、2026-10-06）。17×17 集合体と 69 群で
そのまま回すと、計算時間が実用にならない。

- 速さ

## 検討した案

計画には、CMFD 以外の加速法を比べた記録はない。

## 決定と理由

MOC の固有値計算を、ピンセル単位の粗メッシュによる CMFD で加速する。

## 確認

`cmfd_accelerates_a_pin_cell`、`cmfd_accelerates_a_lattice_with_a_margin`、`cmfd_accelerates_a_lattice_with_vacuum_sides`、
`sweep_currents_balance_each_coarse_cell`。

## 結果

- CMFD の拡散係数の定め方は [ADR-0062](0062-use-optimally-diffusive-cmfd.md)、CMFD の後の更新は [ADR-0063](0063-update-angular-flux-after-cmfd.md)
- M3 で格子を長方形の格子に一般化したとき、粗メッシュのセルをピンセルから格子のセルに変えた（[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）。CMFD で加速する
  ことを含めて、[ADR-0083](0083-generalize-to-rectilinear-lattice.md) にまとめ直した
