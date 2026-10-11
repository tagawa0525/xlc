---
status: accepted
date: 2026-10-08
supersedes: [ADR-0059, ADR-0061]
superseded-by:
issues: []
requires: []
---

# ADR-0083: 格子を列ごとの幅と行ごとの高さを持つ長方形の格子に一般化し、セルの形を増やす

## 背景

M2 の形状は、正方格子に並べたセルと集合体の外周の水ギャップである（[ADR-0059](0059-use-square-lattice-geometry.md)）。M3 では、シンブルプラグ、支持格子、
バッフル、反射体の水を表し、集合体の外の反射体を解く。M2 の CMFD は、ピンセル単位の粗メッシュで MOC を加速していた
（[ADR-0061](0061-accelerate-moc-with-cmfd.md)）。バッフル（2.85 cm）のようにピンのピッチに合わない厚さがある。

## 検討した案

- 正方格子のまま、水ギャップを特別な領域として持つ（[ADR-0059](0059-use-square-lattice-geometry.md)）: ピンのピッチに合わない厚さを表せない
- 長方形の格子に一般化し、余白も薄い行と列にする（採用）

## 決定と理由

セルの形を増やす: 中実の棒（シンブルプラグ）、セルの外周の枠（支持格子。VERA と同じく質量を等体積の枠にする）、
ピンでない長方形の領域（バッフル、反射体の水）。集合体の外の反射体を解くため、格子の幅が場所で違ってよいようにする。

- 格子は、列ごとの幅と行ごとの高さを持つ長方形の格子に一般化する。ピンセルは幅と高さがそのピッチに等しいマスに置き、
  ほかのマスには均質な長方形のセル（小さな長方形に分けられる）を置く。バッフルのようにピンのピッチに合わない厚さを、
  近似なしで表すため
- 集合体の周りの余白（水のギャップ）は特別な領域をやめ、薄い行と列として表す。余白と反射体を同じ仕組みで扱い、形状の
  場合分けを減らすため
- トラッキングは長方形の範囲の cyclic tracking にする
- MOC の固有値計算は CMFD で加速する（[ADR-0061](0061-accelerate-moc-with-cmfd.md) から引き継ぐ。ピンセルでべき乗法に 400 回程度かかり、17×17 集合体と
  69 群ではそのままでは計算時間が実用にならないため）。粗メッシュのセルは、ピンセルに代えて格子のセルとする（余白の
  薄いセルで加速が乱れるなら、隣のセルと合わせる）
- ピンセルの領域は、同心円のリングに加えて方位角方向に分割できる（[ADR-0059](0059-use-square-lattice-geometry.md) から引き継ぐ）。平坦線源の近似の誤差を
  領域の細分で収束させて確かめるため

## 確認

`columns_and_rows_may_have_their_own_sizes`、`rectangle_track_volumes_match_analytic`、`rectangle_tracks_close_by_reflection`、
`framed_cell_track_volumes_match_analytic`、`margin_regions_and_volumes`、`segments_cross_the_margin`。CMFD は
`cmfd_accelerates_a_pin_cell`、`cmfd_accelerates_a_lattice_with_a_margin`、`cmfd_accelerates_a_lattice_with_vacuum_sides`、
`sweep_currents_balance_each_coarse_cell`。

## 帰結

支持格子の枠の遮蔽は [ADR-0086](0086-shield-structural-lumps.md)、[ADR-0088](0088-add-frame-ring-to-equivalent-cell.md)。
