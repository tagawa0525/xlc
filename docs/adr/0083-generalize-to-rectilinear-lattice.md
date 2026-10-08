---
status: accepted
date: 2026-10-08
supersedes: [ADR-0059]
superseded-by:
issues: []
---

# ADR-0083: 格子を列ごとの幅と行ごとの高さを持つ長方形の格子に一般化し、セルの形を増やす

## 背景

M2 の形状は、正方格子に並べたセルと集合体の外周の水ギャップである（[ADR-0059](0059-use-square-lattice-geometry.md)）。M3 では、シンブルプラグ、支持格子、
バッフル、反射体の水を表し、集合体の外の反射体を解く。バッフル（2.85 cm）のようにピンのピッチに合わない厚さがある。

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
- トラッキングは長方形の範囲の cyclic tracking にし、CMFD の粗メッシュのセルは格子のセルとする（余白の薄いセルで加速が
  乱れるなら、隣のセルと合わせる）
- ピンセルの方位角方向の分割（[ADR-0059](0059-use-square-lattice-geometry.md)）は引き続き使う

## 確認

`columns_and_rows_may_have_their_own_sizes`、`rectangle_track_volumes_match_analytic`、`rectangle_tracks_close_by_reflection`、
`framed_cell_track_volumes_match_analytic`、`margin_regions_and_volumes`、`segments_cross_the_margin`。

## 結果

支持格子の枠の遮蔽は [ADR-0086](0086-shield-structural-lumps.md)、[ADR-0088](0088-add-frame-ring-to-equivalent-cell.md)。
