---
status: superseded
date: 2026-10-08
supersedes: []
superseded-by: ADR-0083
issues: []
requires: []
---

# ADR-0059: 集合体の形状は、正方格子に並べたセルと集合体の外周の水ギャップとする

## 背景

M2（2026-10-06）で、17×17 集合体と C5G7 を解く形状を決める必要があった。

## 検討した案

計画には、ほかの形状を比べた記録はない。

## 決定と理由

形状は正方格子に並べたセル（同心円で区切ったピンセル、または均質なセル）と、集合体の外周の水ギャップとする。
ピンセルの領域は、同心円のリングに加えて方位角方向に分割できるようにする。平坦線源の近似の誤差を領域の細分で
収束させて確かめるため。

## 確認

`lattice_track_volumes_match_analytic`、`sector_cell_track_volumes_match_analytic`、`gap_surrounds_the_lattice`。C5G7 との比較
（examples/c5g7/verification.md）。

## 帰結

M3 で、列ごとの幅と行ごとの高さを持つ長方形の格子に一般化し、水ギャップを薄い行と列にした（[ADR-0083](0083-generalize-to-rectilinear-lattice.md)）。
方位角方向の分割は引き続き使う。
