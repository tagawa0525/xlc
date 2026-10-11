---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0024: 集合体全体の輸送計算は 2D MOC とする

## 背景

格子計算の輸送計算の手法を決める必要があった。

## 検討した案

計画には、ほかの手法を比べた記録はない。

## 決定と理由

集合体全体の輸送計算は 2D MOC（特性曲線法）とする。円柱形の燃料ピンを近似なしで扱える現在の主流の手法だからである。

## 確認

M0 で、均質なセルの解析解（無限媒質の多群固有値）との一致と、トラック間隔・角度の収束性を確かめた
（`homogeneous_cell_matches_infinite_medium`、`refinement_converges`）。M2 で C5G7 と VERA Problem 1〜2 と比べた
（examples/c5g7/verification.md、examples/vera/verification.md）。

## 帰結

トラッキングは [ADR-0043](0043-track-pin-cells-cyclically.md)、[ADR-0060](0060-track-whole-assembly-cyclically.md)、反復の加速は [ADR-0083](0083-generalize-to-rectilinear-lattice.md)（M2 では [ADR-0061](0061-accelerate-moc-with-cmfd.md)）。
