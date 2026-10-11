---
status: accepted
date: 2026-10-06
supersedes: []
superseded-by:
issues: []
requires: []
---

# ADR-0043: ピンセルの MOC は外周を鏡面反射とし、cyclic tracking で解く

## 背景

M0 の格子計算は、ピンセルの 2D MOC である（[ADR-0024](0024-solve-lattice-with-2d-moc.md)）。正方形の外周の境界条件と、トラックの置き方を決める必要が
あった。

## 検討した案

計画には、ほかの案を比べた記録はない。

## 決定と理由

ピンセルの 2D MOC を、正方形の外周で鏡面反射とする。鏡面反射を厳密に扱うため、方位角はトラックが外周でつながる
角度に補正する（cyclic tracking）。

## 確認

均質なセルでの解析解（無限媒質の多群固有値）との一致と、トラック間隔・角度の収束性
（`homogeneous_cell_matches_infinite_medium`、`refinement_converges`、`corrected_angles_are_symmetric`、
`following_links_returns_to_start`）。

## 結果

M2 で集合体全体に広げた（[ADR-0060](0060-track-whole-assembly-cyclically.md)）。
