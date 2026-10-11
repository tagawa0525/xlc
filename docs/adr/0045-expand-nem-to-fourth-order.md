---
status: superseded
date: 2026-10-09
supersedes: []
superseded-by: ADR-0102
issues: []
requires: []
---

# ADR-0045: NEM は4次多項式で展開し、横方向の漏れは2次で近似する

## 背景

炉心計算は 3D の2群 NEM とした（[ADR-0027](0027-solve-core-with-nem.md)）。展開の次数と横方向の漏れの近似を決める必要があった。

## 検討した案

計画には、ほかの次数を比べた記録はない。

## 決定と理由

3D の2群 NEM を、4次多項式展開、横方向の漏れは2次近似で実装する。

## 確認

均質な直方体炉心の解析解（2群拡散方程式、B² = Σ(π/L)²）との一致（`bare_cuboid_matches_analytic_buckling`、
`assembly_sized_nodes_converge`）。

## 帰結

CMFD 加速は M3 で入れる（[ADR-0028](0028-accelerate-core-with-nonlinear-cmfd.md)）。
