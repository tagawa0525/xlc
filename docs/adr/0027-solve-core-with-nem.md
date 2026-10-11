---
status: superseded
date: 2026-10-09
supersedes: []
superseded-by: ADR-0102
issues: []
requires: []
---

# ADR-0027: 炉心計算は 3D NEM の2群の拡散計算とし、1集合体を径方向に1ノードとする

## 背景

炉心計算の手法を決める必要があった。

## 検討した案

計画には、ほかの手法を比べた記録はない。

## 決定と理由

3D NEM（ノード展開法）で2群の拡散計算を行う。1集合体を径方向に1ノード（必要なら2×2）とする。商用コードの標準的な
手法で、ADF とピン出力の再構成を含めて確立しているためである。

## 確認

M0 で均質な直方体炉心の解析解と比べた（`bare_cuboid_matches_analytic_buckling`）。M3 で VERA Problem 3〜5 と比べる。

## 結果

展開の次数は [ADR-0045](0045-expand-nem-to-fourth-order.md)、加速は [ADR-0028](0028-accelerate-core-with-nonlinear-cmfd.md)。
